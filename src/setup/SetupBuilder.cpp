// src/setup/SetupBuilder.cpp — see include/simc/setup/SetupBuilder.h.
// Builds the RunSetup (everything the event loop needs that does not change per event) from the input parameters.
// Called by driver/SimcRun.cpp before the pipeline stages are created; the stages only read the result.
//
// Port of the Fortran start-up sequence (simc.f l.53-101 without calculate_central):
// dbase.f:dbase_read post-read logic, init.f:target_init and extreme_trip_thru_target (through the target model's ITargetLossProvider, default material::HallCCryoTarget), theory_init,
// limits_init, radc_init, simc.f:94-101. Line numbers refer to simc_gfortran/{dbase,init,simc}.f.
//
// Numerics: the Fortran is compiled with -fdefault-real-8, so every literal is a double and
// every expression below keeps the Fortran operation order (left to right as written).
// The phases live in SetupReaction.cpp, SetupTheory.cpp (with the data-file readers) and SetupLimits.cpp
// (src/setup/SetupInternal.h); this file keeps the sequence, the dbase.f:689-958 STOP checks, radc_init and simc.f:94-101.
#include "simc/setup/SetupBuilder.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <fstream>
#include <sstream>

#include "simc/core/SimcConstants.h"
#include "simc/core/SpectrometerId.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/core/Units.h"
#include "SetupInternal.h"

namespace simc {

using setup_detail::nint;

// =============================================================================================
// BuildRunSetup
// =============================================================================================
RunSetup BuildRunSetup(const SimcParameters& p, const SetupOptions& opt, ITargetLossProvider* losses) {
    RunSetup s = RunSetup::FromParameters(p, opt.run_name);
    s.built = true;
    s.fortran_compat = opt.fortran_compat;   // D-031 run option (docs/translation/markdown/FORTRAN_DIFFERENCES.md)
    auto stop = [&](const std::string& why) { s.stops.push_back("STOP " + why); };
    auto msg = [&](const std::string& m) { s.messages.push_back(m); };
    // D-032: a retired arm code (SHMS LSA tune) is an input error. The input readers reject it already; this covers
    // parameters built in code. The Fortran STOPs at dbase.f:895-897/915-917.
    for (const std::string& e : ValidateArmCodes(p)) stop(e);
    if (!s.stops.empty()) return s;
    // the phases (each `false` = a Fortran STOP recorded in s.stops: return the setup as it is)
    if (!setup_detail::BuildReactionAndKinematics(s, p)) return s;
    // init.f:target_init (l.1-88), called at dbase.f:551, through the target model's ITargetLossProvider (the pipeline
    // target from SimcRun; default material::HallCCryoTarget, the Fortran's target model). The same provider object
    // serves extreme_trip_thru_target in limits_init, like the single Fortran process.
    material::HallCCryoTarget default_target;
    if (!losses) losses = &default_target;
    double benharEmMax = 0.0;   // Emval(numEm)
    if (!setup_detail::BuildTargetAndTheory(s, p, opt, losses, benharEmMax)) return s;
    if (!setup_detail::BuildLimits(s, p, losses, benharEmMax)) return s;
    const int A = nint(p.targ.A);
    RunSwitches& sw = s.sw;
    TargetInfo& t = s.targ;
    RadSetup& rad = s.rad;

    // =========================================================================================
    // dbase.f:689-958 announcements: only the STOP branches matter for the setup
    // =========================================================================================
    const ReactionFlags& ff = s.flags;
    auto ngammaOk = [&]() { return p.pizero_ngamma == 1 || p.pizero_ngamma == 2; };
    if (ff.doing_eep) {
        if (!(ff.doing_hyd_elast || ff.doing_deuterium || ff.doing_heavy || ff.doing_nuc_elast)) {   // upstream dbase.f:700
            stop("I don't have ANY idea what (e,e'p) we're doing!!!");
            return s;
        }
    } else if (ff.doing_semi) {
        if (ff.doing_semipi) {
            if ((ff.doing_hydsemi || ff.doing_deutsemi || ff.doing_hesemi) && ff.doing_pizero && !ngammaOk()) {
                stop("pizero_ngamma not set correctly(must be 1 or 2), stopping");
                return s;
            }
        } else if (!ff.doing_semika) {
            stop("I don't have ANY idea what (e,e'm)X we're doing!!!");
            return s;
        }
    } else if (ff.doing_delta) {
        // FORTRAN-BUG(D-018, D-048 F1): original Fortran dbase.f:780-788:
        // `if (doing_hyddelta) ... else if (doing_deutpi) ... else if (doing_hepi) ... else stop 'I don''t have ANY idea
        // what (e,e''p)pi we''re doing!!!'` tests the PION flags doing_deutpi / doing_hepi, never set for doing_delta, so
        // every D or A (e,e'p)pi run STOPs; corrected because the flags meant are doing_deutdelta / doing_hedelta (as in
        // dbase.f:187-189 and the D/A delta branches of event.f/init.f). The D/A branches then run; they were never run
        // in the Fortran (no reference). fortran_compat keeps the STOP.
        const bool known = s.fortran_compat ? (ff.doing_hyddelta || ff.doing_deutpi || ff.doing_hepi)
                                            : (ff.doing_hyddelta || ff.doing_deutdelta || ff.doing_hedelta);
        if (!known) {
            stop("I don't have ANY idea what (e,e'p)pi we're doing!!!");
            return s;
        }
    } else if (ff.doing_pion) {
        if (!(ff.doing_hydpi || ff.doing_deutpi || ff.doing_hepi)) {
            stop("I don't have ANY idea what (e,e'pi) we're doing!!!");
            return s;
        }
        if (ff.doing_pizero && !ngammaOk()) {
            stop("pizero_ngamma not set correctly(must be 1 or 2), stopping");
            return s;
        }
        // FORTRAN-BUG(D-038 F2): original Fortran dbase.f:141-142 stops pi- production from hydrogen
        // (`nint(targ%A).eq.1 .and. which_pion.eq.1 .and. .not. doing_pizero`) but accepts pi0 "from the neutron"
        // (which_pion = 1 with doing_pizero) on hydrogen, which then runs with Mtar_struck = Mn on a proton target;
        // corrected because hydrogen has no neutron: the configuration is rejected at start-up (std::invalid_argument,
        // exit status 2). fortran_compat keeps the Fortran (the run proceeds).
        if (!s.fortran_compat && ff.doing_pizero && A == 1 && p.which_pion == 1)
            throw std::invalid_argument("doing_pizero = 1 with which_pion = 1 (pi0 production from a NEUTRON) on a hydrogen "
                                        "target (A = 1): hydrogen has no neutron (the Fortran accepts it and uses "
                                        "Mtar_struck = Mn, dbase.f:141-142, D-038 F2). Use which_pion = 0 (pi0 from the "
                                        "proton).");
        // FORTRAN-BUG(D-038 F1, D-018): original Fortran dbase.f:823-826:
        // `if(doing_hydpi .and. which_pion.eq.3) stop 'Can''t do Delta+ production from neutron for pi0'` does not test
        // doing_pizero, so charged H(e,e'pi-)Delta++ (set up at dbase.f:361-364, weighted 0.55 at event.f:1479-1484;
        // infiles/fpi12_deltapp.inp) can never run; corrected because the message and the physics only exclude the pi0
        // (gamma* n -> pi0 Delta0 needs a neutron): the STOP is limited to doing_pizero. fortran_compat keeps the STOP.
        if (p.which_pion == 3 && ff.doing_hydpi && (s.fortran_compat || ff.doing_pizero)) {
            stop("Can't do Delta+ production from neutron for pi0");
            return s;
        }
    } else if (ff.doing_kaon) {
        if (!(ff.doing_hydkaon || ff.doing_deutkaon || ff.doing_hekaon)) {
            stop("I don't have ANY idea what (e,e'K) we're doing!!!");
            return s;
        }
        const int w = p.which_kaon;
        if (!(w == 0 || w == 1 || w == 2 || w == 10 || w == 11 || w == 12)) {
            stop("I don't have ANY idea what (e,e'K) we're doing!!!");
            return s;
        }
    } else if (ff.doing_rho) {
        if (ff.doing_hydrho) {
            if (p.targ.A >= 3) {
                stop("A(e,e'rho) Not yet implemented");
                return s;
            }
        } else if (!(ff.doing_deutrho || ff.doing_herho)) {
            stop("I don't have ANY idea what (e,e'rho) we're doing!!!");
            return s;
        }
    } else if (!ff.doing_phsp) {
        stop("I don't have ANY idea what we're doing!!!");
        return s;
    }
    // dbase.f:895-897, 915-917 (retired arm code: "This spectrometer option no longer available" + STOP) is checked
    // at the start of BuildRunSetup through ValidateArmCodes (D-032).
    if (ff.doing_phsp) msg("**WARNING: phsp option is not really implemented - buyer beware!!");

    // =========================================================================================
    // init.f:radc_init (l.576-654), called at simc.f:84
    // =========================================================================================
    if (sw.extrad_flag == 0) {
        if (sw.rad_flag == 0) {
            sw.extrad_flag = 3;
        } else if (sw.rad_flag == 1 || sw.rad_flag == 2 || sw.rad_flag == 3) {
            sw.extrad_flag = 1;
        }
    } else if (sw.extrad_flag < 0) {
        stop("Imbecile! check your stupid setting of EXTRAD_FLAG");
        return s;
    }
    rad.etatzai = (12.0 + (t.Z + 1.) / (t.Z * t.L1 + t.L2)) / 9.0;
    sw.produce_output = p.debug[0];
    sw.exponentiate = p.use_expon == 1;
    sw.include_hard = true;
    sw.calculate_spence = true;

    // =========================================================================================
    // simc.f:94-101 target factor and luminosity
    // =========================================================================================
    // FORTRAN-QUIRK: `3.75914d+6` is a DOUBLE PRECISION literal, which -fdefault-real-8 promotes to real*16, so
    // the Fortran evaluates the whole expression in quad precision (objdump of simc.o: __divtf3/__multf3) and rounds
    // the result to real*8. The double-precision evaluation differs by 1 ulp for the D(e,e'p) capture input.
    {
        __extension__ typedef __float128 quad;
        quad q = static_cast<quad>(t.mass_amu) / static_cast<quad>(3.75914e+6);
        q = q / static_cast<quad>(t.abundancy / 100.);
        q = q * static_cast<quad>(std::abs(std::cos(t.angle)));
        q = q / static_cast<quad>(1000. * t.thick);
        s.targetfac = static_cast<double>(q);
    }
    s.luminosity = p.EXPER.charge / s.targetfac;
    return s;
}

}  // namespace simc
