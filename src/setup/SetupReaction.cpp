// src/setup/SetupReaction.cpp — BuildRunSetup phase 1 (see SetupInternal.h): reaction, kinematics and target input,
// dbase.f:124-547 post-read logic: reaction flags and Mh, run switches, kinematics_main, target masses, struck/recoil
// masses and sign_hadron, target geometry, offsets, SPedge, radiation tails, generation angles, flags.
// Numerics: the Fortran is compiled with -fdefault-real-8, so every literal is a double and every expression keeps
// the Fortran operation order (left to right as written).
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
namespace setup_detail {

namespace K = constants;

bool BuildReactionAndKinematics(RunSetup& s, const SimcParameters& p) {
    auto stop = [&](const std::string& why) { s.stops.push_back("STOP " + why); };
    auto msg = [&](const std::string& m) { s.messages.push_back(m); };
    const int A = nint(p.targ.A);
    double mrec_guess = 0.0;   // dbase.f local 'Mrec_guess' (static zero under -fno-automatic)

    // ---- dbase.f:124-222 flags (DeriveReactionFlags) and Mh. Mh is chosen BEFORE doing_phsp clears
    // the reaction flags (dbase.f:227-234), so evaluate the branch on the pre-phsp flags.
    std::vector<std::string> flag_stops;
    s.flags = DeriveReactionFlags(p, &flag_stops);
    if (!flag_stops.empty()) {
        s.stops = flag_stops;
        return false;
    }
    SimcParameters pre = p;
    pre.doing_phsp = false;
    const ReactionFlags f0 = DeriveReactionFlags(pre);
    // doing_rho together with another reaction flag. dbase.f:136-215 takes the first branch of its if-chain
    // (pion, kaon, delta, semi) but leaves doing_rho set, so event.f:420-421 still calls rho_decay on that hadron
    // (and complete_main uses the other cross section): no meaningful run. Refused in both modes (input error).
    // FORTRAN-BUG(D-057 F16): original Fortran dbase.f:209-214, 875-876: no doing_*rho flag for nint(A) >= 4, so the run
    // STOPs ('I don't have ANY idea what (e,e'rho) we're doing!!!', exit status 0) after printing E_bind = -1e10;
    // refused in both modes with a clear message (there is no rho model for A >= 4). The refusal is made after the
    // target-geometry STOPs of dbase.f:461/465 (below), which the Fortran reaches first (shipped rhotest.inp, legacy names:
    // targ%A = 0 and targ%can = 0 -> 'bad targ.can value', reference/P15.F1_report.md §3).
    if (p.doing_rho && !p.doing_phsp && (p.doing_pion || p.doing_kaon || p.doing_delta || p.doing_semi))
        throw std::invalid_argument(
            "doing_rho = 1 together with doing_pion/doing_kaon/doing_delta/doing_semi: the Fortran runs the other "
            "reaction and still decays its hadron as a rho (dbase.f:136-215, event.f:420-421). Set only one reaction flag.");
    if (f0.doing_pion) {                         // dbase.f:138-139
        s.Mh = K::Mpi;
        if (f0.doing_pizero) s.Mh = K::Mpi0;
    } else if (f0.doing_kaon) {                  // dbase.f:163
        s.Mh = K::Mk;
    } else if (f0.doing_delta) {                 // dbase.f:182-185
        s.Mh = K::Mp;
        // FORTRAN-BUG(D-048 F6): original Fortran dbase.f:184-186: `write(6,*) 'WARNING: Delta cross section model",`
        // + continuation `" only set up for proton target!'`: a quote typo makes it ONE string, printed with `",`, the
        // 72 blanks padding the 132-column line, the 7 blanks of the continuation and `"`; corrected because the
        // message is meant as one sentence. fortran_compat prints the Fortran text.
        if (A >= 2)
            msg(s.fortran_compat ? "WARNING: Delta cross section model\"," + std::string(79, ' ') +
                                       "\" only set up for proton target!"
                                 : "WARNING: Delta cross section model only set up for proton target!");
    } else if (f0.doing_semi) {                  // dbase.f:191-205
        if (f0.doing_semipi) {
            s.Mh = K::Mpi;
            if (f0.doing_pizero) s.Mh = K::Mpi0;
        } else if (f0.doing_semika) {
            s.Mh = K::Mk;
        }
        if (f0.doing_hydsemi && p.do_fermi) msg("WARNING: Cannot do Fermi motion for Hydrogen!");
    } else if (f0.doing_rho) {                   // dbase.f:208
        s.Mh = K::Mrho;
    } else {                                     // dbase.f:216
        s.Mh = K::Mp;
    }
    s.Mh2 = s.Mh * s.Mh;                         // dbase.f:224-225
    s.Mh2_final = s.Mh2;
    // FORTRAN-BUG(D-037 F9): original Fortran dbase.f:1147 reads doing_decay without tying it to the reaction, and the
    // decay of a hadron that is neither a pion nor a kaon executes a bare `stop` in project.f:87-91 / transp.f:160-165,
    // 253-258 (exit status 0, no output; with the default ctau = 0 at the first project call). Corrected: such a
    // configuration is rejected here with a clear message (and project/transp never check a non-pi/K particle for
    // decay). fortran_compat keeps the Fortran (simc::FortranStop at the first decay).
    // (the same pi/K test as optics::decaysInFlight and project.f:78-79)
    // doing_rho: the detected hadron is the pi+/pi- of rho_decay.f (Mh2 = Mpi2 after the decay, rho_decay.f:148-151),
    // which is what project.f/transp.f see at the decay test.
    const bool hadronDecays = std::fabs(std::sqrt(s.Mh2) - K::Mpi) < 2 || std::fabs(std::sqrt(s.Mh2) - K::Mk) < 2 ||
                              (f0.doing_rho && !p.doing_phsp);
    // D-062 G1: with doing_pizero the detected particles are the decay photons, which the corrected mode tracks as
    // massless and stable: doing_decay = 1 is accepted and has no effect (the Fortran STOPs at the first photon "decay").
    if (p.doing_decay && !s.fortran_compat && !hadronDecays && !f0.doing_pizero) {
        char b[320];
        std::snprintf(b, sizeof b,
                      "doing_decay = 1, but the hadron (mass %.5f MeV) is neither a pion nor a kaon: only pi and K "
                      "decay in flight (the Fortran stops at the first decay, project.f:87-91). Set doing_decay = 0.",
                      s.Mh);
        throw std::invalid_argument(b);
    }
    const ReactionFlags& f = s.flags;            // post-phsp flags from here on

    // ---- switches that later code needs (convert_to_logical values)
    RunSwitches& sw = s.sw;
    sw.using_Eloss = p.using_Eloss;
    sw.correct_Eloss = p.correct_Eloss;
    sw.using_Coulomb = p.using_Coulomb;
    sw.mc_smear = p.mc_smear;
    sw.using_rad = p.using_rad;
    sw.rad_flag = p.rad_flag;
    sw.extrad_flag = p.extrad_flag;
    // FORTRAN-OPEN(D-031 flag 1, D-009): only rad_flag is zeroed; using_rad stays on in acceptance (phsp) mode.
    if (p.doing_phsp) sw.rad_flag = 0;           // dbase.f:227-228

    // ---- dbase.f:238-263 kinematics_main
    s.Ebeam = p.Ebeam;
    s.dEbeam = p.Ebeam * p.dEbeam / units::kPercent;
    auto arm = [&](const SimcParameters::Arm& a, SpecArm& o) {
        o.P = a.P;
        o.theta = std::abs(a.theta) / units::kDegPerRad;
        o.cos_th = std::cos(o.theta);
        o.sin_th = std::sin(o.theta);
        o.offset.x = a.offset.x;
        o.offset.y = a.offset.y;
        o.offset.z = a.offset.z;
        o.offset.xptar = a.offset.xptar;   // converted below (dbase.f:472-475)
        o.offset.yptar = a.offset.yptar;
    };
    arm(p.spec.e, s.spec.e);
    arm(p.spec.p, s.spec.p);
    // dbase.f:247-263 through the spectrometer table (D-032): beam-right arms (HMS, HRS-R, calo HMS side) get
    // phi = 3pi/2, the others pi/2. (The Fortran lists the retired code 6 with pi/2; it is rejected above.)
    auto armPhi = [&](int code, double& phi, const char* which) {
        const std::optional<SpectrometerId> id = SpectrometerFromCode(code);
        if (id && IsBeamRight(*id)) {
            phi = 3 * K::pi / 2.;
        } else if (id) {
            phi = K::pi / 2.;
        } else {
            stop(std::string("I dont know what phi should be for the ") + which + " arm");
            return false;
        }
        return true;
    };
    if (!armPhi(p.electron_arm, s.spec.e.phi, "electron")) return false;   // dbase.f:247-254
    if (!armPhi(p.hadron_arm, s.spec.p.phi, "hadron")) return false;       // dbase.f:256-263

    // ---- dbase.f:267-285 target masses
    TargetInfo& t = s.targ;
    t.A = p.targ.A;
    t.Z = p.targ.Z;
    t.mass_amu = p.targ.mass_amu;
    t.mrec_amu = p.targ.mrec_amu;
    t.rho = p.targ.rho;
    t.thick = p.targ.thick;
    t.angle = p.targ.angle;
    t.abundancy = p.targ.abundancy;
    t.zoffset = p.targ.zoffset;
    t.xoffset = p.targ.xoffset;
    t.yoffset = p.targ.yoffset;
    t.fr1 = p.targ.fr1;
    t.fr2 = p.targ.fr2;
    t.fr_pattern = p.targ.fr_pattern;
    t.can = p.targ.can;
    t.N = t.A - t.Z;
    t.M = t.mass_amu * K::amu;
    t.Mrec = t.mrec_amu * K::amu;
    if (A == 1) {
        if (std::abs(t.M - K::Mp) > 0.1) msg("WARNING: forcing target mass to Mp!!!");
        if (std::abs(t.Mrec) > 0.1) msg("WARNING: forcing recoil mass to zero!!!");
        t.M = K::Mp;
        t.Mrec = 0.;
    } else if (A == 2) {
        if (std::abs(t.M - K::Md) > 0.1) msg("WARNING: forcing target mass to Md!!!");
        t.M = K::Md;
    } else {
        mrec_guess = t.M - K::Mp;   // approx A-1 mass
        if (std::abs(t.Mrec - mrec_guess) > 100.) {
            msg("targ.Mrec replaced by the expectation M - Mp (dbase.f:283-284)");
            t.Mrec = mrec_guess;
        }
    }

    // ---- dbase.f:290-439 struck / recoil masses and sign_hadron (post-phsp flags:
    // FORTRAN-OPEN(D-031 flag 1, D-009/D-028): with doing_phsp no branch is taken, Mtar_struck = Mrec_struck = 0, which
    // makes Egamma_tot_max one proton mass too low (init.f:411-413). Kept as the Fortran until the user decides).
    if (f.doing_eep) {
        t.Mtar_struck = K::Mp;
        t.Mrec_struck = 0.0;
        s.sign_hadron = 1.0;
    } else if (f.doing_delta) {
        t.Mtar_struck = K::Mp;
        t.Mrec_struck = K::Mpi;
        s.sign_hadron = 1.0;
    } else if (f.doing_semi) {
        t.Mtar_struck = K::Mp;
        t.Mrec_struck = K::Mp;
        s.sign_hadron = p.doing_hplus ? 1.0 : -1.0;
    } else if (f.doing_rho) {
        t.Mtar_struck = K::Mp;
        t.Mrec_struck = K::Mp;
        s.sign_hadron = p.doing_hplus ? 1.0 : -1.0;
    } else if (f.doing_pion) {
        const int w = p.which_pion;
        if (f.doing_pizero) {
            if (w == 0) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::Mp; s.sign_hadron = 1.0; }
            else if (w == 1) { t.Mtar_struck = K::Mn; t.Mrec_struck = K::Mn; s.sign_hadron = -1.0; }
            else if (w == 2) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::MDelta; s.sign_hadron = 1.0; }
            else if (w == 3) { t.Mtar_struck = K::Mn; t.Mrec_struck = K::MDelta; s.sign_hadron = -1.0; }
        } else {
            if (w == 0) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::Mn; s.sign_hadron = 1.0; }
            else if (w == 1) { t.Mtar_struck = K::Mn; t.Mrec_struck = K::Mp; s.sign_hadron = -1.0; }
            else if (w == 2) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::MDelta; s.sign_hadron = 1.0; }
            else if (w == 3) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::MDelta; s.sign_hadron = -1.0; }
            else if (w == 10) {
                t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec;
                mrec_guess = t.M - K::Mp + K::Mn; s.sign_hadron = 1.0;
            } else if (w == 11) {
                t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec;
                mrec_guess = t.M - K::Mn + K::Mp; s.sign_hadron = -1.0;
            } else {
                stop("Bad value for which_pion");
                return false;
            }
        }
        if (w >= 10) {   // coherent pion production (dbase.f:382-391)
            if (std::abs(t.Mrec_struck - mrec_guess) > 100.) {
                t.Mrec_struck = mrec_guess;
                msg("For coherent pion production the recoil mass was replaced by M_A + M_N(struck) - M_N(recoil)");
            }
            t.Mrec = 0.;
        }
    } else if (f.doing_kaon) {
        const int w = p.which_kaon;
        if (w == 0) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::Mlambda; s.sign_hadron = 1.0; }
        else if (w == 1) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::Msigma0; s.sign_hadron = 1.0; }
        else if (w == 2) { t.Mtar_struck = K::Mn; t.Mrec_struck = K::Msigma_minus; s.sign_hadron = 1.0; }
        else if (w == 10) {
            t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec; mrec_guess = t.M - K::Mp + K::Mlambda; s.sign_hadron = 1.0;
        } else if (w == 11) {
            t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec; mrec_guess = t.M - K::Mp + K::Msigma0; s.sign_hadron = 1.0;
        } else if (w == 12) {
            t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec; mrec_guess = t.M - K::Mn + K::Msigma_minus; s.sign_hadron = 1.0;
        } else {
            stop("Bad value for which_kaon");
            return false;
        }
        if (w >= 10) {   // dbase.f:427-436
            if (std::abs(t.Mrec_struck - mrec_guess) > 100.) {
                t.Mrec_struck = mrec_guess;
                msg("For coherent kaon production the recoil mass was replaced by M_A - M_N + M_Y");
            }
            t.Mrec = 0.;
        }
    }

    // ---- dbase.f:443-448 deuterium recoil mass
    if (A == 2) {
        mrec_guess = K::Mp + K::Mn - t.Mtar_struck;
        if (std::abs(t.Mrec - mrec_guess) > 0.1) msg("WARNING: forcing recoil mass to M_Nucleon");
        t.Mrec = mrec_guess;
    }

    // ---- dbase.f:451-468 target geometry
    t.thick = t.thick / units::kMgPerG;
    t.length = t.thick / t.rho;
    t.angle = t.angle / units::kDegPerRad;
    s.targ_Bangle = p.targ_Bangle / units::kDegPerRad;
    s.targ_Bphi = p.targ_Bphi / units::kDegPerRad;
    if (t.Z < 2.4) {   // liquid target
        if (std::abs(t.angle) > 0.0001) {
            t.angle = 0.0;
            msg("Forcing target angle to zero for cryotarget.");
        }
        if (t.can != 1 && t.can != 2 && t.can != 3) {
            stop("bad targ.can value");
            return false;
        }
    }
    if (std::sin(t.angle) > 0.85) {
        stop("BAD targ.angle (0 is perp. to beam, +ve is rotated towards SOS)");
        return false;
    }
    // D-057 F16 (see above): rho with nint(A) outside 1..3, refused after the dbase.f:461/465 STOPs
    if (p.doing_rho && !p.doing_phsp && (A >= 4 || A < 1))
        throw std::invalid_argument("doing_rho = 1 with targ%A = " + std::to_string(p.targ.A) +
                                    ": exclusive rho production is only set up for H, D and 3He (dbase.f:209-214; the "
                                    "Fortran STOPs at dbase.f:876 'I don't have ANY idea what (e,e'rho) we're doing!!!')");

    // ---- dbase.f:472-475 offsets mrad -> rad
    s.spec.e.offset.xptar = s.spec.e.offset.xptar / units::kMradPerRad;
    s.spec.e.offset.yptar = s.spec.e.offset.yptar / units::kMradPerRad;
    s.spec.p.offset.xptar = s.spec.p.offset.xptar / units::kMradPerRad;
    s.spec.p.offset.yptar = s.spec.p.offset.yptar / units::kMradPerRad;

    // ---- dbase.f:479-488 SPedge (min_max_init ±1e10 defaults already in params)
    auto spedge = [](const SimcParameters::ArmEdge& in, ArmCuts& o) {
        o.delta.min = in.delta.min;
        o.delta.max = in.delta.max;
        if (o.delta.min <= -100.0) o.delta.min = -99.99;
        o.yptar.min = in.yptar.min / units::kMradPerRad;
        o.yptar.max = in.yptar.max / units::kMradPerRad;
        o.xptar.min = in.xptar.min / units::kMradPerRad;
        o.xptar.max = in.xptar.max / units::kMradPerRad;
    };
    spedge(p.SPedge.e, s.SPedge.e);
    spedge(p.SPedge.p, s.SPedge.p);

    // ---- dbase.f:494-518 radiation tails
    const int ot = p.one_tail;
    sw.doing_tail[0] = ot == 0 || ot == 1 || ot == -2 || ot == -3;
    sw.doing_tail[1] = ot == 0 || ot == 2 || ot == -3 || ot == -1;
    sw.doing_tail[2] = ot == 0 || ot == 3 || ot == -1 || ot == -2;
    if (!sw.doing_tail[0] || !sw.doing_tail[1]) msg("**WARNING: Shutting off tails 1 or 2 not well defined");
    if (std::abs(ot) > 3 && p.using_rad) {
        stop("Moron! one_tail>3 turns radiation off, but using_rad wants it on.");
        return false;
    }
    if (!p.using_rad) sw.doing_tail = {false, false, false};
    sw.hardwired_rad = p.Egamma_gen_max > 0.01;

    // ---- dbase.f:521-528 gen angles /1000 (only the ±1e10 min_max_init values at this point)
    GenLimits& gen = s.gen;
    gen.xwid = p.gen.xwid;
    gen.ywid = p.gen.ywid;
    gen.e.yptar.min = gen.e.yptar.min / units::kMradPerRad;
    gen.e.yptar.max = gen.e.yptar.max / units::kMradPerRad;
    gen.p.yptar.min = gen.p.yptar.min / units::kMradPerRad;
    gen.p.yptar.max = gen.p.yptar.max / units::kMradPerRad;
    gen.e.xptar.min = gen.e.xptar.min / units::kMradPerRad;
    gen.e.xptar.max = gen.e.xptar.max / units::kMradPerRad;
    gen.p.xptar.min = gen.p.xptar.min / units::kMradPerRad;
    gen.p.xptar.max = gen.p.xptar.max / units::kMradPerRad;

    // ---- dbase.f:531-547 flags
    sw.using_E_arm_montecarlo = p.spect_mode != 1 && p.spect_mode != -1;
    sw.using_P_arm_montecarlo = p.spect_mode != 1 && p.spect_mode != -2;
    // D-065 S2: the Fortran accepts any spect_mode; every value other than 1/-1/-2 silently means "both arms on".
    // Corrected mode says so (log only, the run continues as in the Fortran); fortran_compat prints nothing.
    if (!s.fortran_compat && p.spect_mode != 0 && p.spect_mode != 1 && p.spect_mode != -1 && p.spect_mode != -2)
        msg("WARNING (simc++): spect_mode = " + std::to_string(p.spect_mode) +
            " is not one of 0, 1, -1, -2: running with BOTH arm Monte Carlos ON, as the Fortran does (dbase.f:531-532)");
    s.cuts.Em.min = p.cuts.Em.min;
    s.cuts.Em.max = p.cuts.Em.max;
    if (f.doing_pion || f.doing_kaon || f.doing_delta || (s.cuts.Em.min == s.cuts.Em.max)) {
        s.cuts.Em.min = -1.e6;
        s.cuts.Em.max = 1.e6;
    }
    if (std::abs(p.deForest_flag) > 1) {
        stop("Idiot! check setting of deForest_flag");
        return false;
    }
    if (sw.correct_Eloss && !sw.using_Eloss) {
        msg("using_Eloss is false, SO WE ARE FORCING CORRECT_ELOSS TO BE FALSE!!!");
        sw.correct_Eloss = false;
    }
    if (nint(t.Z) == 1) sw.using_Coulomb = false;   // no coulomb corr. for Z=1
    return true;
}

}  // namespace setup_detail
}  // namespace simc
