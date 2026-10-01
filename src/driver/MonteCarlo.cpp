// src/driver/MonteCarlo.cpp — see include/simc/driver/MonteCarlo.h.
// Fortran: simc.f:montecarlo (l.1310-1853). Operation order kept for bitwise results.
#include "simc/driver/MonteCarlo.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "simc/core/Json.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace driver {

namespace {
constexpr double beta_electron = 1.;   // simc.f:1329
constexpr double zero = 0.0e0;         // simc.f:1347
}  // namespace

MonteCarlo::ArmChoice MonteCarlo::DefaultArm(int code) {
    ArmChoice c;
    switch (code) {
        case 1:
            c.type = "hms";
            break;
        case 5:
        case 6:   // 6 = SHMS LSA: mc_shms STOPs on the first call (FortranStop)
            c.type = "shms";
            c.config["spectr"] = code;
            break;
        case 2:
        case 3:
        case 4:
            throw NotPorted("montecarlo: arm code " + std::to_string(code) + " (SOS/HRS) not yet ported (P14)");
        case 7:
        case 8:
            throw NotPorted("montecarlo: arm code " + std::to_string(code) + " (mc_calo) not ported");
        default:
            throw NotPorted("montecarlo: unknown arm code " + std::to_string(code));
    }
    return c;
}

ISpectrometer* MonteCarlo::makeArm(int code, const ArmChoice& choice) {
    std::string name = choice.type;
    Json cfg = Json::object();
    if (name.empty()) {
        const ArmChoice d = DefaultArm(code);
        name = d.type;
        cfg = d.config;
    } else {
        // an explicit type keeps the default config of the arm code when it is the same type (e.g. SHMS "spectr")
        try {
            const ArmChoice d = DefaultArm(code);
            if (d.type == name) cfg = d.config;
        } catch (const NotPorted&) {
        }
    }
    if (choice.config.is_object())
        for (auto it = choice.config.begin(); it != choice.config.end(); ++it) cfg[it.key()] = it.value();
    auto it = specs_.find(name);
    if (it != specs_.end()) return it->second.get();   // both arms the same type: one COMMON block, one object
    std::unique_ptr<ISpectrometer> s = Registry<ISpectrometer>::create(name);
    s->configure(cfg, setup_.params);
    StageContext ctx(setup_);
    ctx.target = &target_;
    s->init(ctx);
    ISpectrometer* raw = s.get();
    specs_[name] = std::move(s);
    return raw;
}

MonteCarlo::MonteCarlo(const RunSetup& setup, ITarget& target) : MonteCarlo(setup, target, ArmChoice{}, ArmChoice{}) {}

MonteCarlo::MonteCarlo(const RunSetup& setup, ITarget& target, const ArmChoice& hadron, const ArmChoice& electron)
    : setup_(setup), target_(target), log_(&std::cout) {
    if (setup.params.using_tgt_field || setup.flags.using_tgt_field)
        throw NotPorted("montecarlo: using_tgt_field (track_from_tgt / track_to_tgt, trg_track.f) not ported");
    // Fortran order of the calls: hadron arm first, then electron arm.
    pArm_ = makeArm(setup.params.hadron_arm, hadron);
    eArm_ = makeArm(setup.params.electron_arm, electron);
}

// Fortran: simc.f:montecarlo
bool MonteCarlo::run(const Event& orig, MainEvent& main, Event& recon, IRandom& rng) {
    const RunSetup& S = setup_;
    const SpecArm& sp = S.spec.p;
    const SpecArm& se = S.spec.e;
    const bool mc_smear = S.sw.mc_smear;
    info_ = MonteCarloInfo{};

    // Prepare the event for the Monte Carlo's and/or spectrometer cuts
    bool success = false;
    main.ntup.resfac = 0.0;   // resfac (see simulate.inc)
    double fry, frx;
    if (S.params.correct_raster) {
        fry = -main.target.rastery;
        frx = -main.target.rasterx;   // frx positive for beam left
    } else {
        fry = 0.0;
        frx = 0.0;
    }
    (void)fry;   // passed to SOS/HRS/calo only (not ported); HMS/SHMS get xtar_init
    (void)frx;

    // Beam multiple scattering
    double dang_in[2], dangles[2];
    if (mc_smear) {
        target_.multipleScattering(orig.Ein, beta_electron, main.target.teff[0], dang_in, rng);
    } else {
        dang_in[0] = 0.0;   // yp offset, goes directly to yp of both arms
        dang_in[1] = 0.0;   // xp offset, mult. by cos_th for xp of both arms
    }

    //____ P arm ________________
    if (S.sw.using_Eloss) {
        const double d = orig.p.E - main.target.Eloss[2];
        main.SP.p.delta = (std::sqrt(std::fabs(d * d - S.Mh2)) - sp.P) / sp.P * 100.;
    } else {
        main.SP.p.delta = orig.p.delta;
    }
    // ... multiple scattering
    if (mc_smear) {
        const double beta = orig.p.P / orig.p.E;
        target_.multipleScattering(orig.p.P, beta, main.target.teff[2], dangles, rng);
    } else {
        dangles[0] = 0.0;
        dangles[1] = 0.0;
    }
    main.SP.p.yptar = orig.p.yptar + dangles[0] + dang_in[0];
    main.SP.p.xptar = orig.p.xptar + dangles[1] + dang_in[1] * sp.cos_th;

    if (S.sw.using_P_arm_montecarlo) {
        // ... change to P arm spectrometer coordinates (TRANSPORT system)
        if (std::fabs(std::cos(sp.phi)) > 0.0001 && log_) {
            *log_ << " y_P_arm, z_P_arm will be incorrect if spec.p.phi <> pi/2 or 3*pi/2\n";
            *log_ << " spec%p%phi= " << sp.phi << " = " << sp.phi * 180 / constants::pi << " degrees\n";
        }
        const double delta_P_arm = main.SP.p.delta;
        double x_P_arm = -main.target.y;
        double y_P_arm = -main.target.x * sp.cos_th - main.target.z * sp.sin_th * std::sin(sp.phi);
        double z_P_arm = main.target.z * sp.cos_th + main.target.x * sp.sin_th * std::sin(sp.phi);
        // %.. Apply spectrometer offset (using spectrometer coordinate system).
        x_P_arm = x_P_arm - sp.offset.x;
        y_P_arm = y_P_arm - sp.offset.y;
        z_P_arm = z_P_arm - sp.offset.z;
        const double dx_P_arm = main.SP.p.xptar - sp.offset.xptar;
        const double dy_P_arm = main.SP.p.yptar - sp.offset.yptar;
        // (using_tgt_field: track_from_tgt — rejected in the constructor)
        // ........ drift this position back to z=0, the plane through the target center
        x_P_arm = x_P_arm - z_P_arm * dx_P_arm;
        y_P_arm = y_P_arm - z_P_arm * dy_P_arm;
        z_P_arm = 0.0;
        const double xtar_init_P = x_P_arm;
        main.SP.p.z = y_P_arm;

        // ... Monte Carlo through the P arm (simc.f:1459-1467: m2 = Mh2, pathlen = 0)
        SpecInput in;
        in.p_spec = sp.P;
        in.th_spec = sp.theta;
        in.delta = delta_P_arm;
        in.x = x_P_arm;
        in.y = y_P_arm;
        in.z = z_P_arm;
        in.dxdz = dx_P_arm;
        in.dydz = dy_P_arm;
        in.m2 = S.Mh2;
        in.ms_flag = mc_smear;
        in.wcs_flag = mc_smear;
        in.decay_flag = S.params.doing_decay;
        in.fry = xtar_init_P;
        info_.p_called = true;
        const SpecResult r = pArm_->transport(in, rng);
        // the resmult argument IS ntup%resfac: written only when the routine writes it
        if (r.resmult_written) main.ntup.resfac = r.resmult;
        main.ntup.decdist = r.decdist;
        main.ntup.Mh2_final = r.Mh2_final;
        info_.p_stop = r.stop_reason;
        if (!r.ok) {
            info_.fail_arm = 1;
            return false;
        }
        main.RECON.p.delta = r.delta_recon;
        main.RECON.p.yptar = r.dydz_recon;
        main.RECON.p.xptar = r.dxdz_recon;
        main.RECON.p.z = r.y_recon;
        main.FP.p.x = r.x_fp;
        main.FP.p.dx = r.dx_fp;
        main.FP.p.y = r.y_fp;
        main.FP.p.dy = r.dy_fp;
        main.FP.p.path = r.pathlen;
        // (doing_pizero: calorimeter only — not reachable here)
    } else {
        // CASE 2: Not using the detailed Monte Carlo, just copy the IN event to the OUT record
        main.RECON.p.delta = main.SP.p.delta;
        main.RECON.p.yptar = main.SP.p.yptar;
        main.RECON.p.xptar = main.SP.p.xptar;
    }

    // Fill recon quantities.
    recon.p.delta = main.RECON.p.delta;
    recon.p.yptar = main.RECON.p.yptar;
    recon.p.xptar = main.RECON.p.xptar;
    recon.p.z = main.RECON.p.z;
    recon.p.P = sp.P * (1. + recon.p.delta / 100.);
    recon.p.E = std::sqrt(recon.p.P * recon.p.P + S.Mh2);
    {
        const double dx_tmp = recon.p.xptar + sp.offset.xptar;
        const double dy_tmp = recon.p.yptar + sp.offset.yptar;
        PhysicsAngles(sp.theta, sp.phi, dx_tmp, dy_tmp, recon.p.theta, recon.p.phi, log_);
    }
    // ... correct for energy loss - use most probable (last flag = 4)
    if (S.sw.correct_Eloss) {
        double eloss_P_arm, r;
        target_.tripThruTarget(ITarget::kHadron, zero, recon.p.E, recon.p.theta, S.Mh, ITarget::kMostProbable, rng,
                               eloss_P_arm, r);
        recon.p.E = recon.p.E + eloss_P_arm;
        recon.p.E = std::max(recon.p.E, std::sqrt(S.Mh2 + 0.000001));
        recon.p.P = std::sqrt(recon.p.E * recon.p.E - S.Mh2);
    }

    //___ E arm ______________
    main.SP.e.delta = 100 * (orig.e.E - main.target.Eloss[1] - main.target.Coulomb - se.P) / se.P;
    if (mc_smear) {
        target_.multipleScattering(orig.e.P, beta_electron, main.target.teff[1], dangles, rng);
    } else {
        dangles[0] = 0.0;
        dangles[1] = 0.0;
    }
    main.SP.e.yptar = orig.e.yptar + dangles[0] + dang_in[0];
    main.SP.e.xptar = orig.e.xptar + dangles[1] + dang_in[1] * se.cos_th;

    if (S.sw.using_E_arm_montecarlo) {
        if (std::fabs(std::cos(se.phi)) > 0.0001 && log_) {
            *log_ << " y_E_arm, z_E_arm will be incorrect if spec.e.phi <> pi/2 or 3*pi/2\n";
            *log_ << " spec.e.phi= " << se.phi << " = " << se.phi * 180 / constants::pi << " degrees\n";
        }
        const double delta_E_arm = main.SP.e.delta;
        double x_E_arm = -main.target.y;
        double y_E_arm = -main.target.x * se.cos_th - main.target.z * se.sin_th * std::sin(se.phi);
        double z_E_arm = main.target.z * se.cos_th + main.target.x * se.sin_th * std::sin(se.phi);
        x_E_arm = x_E_arm - se.offset.x;
        y_E_arm = y_E_arm - se.offset.y;
        z_E_arm = z_E_arm - se.offset.z;
        const double dx_E_arm = main.SP.e.xptar - se.offset.xptar;
        const double dy_E_arm = main.SP.e.yptar - se.offset.yptar;
        x_E_arm = x_E_arm - z_E_arm * dx_E_arm;
        y_E_arm = y_E_arm - z_E_arm * dy_E_arm;
        z_E_arm = 0.0;
        const double xtar_init_E = x_E_arm;
        main.SP.e.z = y_E_arm;

        // simc.f:1726-1755: m2 = me2 (the PARAMETER is passed), decay .false., resmult -> tmpfact
        SpecInput in;
        in.p_spec = se.P;
        in.th_spec = se.theta;
        in.delta = delta_E_arm;
        in.x = x_E_arm;
        in.y = y_E_arm;
        in.z = z_E_arm;
        in.dxdz = dx_E_arm;
        in.dydz = dy_E_arm;
        in.m2 = constants::Me2;
        in.ms_flag = mc_smear;
        in.wcs_flag = mc_smear;
        in.decay_flag = false;
        in.fry = xtar_init_E;
        info_.e_called = true;
        const SpecResult r = eArm_->transport(in, rng);
        // FORTRAN-QUIRK: tmpfact is a SAVEd local, stale if the routine did not write it
        if (r.resmult_written) tmpfact_ = r.resmult;
        main.ntup.resfac = main.ntup.resfac + tmpfact_;
        info_.e_stop = r.stop_reason;
        if (!r.ok) {
            info_.fail_arm = 2;
            return false;
        }
        main.RECON.e.delta = r.delta_recon;
        main.RECON.e.yptar = r.dydz_recon;
        main.RECON.e.xptar = r.dxdz_recon;
        main.RECON.e.z = r.y_recon;
        main.FP.e.x = r.x_fp;
        main.FP.e.dx = r.dx_fp;
        main.FP.e.y = r.y_fp;
        main.FP.e.dy = r.dy_fp;
        main.FP.e.path = r.pathlen;
    } else {
        main.RECON.e.delta = main.SP.e.delta;
        main.RECON.e.yptar = main.SP.e.yptar;
        main.RECON.e.xptar = main.SP.e.xptar;
    }

    // Fill recon quantities
    recon.e.delta = main.RECON.e.delta;
    recon.e.yptar = main.RECON.e.yptar;
    recon.e.xptar = main.RECON.e.xptar;
    recon.e.z = main.RECON.e.z;
    recon.e.P = se.P * (1. + recon.e.delta / 100.);
    recon.e.E = recon.e.P;
    {
        const double dx_tmp = recon.e.xptar + se.offset.xptar;
        const double dy_tmp = recon.e.yptar + se.offset.yptar;
        PhysicsAngles(se.theta, se.phi, dx_tmp, dy_tmp, recon.e.theta, recon.e.phi, log_);
    }
    // ... correct for energy loss (Coulomb deceleration is not corrected, simc.f:1840-1842)
    if (S.sw.correct_Eloss) {
        double eloss_E_arm, r;
        target_.tripThruTarget(ITarget::kElectron, zero, recon.e.E, recon.e.theta, constants::Me,
                               ITarget::kMostProbable, rng, eloss_E_arm, r);
        recon.e.E = recon.e.E + eloss_E_arm;
    }
    // recon%e%E = recon%e%E (simc.f:1842, no-op)
    recon.e.P = recon.e.E;

    // Made it!
    success = true;
    info_.success = true;
    return success;
}

}  // namespace driver
}  // namespace simc
