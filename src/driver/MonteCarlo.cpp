// src/driver/MonteCarlo.cpp — see include/simc/driver/MonteCarlo.h.
// Fortran: simc.f:montecarlo (l.1310-1853). Operation order kept for bitwise results.
#include "simc/driver/MonteCarlo.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "simc/core/Json.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SpectrometerId.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/driver/PizeroCalorimeter.h"
#include "simc/driver/TargetFieldAngles.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace driver {

namespace {
constexpr double beta_electron = 1.;   // simc.f:1329
constexpr double zero = 0.0e0;         // simc.f:1347
}  // namespace

// The default plug-in of an arm code, through the spectrometer table (D-032).
MonteCarlo::ArmChoice MonteCarlo::DefaultArm(int code) {
    const std::optional<SpectrometerId> id = SpectrometerFromCode(code);
    if (!id) throw std::invalid_argument("montecarlo: " + ArmCodeError(code, "spectrometer"));
    ArmChoice c;
    switch (*id) {
        case SpectrometerId::HMS:
        case SpectrometerId::SOS:
        case SpectrometerId::HRSR:
        case SpectrometerId::HRSL:
        case SpectrometerId::SHMS:
            c.type = Key(*id);   // plug-ins "hms", "sos" (P14.1), "hrsr", "hrsl" (P14.2), "shms"
            break;
        case SpectrometerId::CALO_HMS_SIDE:
        case SpectrometerId::CALO_SHMS_SIDE:
            c.type = "calo";   // one plug-in for both sides (mc_calo ignores spect_side, D-062 G7; P14.3)
            break;
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
    // Both arms with the same type AND the same config share one object: one Fortran COMMON block of counters (e.g.
    // HMS for both arms). Different configs (e.g. two ideal_box arms with different acceptances) get two objects.
    const std::string key = name + "\n" + cfg.dump();
    auto it = specs_.find(key);
    if (it != specs_.end()) return it->second.get();
    std::unique_ptr<ISpectrometer> s = Registry<ISpectrometer>::create(name);
    s->configure(cfg, setup_.params);
    StageContext ctx(setup_);
    ctx.target = &target_;
    s->init(ctx);
    ISpectrometer* raw = s.get();
    specs_[key] = std::move(s);
    return raw;
}

MonteCarlo::MonteCarlo(const RunSetup& setup, ITarget& target) : MonteCarlo(setup, target, ArmChoice{}, ArmChoice{}) {}

MonteCarlo::MonteCarlo(const RunSetup& setup, ITarget& target, const ArmChoice& hadron, const ArmChoice& electron)
    : setup_(setup), target_(target), log_(&std::cout) {
    // Fortran order of the calls: hadron arm first, then electron arm.
    pArm_ = makeArm(setup.params.hadron_arm, hadron);
    eArm_ = makeArm(setup.params.electron_arm, electron);
    if (setup.params.using_tgt_field || setup.flags.using_tgt_field) initTargetField();
}

// Fortran: simc.f:120-155 (after the random-number initialisation, before the event loop): the arm angles of the
// polarized-target field and trgInit (task P15.2b).
void MonteCarlo::initTargetField() {
    const RunSetup& S = setup_;
    const ArmFieldAngles a = TargetFieldArmAngles(S.targ_Bangle, S.targ_Bphi, S.spec.e.theta, S.spec.e.phi,
                                                  S.spec.p.theta, S.spec.p.phi, S.fortran_compat);
    fieldLog_ = a.log;
    // tgt_field_file: '0' = zero field, ' ' = uniform test field, else a file in the Fortran working directory
    std::string file = S.params.tgt_field_file;
    while (!file.empty() && file.back() == ' ') file.pop_back();
    const bool blank = S.params.tgt_field_file.find_first_not_of(' ') == std::string::npos;
    if (S.params.tgt_field_file.empty())
        throw std::invalid_argument("using_tgt_field = 1 needs tgt_field_file (a field map such as 'trg_field_map.dat', "
                                    "'0' for no field or ' ' for the uniform test field; the Fortran default is unset)");
    std::string path = blank ? std::string(" ") : file;
    if (!blank && file != "0" && !std::filesystem::path(file).is_absolute()) {
        namespace fs = std::filesystem;
        std::vector<fs::path> cands;
        if (!S.data_dir.empty()) cands.push_back(fs::path(S.data_dir) / file);
        cands.push_back(fs::path(file));
        cands.push_back(fs::path(SIMC_FIELD_DATA_DIR) / file);
        path = cands.back().string();
        for (const auto& c : cands)
            if (fs::is_regular_file(c)) {
                path = c.string();
                break;
            }
    }
    field_ = std::make_unique<optics::TargetField>();
    field_->init(path, a.earm * constants::degrad, 0., a.parm * constants::degrad, 0.);
}

// Fortran: simc.f:1573-1590 / 1777-1794: the spectrometer side decides the sign of stheta (right: HMS, HRSR; left: SOS,
// HRSL, SHMS); any other arm STOPs (T17).
double MonteCarlo::fieldStheta(int arm, double stheta) const {
    if (Is(arm, SpectrometerId::HMS) || Is(arm, SpectrometerId::HRSR)) return stheta;
    if (Is(arm, SpectrometerId::SOS) || Is(arm, SpectrometerId::HRSL) || Is(arm, SpectrometerId::SHMS)) return -stheta;
    if (log_) *log_ << " Target field reconstruction not set up for your spectrometer\n";
    throw FortranStop("Target field reconstruction not set up for your spectrometer (simc.f:1583-1585, 1787-1789)");
}

// Fortran: simc.f:montecarlo
bool MonteCarlo::run(const Event& orig, MainEvent& main, Event& recon, IRandom& rng) {
    const RunSetup& S = setup_;
    const SpecArm& sp = S.spec.p;
    const SpecArm& se = S.spec.e;
    const bool mc_smear = S.sw.mc_smear;
    info_ = MonteCarloInfo{};
    // COMMON /gnrl/ Mh, Mh2: the set-up values, or the per-try values of doing_rho (the decay pion, task P13.2)
    const double Mh = main.gnrl_per_event ? main.gnrl.Mh : S.Mh;
    const double Mh2 = main.gnrl_per_event ? main.gnrl.Mh2 : S.Mh2;

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
    // fry is the recon argument of mc_sos/mc_hrs* (SpecInput::raster_fry); HMS/SHMS get xtar_init (SpecInput::fry)

    // Beam multiple scattering
    // dang_in.dyptar: yp offset, goes directly to yp of both arms; dang_in.dxptar: xp offset, mult. by cos_th
    ITarget::ScatteringAngles dang_in, dangles;   // 0 unless mc_smear
    if (mc_smear) dang_in = target_.multipleScattering(orig.Ein, beta_electron, main.target.teff[0], rng);

    //____ P arm ________________
    if (S.sw.using_Eloss) {
        const double d = orig.p.E - main.target.Eloss[2];
        main.SP.p.delta = (std::sqrt(std::fabs(d * d - Mh2)) - sp.P) / sp.P * 100.;
    } else {
        main.SP.p.delta = orig.p.delta;
    }
    // ... multiple scattering
    if (mc_smear) {
        const double beta = orig.p.P / orig.p.E;
        dangles = target_.multipleScattering(orig.p.P, beta, main.target.teff[2], rng);
    } else {
        dangles = ITarget::ScatteringAngles{};
    }
    main.SP.p.yptar = orig.p.yptar + dangles.dyptar + dang_in.dyptar;
    main.SP.p.xptar = orig.p.xptar + dangles.dxptar + dang_in.dxptar * sp.cos_th;

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
        double dx_P_arm = main.SP.p.xptar - sp.offset.xptar;
        double dy_P_arm = main.SP.p.yptar - sp.offset.yptar;
        if (field_) {   // do target field tracking - GAW (simc.f:1425-1432); the returned ok is ignored (T11)
            bool ok_from = true;
            field_->trackFromTgt(x_P_arm, y_P_arm, z_P_arm, dx_P_arm, dy_P_arm,
                                 S.sign_hadron * sp.P * (1 + delta_P_arm / 100.), Mh, 1, ok_from);
        }
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
        in.m2 = Mh2;
        in.ms_flag = mc_smear;
        in.wcs_flag = mc_smear;
        in.decay_flag = S.params.doing_decay;
        in.fry = xtar_init_P;
        in.raster_fry = fry;   // mc_sos/mc_hrs* (simc.f:1472, 1477, 1482)
        in.raster_frx = frx;
        info_.p_called = true;
        const int harm = S.params.hadron_arm;
        const bool caloArm = Is(harm, SpectrometerId::CALO_HMS_SIDE) || Is(harm, SpectrometerId::CALO_SHMS_SIDE);
        if (caloArm && S.flags.doing_pizero) {
            // simc.f:1489-1556: the two pi0 photons through the calorimeter (PizeroCalorimeter, D-062 G1/G3/G4)
            const PizeroCalorimeter pc(harm, S.params.pizero_ngamma, S.fortran_compat);
            const PizeroCaloResult g = pc.detect(*pArm_, in, main.ntup.gamma1, main.ntup.gamma2, rng);
            if (field_) fieldStheta(harm, 0.0);   // simc.f:1583-1585: STOP for a calorimeter arm (T17)
            main.ntup.xcal_gamma1 = g.xcal_gamma1;
            main.ntup.ycal_gamma1 = g.ycal_gamma1;
            main.ntup.xcal_gamma2 = g.xcal_gamma2;
            main.ntup.ycal_gamma2 = g.ycal_gamma2;
            info_.p_stop = g.stop_reason;
            if (!g.ok) {
                info_.fail_arm = 1;
                return false;
            }
            // simc.f:1597-1609: RECON%p%z = y_P_arm, FP%p = the last photon hit, FP%p%path = pathlen; no pi0
            // reconstruction: RECON%p = SP%p
            main.RECON.p.z = y_P_arm;
            main.FP.p.x = g.fp.x_fp;
            main.FP.p.dx = g.fp.dx_fp;
            main.FP.p.y = g.fp.y_fp;
            main.FP.p.dy = g.fp.dy_fp;
            main.FP.p.path = g.pathlen;
            main.RECON.p.delta = main.SP.p.delta;
            main.RECON.p.yptar = main.SP.p.yptar;
            main.RECON.p.xptar = main.SP.p.xptar;
        } else {
            const SpecResult r = pArm_->transport(in, rng);
            // the resmult argument IS ntup%resfac: written only when the routine writes it
            if (r.resolution_factor) main.ntup.resfac = *r.resolution_factor;
            main.ntup.decdist = r.decay_distance;
            main.ntup.Mh2_final = r.Mh2_final;
            info_.p_stop = r.stop_reason;
            const double stheta_field = field_ ? fieldStheta(harm, std::sin(sp.theta)) : 0.0;   // (STOP check first)
            if (!r.ok) {   // (track_to_tgt returns at once for ok = .false.)
                info_.fail_arm = 1;
                return false;
            }
            double delta_rec = r.delta_recon, y_rec = r.y_recon, dx_rec = r.dxdz_recon, dy_rec = r.dydz_recon;
            if (field_) {   // simc.f:1573-1590: reconstruction through the target field
                const double phad = sp.P * (1. + delta_rec / 100.0);
                const double ctheta = std::cos(sp.theta);
                double mom = S.sign_hadron * phad;
                bool ok = true;
                field_->trackToTgt(delta_rec, y_rec, dx_rec, dy_rec, frx, -fry, mom, std::sqrt(r.m2_final), ctheta,
                                   stheta_field, 1, ok, harm,
                                   [this](double f, double& d, double& yp, double& xp, double& yt) {
                                       pArm_->reconstructLastTrack(f, d, yp, xp, yt);
                                   });
                if (!ok) {
                    info_.p_stop = "target_field";
                    info_.fail_arm = 1;
                    return false;
                }
            }
            main.RECON.p.delta = delta_rec;
            main.RECON.p.yptar = dy_rec;
            main.RECON.p.xptar = dx_rec;
            main.RECON.p.z = y_rec;
            main.FP.p.x = r.x_fp;
            main.FP.p.dx = r.dx_fp;
            main.FP.p.y = r.y_fp;
            main.FP.p.dy = r.dy_fp;
            main.FP.p.path = r.pathlen;
            if (S.flags.doing_pizero) {   // simc.f:1605-1609 (no pi0 reconstruction yet)
                main.RECON.p.delta = main.SP.p.delta;
                main.RECON.p.yptar = main.SP.p.yptar;
                main.RECON.p.xptar = main.SP.p.xptar;
            }
        }
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
    recon.p.E = std::sqrt(recon.p.P * recon.p.P + Mh2);
    {
        const double dx_tmp = recon.p.xptar + sp.offset.xptar;
        const double dy_tmp = recon.p.yptar + sp.offset.yptar;
        PhysicsAngles(sp.theta, sp.phi, dx_tmp, dy_tmp, recon.p.theta, recon.p.phi, log_);
    }
    // ... correct for energy loss - use most probable (last flag = 4)
    if (S.sw.correct_Eloss) {
        const double eloss_P_arm =
            target_.tripThruTarget(ITarget::kHadron, zero, recon.p.E, recon.p.theta, Mh, ITarget::kMostProbable, rng)
                .eloss;
        recon.p.E = recon.p.E + eloss_P_arm;
        recon.p.E = std::max(recon.p.E, std::sqrt(Mh2 + 0.000001));
        recon.p.P = std::sqrt(recon.p.E * recon.p.E - Mh2);
    }

    //___ E arm ______________
    main.SP.e.delta = 100 * (orig.e.E - main.target.Eloss[1] - main.target.Coulomb - se.P) / se.P;
    if (mc_smear) {
        dangles = target_.multipleScattering(orig.e.P, beta_electron, main.target.teff[1], rng);
    } else {
        dangles = ITarget::ScatteringAngles{};
    }
    main.SP.e.yptar = orig.e.yptar + dangles.dyptar + dang_in.dyptar;
    main.SP.e.xptar = orig.e.xptar + dangles.dxptar + dang_in.dxptar * se.cos_th;

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
        double dx_E_arm = main.SP.e.xptar - se.offset.xptar;
        double dy_E_arm = main.SP.e.yptar - se.offset.yptar;
        if (field_) {   // GAW - project to z=0 to compare with reconstructed target positions (simc.f:1693-1701)
            bool ok_from = true;
            field_->trackFromTgt(x_E_arm, y_E_arm, z_E_arm, dx_E_arm, dy_E_arm, -se.P * (1 + delta_E_arm / 100.),
                                 constants::Me, -1, ok_from);
        }
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
        if (S.fortran_compat &&
            (Is(S.params.electron_arm, SpectrometerId::CALO_HMS_SIDE) || Is(S.params.electron_arm, SpectrometerId::CALO_SHMS_SIDE)))
            // FORTRAN-BUG(D-062 G2): original Fortran simc.f:1761-1765: the calorimeter electron arm passes `doing_decay`
            // (every other electron arm .false.), so with doing_decay = 1 the electron can "decay" and project.f STOPs;
            // corrected: .false. as for the other electron arms. fortran_compat keeps doing_decay.
            in.decay_flag = S.params.doing_decay;
        in.fry = xtar_init_E;
        in.raster_fry = fry;   // mc_sos/mc_hrs* (simc.f:1738, 1743, 1748)
        in.raster_frx = frx;
        info_.e_called = true;
        const SpecResult r = eArm_->transport(in, rng);
        // FORTRAN-DIFF-CANDIDATE(D-015): simc.f:1334 tmpfact is a SAVEd local; when the arm routine did not reach the hut
        // it keeps the previous event's value, which is still added to ntup%resfac. Question: intended?
        if (r.resolution_factor) tmpfact_ = *r.resolution_factor;
        main.ntup.resfac = main.ntup.resfac + tmpfact_;
        info_.e_stop = r.stop_reason;
        const double stheta_field = field_ ? fieldStheta(S.params.electron_arm, std::sin(se.theta)) : 0.0;
        if (!r.ok) {
            info_.fail_arm = 2;
            return false;
        }
        double delta_rec = r.delta_recon, y_rec = r.y_recon, dx_rec = r.dxdz_recon, dy_rec = r.dydz_recon;
        if (field_) {   // simc.f:1777-1794
            const double pelec = se.P * (1. + delta_rec / 100.0);
            const double ctheta = std::cos(se.theta);
            double mom = -1.0 * pelec;
            bool ok = true;
            field_->trackToTgt(delta_rec, y_rec, dx_rec, dy_rec, frx, -fry, mom, std::sqrt(constants::Me2), ctheta,
                               stheta_field, -1, ok, S.params.electron_arm,
                               [this](double f, double& d, double& yp, double& xp, double& yt) {
                                   eArm_->reconstructLastTrack(f, d, yp, xp, yt);
                               });
            if (!ok) {
                info_.e_stop = "target_field";
                info_.fail_arm = 2;
                return false;
            }
        }
        main.RECON.e.delta = delta_rec;
        main.RECON.e.yptar = dy_rec;
        main.RECON.e.xptar = dx_rec;
        main.RECON.e.z = y_rec;
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
        const double eloss_E_arm = target_
                                       .tripThruTarget(ITarget::kElectron, zero, recon.e.E, recon.e.theta,
                                                       constants::Me, ITarget::kMostProbable, rng)
                                       .eloss;
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
