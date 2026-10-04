// src/material/HallCCryoTarget.cpp — see include/simc/material/HallCCryoTarget.h.
#include "simc/material/HallCCryoTarget.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/SpectrometerId.h"
#include "simc/material/FortranFormat.h"
#include "simc/material/MaterialConstants.h"
#include "simc/material/SimcMaterial.h"
#include "simc/stages/IRandom.h"
#include "simc/core/Units.h"

namespace simc {
namespace material {

namespace {

/// RNG for calls that must not consume random numbers (typeflags 2-4).
class NoRandom : public IRandom {
public:
    std::string name() const override { return "none"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}
    double uniform() override { throw std::logic_error("hallc_cryo: unexpected random-number request"); }
    double gauss1(double) override { throw std::logic_error("hallc_cryo: unexpected random-number request"); }
    void seed(std::uint64_t) override {}
};

}  // namespace

HallCCryoTarget::HallCCryoTarget() : material_(std::make_unique<SimcMaterial>()), log_(&std::cout) {}

void HallCCryoTarget::setLog(std::ostream* os) {
    log_ = os;
    if (auto* sm = dynamic_cast<SimcMaterial*>(material_.get())) sm->setLog(os);
}

// Fortran: dbase.f:451-468 (targ conversions and checks), 239-244 (spectrometer angles),
//          547 (no Coulomb correction for Z=1), 551 (target_init arguments).
void HallCCryoTarget::configure(const Json& config, const SimcParameters& params) {
    if (config.is_object() && config.contains("verbose") && !config.at("verbose").get<bool>()) log_ = nullptr;
    const auto& t = params.targ;
    targ_ = TargetInfo{};
    targ_.A = t.A;
    targ_.Z = t.Z;
    targ_.N = t.A - t.Z;
    targ_.mass_amu = t.mass_amu;
    targ_.mrec_amu = t.mrec_amu;
    targ_.rho = t.rho;
    targ_.xoffset = t.xoffset;
    targ_.yoffset = t.yoffset;
    targ_.fr_pattern = t.fr_pattern;
    targ_.fr1 = t.fr1;
    targ_.fr2 = t.fr2;
    targ_.zoffset = t.zoffset;
    targ_.abundancy = t.abundancy;
    targ_.can = t.can;
    targ_.thick = t.thick / units::kMgPerG;                      // dbase.f:451 mg/cm2 -> g/cm2
    targ_.length = targ_.thick / targ_.rho;              // dbase.f:452
    targ_.angle = t.angle / units::kDegPerRad;           // dbase.f:453
    if (targ_.Z < 2.4) {                                 // liquid target (dbase.f:456-462)
        if (std::fabs(targ_.angle) > 0.0001) {
            targ_.angle = 0.0;
            if (log_) *log_ << " Forcing target angle to zero for cryotarget.\n";
        }
        if (targ_.can != 1 && targ_.can != 2 && targ_.can != 3)
            throw FortranStop("bad targ.can value (dbase.f:461)");
    }
    if (std::sin(targ_.angle) > 0.85) {                  // dbase.f:463-468
        if (log_) {
            *log_ << " BAD targ.angle (0 is perp. to beam, +ve is rotated towards SOS)\n";
            *log_ << " If the target really is >60 degrees from normal to the beam, then\n";
            *log_ << " comment out this error checking code.\n";
        }
        throw FortranStop("BAD targ.angle (dbase.f:463-468)");
    }

    setArms(params.electron_arm, params.hadron_arm);

    double Mh = constants::Mp;
    if (config.is_object() && config.contains("Mh")) Mh = config.at("Mh").get<double>();
    init_.using_Eloss = params.using_Eloss;
    init_.using_Coulomb = params.using_Coulomb;
    if (std::lround(targ_.Z) == 1) init_.using_Coulomb = false;   // dbase.f:547
    init_.the_cent = std::fabs(params.spec.e.theta) / units::kDegPerRad;   // dbase.f:239
    init_.thp_cent = std::fabs(params.spec.p.theta) / units::kDegPerRad;   // dbase.f:242
    init_.Pp_cent = params.spec.p.P;
    init_.Pe_cent = params.spec.e.P;
    init_.Ebeam = params.Ebeam;
    init_.Mh2 = Mh * Mh;

    std::string mat = "simc_material";
    if (config.is_object() && config.contains("material")) mat = config.at("material").get<std::string>();
    material_ = Registry<IMaterialModel>::create(mat);
    material_->configure(config.is_object() && config.contains("material_config") ? config.at("material_config")
                                                                                  : Json::object(),
                         params);
    setLog(log_);
}

void HallCCryoTarget::init(const StageContext& ctx) {
    material_->init(StageContext(ctx.setup));
    targetInit(init_.using_Eloss, init_.using_Coulomb, init_.the_cent, init_.thp_cent, init_.Pp_cent, init_.Mh2,
               init_.Ebeam, init_.Pe_cent);
}

// Fortran: event.f:185  main%target%z = (0.5-grnd())*targ%length+targ%zoffset
double HallCCryoTarget::generateZ(IRandom& rng) { return (0.5 - rng.uniform()) * targ_.length + targ_.zoffset; }

ITarget::Traversal HallCCryoTarget::tripThruTarget(Arm arm, double zpos, double energy, double theta, double mass,
                                                   LossMode mode, IRandom& rng) {
    Traversal t;
    trip(static_cast<int>(arm), zpos, energy, theta, t.eloss, t.radlen, mass, static_cast<int>(mode), rng);
    return t;
}

ITarget::ScatteringAngles HallCCryoTarget::multipleScattering(double p, double beta, double teff, IRandom& rng) {
    double dangles[2];
    targetMusc(p, beta, teff, dangles, rng);
    return {dangles[0], dangles[1]};
}

// Fortran: target.f:trip_thru_target (l.1-306)
void HallCCryoTarget::trip(int narm, double zpos, double energy, double theta, double& Eloss, double& radlen,
                           double mass, int typeflag, IRandom& rng) {
    IMaterialModel& mat = *material_;
    const TargetInfo& targ = targ_;
    double forward_path = 0, side_path;
    double s_target, s_Al, s_kevlar = 0, s_air = 0, s_mylar = 0;
    double Eloss_target, Eloss_Al, Eloss_air, Eloss_kevlar, Eloss_mylar;

    s_Al = 0.0;
    const bool liquid = targ.Z < 2.4;

    if (std::fabs(zpos) > (targ.length / 2. + 1.e-5)) {
        if (log_) {
            *log_ << " call to trip_thru_target has |zpos| > targ.length/2.\n";
            *log_ << " could be numerical error, or could be error in target offset\n";
            *log_ << " zpos= " << fortran::listReal(zpos) << "   targ%length/2.= " << fortran::listReal(targ.length / 2.)
                  << "\n";
        }
    }

    // goto (10,20,30) narm — FORTRAN-QUIRK: out-of-range narm falls through to label 10 (beam).
    if (narm != 2 && narm != 3) {
        // 10: the incoming electron
        s_target = (targ.length / 2. + zpos) / std::fabs(std::cos(targ.angle));
        if (liquid) {
            if (targ.can == 1) {            // beer can (2.8 mil endcap)
                s_Al = s_Al + 0.0028 * inch_cm;
            } else if (targ.can == 2) {     // pudding can (5 mil Al, for now)
                s_Al = s_Al + 0.0050 * inch_cm;
            } else if (targ.can == 3) {     // 2017 target 10 cm cells
                s_Al = s_Al + 0.013;        // avg. of 3 loops
            }
        }
        // ... compute distance in radiation lengths and energy loss
        radlen = s_target / targ.X0_cm + s_Al / X0_cm_Al;
        Eloss_target = mat.energyLoss(s_target, targ.rho, targ.Z, targ.A, energy, mass, typeflag, rng);
        Eloss_Al = mat.energyLoss(s_Al, rho_Al, Z_Al, A_Al, energy, mass, typeflag, rng);
        Eloss = Eloss_target + Eloss_Al;
        return;
    }

    // 20: the scattered electron / 30: the scattered hadron (identical code, different arm and warning)
    // Window materials per spectrometer (target.f:80-131 / 200-247), through the spectrometer table (D-032). The
    // Fortran also lists the retired code 6 with the SHMS (target.f:97/217); the input rejects 6.
    const int arm = (narm == 2) ? electron_arm_ : hadron_arm_;
    const std::optional<SpectrometerId> spec = SpectrometerFromCode(arm);
    if (!spec || *spec == SpectrometerId::CALO_HMS_SIDE || *spec == SpectrometerId::CALO_SHMS_SIDE) {
        // FORTRAN-QUIRK(D-012 F4): the calorimeter arms (7/8) and other codes leave s_air, s_kevlar, s_mylar
        // and forward_path at the previous call's values. Out of scope (P14).
        throw std::invalid_argument("hallc_cryo trip_thru_target: " + std::string(narm == 2 ? "electron" : "hadron") +
                                    "_arm = " + std::to_string(arm) +
                                    " not supported (HMS, SOS, HRS-R, HRS-L, SHMS only; calorimeter arms: P14)");
    }
    switch (*spec) {
        case SpectrometerId::HMS:
            s_Al = 0.016 * inch_cm;
            s_air = 15;
            s_kevlar = 0.015 * inch_cm;
            s_mylar = 0.005 * inch_cm;
            forward_path = (targ.length / 2. - zpos) / std::fabs(std::cos(theta + targ.angle));
            break;
        case SpectrometerId::SOS:
            s_Al = 0.008 * inch_cm;
            s_air = 15;
            s_kevlar = 0.005 * inch_cm;
            s_mylar = 0.003 * inch_cm;
            forward_path = (targ.length / 2. - zpos) / std::fabs(std::cos(theta - targ.angle));
            break;
        case SpectrometerId::HRSR:
        case SpectrometerId::HRSL:
            s_Al = 0.013 * inch_cm;
            s_air = 15;
            s_kevlar = 0. * inch_cm;
            s_mylar = 0.010 * inch_cm;
            forward_path = (targ.length / 2. - zpos) / std::fabs(std::cos(theta - targ.angle));
            break;
        case SpectrometerId::SHMS:
            s_Al = (0.02 + 0.01) * inch_cm;   // 20mil scattering chamber window + 10 mil entrance window
            s_air = 57.27;
            s_kevlar = 0.0;
            s_mylar = 0.0;
            forward_path = (targ.length / 2. - zpos) / std::fabs(std::cos(theta - targ.angle));
            break;
        default:
            break;   // calorimeters: refused above
    }
    s_target = forward_path;

    if (liquid) {
        if (targ.can == 1) {                // beer can
            side_path = 1.325 * inch_cm / std::fabs(std::sin(theta));
            if (forward_path < side_path) {
                s_Al = s_Al + 0.005 * inch_cm / std::fabs(std::cos(theta));
            } else {
                s_target = side_path;
                s_Al = s_Al + 0.005 * inch_cm / std::fabs(std::sin(theta));
            }
        } else if (targ.can == 2) {         // pudding can (5 mil Al, for now)
            const double tt = std::tan(theta);
            const double t = tt * tt;
            const double atmp = 1 + t;
            const double btmp = -2 * zpos * t;
            const double half = targ.length / 2.;
            const double ctmp = zpos * zpos * t - half * half;
            const double z_can = (-btmp + std::sqrt(btmp * btmp - 4. * atmp * ctmp)) / 2. / atmp;
            side_path = (z_can - zpos) / std::fabs(std::cos(theta));
            s_target = side_path;
            const double costmp = z_can / (targ.length / 2.);
            if (std::fabs(costmp) <= 1) {
                thCan_ = std::acos(z_can / (targ.length / 2.));
            } else if (std::fabs(costmp - 1.) <= 0.000001) {
                thCan_ = 0.;   // extreme_trip_thru_target can give z/R SLIGHTLY>1.0
            } else {
                // FORTRAN-QUIRK(D-012 F4): th_can keeps its previous value.
                if (log_) {
                    if (narm == 2)
                        *log_ << " z_can > can radius in target.f !!! " << fortran::listReal(z_can) << "\n";
                    else
                        *log_ << " z_can > can radius in target.f !!! " << fortran::listReal(t) << " "
                              << fortran::listReal(theta) << "\n";
                }
            }
            s_Al = s_Al + 0.0050 * inch_cm / std::fabs(std::sin(target_pi / 2 - (theta - thCan_)));
        } else if (targ.can == 3) {         // 2017 10 cm cryo cells
            const double ecir = 1.315 * 2.54;              // endcap inner radius (cm)
            const double ecor = (1.315 + 0.0071) * 2.54;   // endcap outer radius (cm)
            const double entec = targ.length - ecir;       // entrance to end cap (cm)
            const double twall = ecor - ecir;
            const double tcm = zpos + targ.length / 2.0;   // length of target already traversed in cm
            double tliquid, tal;
            if ((tcm + ecir / std::tan(theta)) < entec) {  // e goes through sidewall
                tliquid = ecir / std::sin(theta);
                tal = twall / std::sin(theta);
            } else {                                       // e goes through end cap
                const double a = (targ.length - ecir - tcm) * std::sin(theta);
                tliquid = (std::sqrt(ecir * ecir - a * a) + (targ.length - ecir - tcm) * std::cos(theta));
                tal = +(std::sqrt(ecor * ecor - a * a) - std::sqrt(ecir * ecir - a * a)) * twall / (ecor - ecir);
            }
            s_Al = s_Al + tal;
            s_target = tliquid;
        }
    }

    // ... compute distance in radiation lengths and energy loss
    radlen = s_target / targ.X0_cm + s_Al / X0_cm_Al + s_air / X0_cm_air + s_kevlar / X0_cm_kevlar +
             s_mylar / X0_cm_mylar;
    Eloss_target = mat.energyLoss(s_target, targ.rho, targ.Z, targ.A, energy, mass, typeflag, rng);
    Eloss_Al = mat.energyLoss(s_Al, rho_Al, Z_Al, A_Al, energy, mass, typeflag, rng);
    Eloss_air = mat.energyLoss(s_air, rho_air, Z_air, A_air, energy, mass, typeflag, rng);
    Eloss_kevlar = mat.energyLoss(s_kevlar, rho_kevlar, Z_kevlar, A_kevlar, energy, mass, typeflag, rng);
    Eloss_mylar = mat.energyLoss(s_mylar, rho_mylar, Z_mylar, A_mylar, energy, mass, typeflag, rng);
    Eloss = Eloss_target + Eloss_Al + Eloss_air + Eloss_kevlar + Eloss_mylar;
}

// Fortran: target.f:extreme_trip_thru_target (l.310-544)
void HallCCryoTarget::extremeTripThruTarget(double ebeam, const Limits& the, const Limits& thp, const Limits& pe,
                                            const Limits& pp, const Limits& z, double m) {
    NoRandom none;
    TargetInfo& targ = targ_;
    constexpr double zero = 0.0e0;
    constexpr double Me = constants::Me;
    double th_corner, th_corner_min, th_corner_max;
    double E1, E2, E3, E4, t1, t2, t3, t4;
    double zz = 0, th1 = 0, th2 = 0;
    double energymin, energymax;
    Limits betap;

    const bool liquid = targ.Z < 2.4;

    // Incoming electron
    trip(1, z.max, ebeam, zero, targ.Eloss[0].max, targ.teff[0].max, Me, 3, none);
    trip(1, z.min, ebeam, zero, targ.Eloss[0].min, targ.teff[0].min, Me, 2, none);

    // Scattered electron
    energymin = pe.min;
    energymax = pe.max;
    if (!liquid) {   // ... solid target is infinitely wide
        trip(2, z.min, energymax, the.max, targ.Eloss[1].max, targ.teff[1].max, Me, 3, none);
        trip(2, z.max, energymin, the.min, targ.Eloss[1].min, targ.teff[1].min, Me, 2, none);
    } else if (targ.can == 1) {   // beer can
        if (z.max >= targ.length / 2.) {
            th_corner_max = target_pi / 2.;
        } else {
            th_corner_max = std::atan(1.25 * inch_cm / (targ.length / 2. - z.max));
        }
        th_corner_min = std::atan(1.25 * inch_cm / (targ.length / 2. - z.min));

        // ... max loss: do we have access to a corner shot?
        if (th_corner_min <= the.max && th_corner_max >= the.min) {
            th_corner = std::max(th_corner_min, the.min);
            zz = targ.length / 2. - 1.25 * inch_cm / std::tan(th_corner);
            th1 = th_corner - .0001;
            th2 = th_corner + .0001;
        } else {
            zz = z.min;
            th1 = the.min;
            th2 = the.max;
        }
        trip(2, zz, energymax, th1, E1, t1, Me, 3, none);
        trip(2, zz, energymax, th2, E2, t2, Me, 3, none);
        targ.Eloss[1].max = std::max(E1, E2);
        targ.teff[1].max = std::max(t1, t2);

        // ........ min loss: try all possibilities
        targ.Eloss[1].min = 1.e10;
        for (int i = 0; i <= 3; ++i) {
            trip(2, z.min + (i / 2) * (z.max - z.min), energymin, the.min + (i % 2) * (the.max - the.min), E1, t1, Me,
                 2, none);
            if (E1 < targ.Eloss[1].min) {
                targ.Eloss[1].min = E1;
                targ.teff[1].min = t1;
            }
        }
    } else if (targ.can == 2 || targ.can == 3) {   // pudding can or cryo 17
        zz = -(targ.length / 2.) / std::tan(the.min);
        zz = std::max(zz, (-targ.length / 2.));
        trip(2, zz, energymax, the.min, targ.Eloss[1].max, targ.teff[1].max, Me, 3, none);

        targ.Eloss[1].min = 1.e10;
        for (int i = 0; i <= 3; ++i) {
            trip(2, z.min + (i / 2) * (z.max - z.min), energymin, the.min + (i % 2) * (the.max - the.min), E1, t1, Me,
                 2, none);
            if (E1 < targ.Eloss[1].min) {
                targ.Eloss[1].min = E1;
                targ.teff[1].min = t1;
            }
        }
    }

    // Scattered proton
    energymin = std::sqrt(pp.min * pp.min + m * m);
    energymax = std::sqrt(pp.max * pp.max + m * m);
    betap.min = pp.min / std::sqrt(pp.min * pp.min + m * m);
    betap.max = pp.max / std::sqrt(pp.max * pp.max + m * m);
    if (!liquid) {   // ... solid target is infinitely wide
        trip(3, z.min, energymin, thp.max, E1, t1, m, 3, none);
        trip(3, z.min, energymax, thp.max, E2, t2, m, 3, none);
        targ.Eloss[2].max = std::max(E1, E2);
        targ.teff[2].max = std::max(t1, t2);

        trip(3, z.max, energymin, thp.min, E1, t1, m, 2, none);
        trip(3, z.max, energymax, thp.min, E2, t2, m, 2, none);
        targ.Eloss[2].min = std::min(E1, E2);
        targ.teff[2].min = std::min(t1, t2);
    } else if (targ.can == 1) {   // beer can
        if (z.max >= targ.length / 2.) {
            th_corner_max = target_pi / 2.;
        } else {
            th_corner_max = std::atan(1.25 * inch_cm / (targ.length / 2. - z.max));
        }
        th_corner_min = std::atan(1.25 * inch_cm / (targ.length / 2. - z.min));
        if (th_corner_min <= thp.max && th_corner_max >= thp.min) {
            th_corner = std::max(th_corner_min, thp.min);
            zz = targ.length / 2. - 1.25 * inch_cm / std::tan(th_corner);
            th1 = th_corner - .0001;
            th2 = th_corner + .0001;
        } else {
            zz = z.min;
            th1 = thp.min;
            th2 = thp.max;
        }
        trip(3, zz, energymin, th1, E1, t1, m, 3, none);
        trip(3, zz, energymin, th2, E2, t2, m, 3, none);
        trip(3, zz, energymax, th1, E3, t3, m, 3, none);
        trip(3, zz, energymax, th2, E4, t4, m, 3, none);
        targ.Eloss[2].max = std::max(std::max(std::max(E1, E2), E3), E4);
        targ.teff[2].max = std::max(std::max(std::max(t1, t2), t3), t4);

        targ.Eloss[2].min = 1.e10;
        for (int i = 0; i <= 3; ++i) {
            trip(3, z.min + (i / 2) * (z.max - z.min), energymin, thp.min + (i % 2) * (thp.max - thp.min), E1, t1, m,
                 2, none);
            if (E1 < targ.Eloss[2].min) {
                targ.Eloss[2].min = E1;
                targ.teff[2].min = t1;
                zz = z.min + (i / 2) * (z.max - z.min);
                th1 = thp.min + (i % 2) * (thp.max - thp.min);
            }
        }
        trip(3, zz, energymax, th1, E1, t1, m, 2, none);
        targ.Eloss[2].min = std::min(targ.Eloss[2].min, E1);
    } else if (targ.can == 2 || targ.can == 3) {   // pudding can or cryo17
        // FORTRAN-QUIRK (suspected bug, target.f:501): the ELECTRON angle the%min is used here.
        zz = -(targ.length / 2.) / std::tan(the.min);
        zz = std::max(zz, (-targ.length / 2.));

        trip(3, zz, energymin, thp.min, E1, t1, m, 3, none);
        trip(3, zz, energymax, thp.min, E2, t2, m, 3, none);
        targ.Eloss[2].max = std::max(E1, E2);
        targ.teff[2].max = std::max(t1, t2);

        targ.Eloss[2].min = 1.e10;
        for (int i = 0; i <= 3; ++i) {
            trip(3, z.min + (i / 2) * (z.max - z.min), energymin, thp.min + (i % 2) * (thp.max - thp.min), E1, t1, m,
                 2, none);
            if (E1 < targ.Eloss[2].min) {
                targ.Eloss[2].min = E1;
                targ.teff[2].min = t1;
                zz = z.min + (i / 2) * (z.max - z.min);
                th1 = thp.min + (i % 2) * (thp.max - thp.min);
            }
        }
        trip(3, zz, energymax, th1, E1, t1, m, 2, none);
        targ.Eloss[2].min = std::min(targ.Eloss[2].min, E1);
    }

    // Extreme multiple scattering. Use nominal beam energy rather than minimum.
    extremeTargetMusc(ebeam, 1.e0, targ.teff[0].max, targ.musc_max[0], targ.musc_nsig_max);
    extremeTargetMusc(pe.min, 1.e0, targ.teff[1].max, targ.musc_max[1], targ.musc_nsig_max);
    extremeTargetMusc(pp.min, betap.min, targ.teff[2].max, targ.musc_max[2], targ.musc_nsig_max);
}

// Fortran: target.f:target_musc (l.548-577)
void HallCCryoTarget::targetMusc(double p, double beta, double teff, double dangles[2], IRandom& rng) {
    constexpr double Es = 13.6;   // MeV
    constexpr double epsilon = 0.088;
    if (p < 25. && log_)
        *log_ << " Momentum passed to target_musc should be in MeV, but p= " << fortran::listReal(p) << "\n";
    // Lynch and Dahl, NIM B58 (1991) p.6-10, Eqn. 6
    const double theta_sigma = Es / p / beta * std::sqrt(teff) * (1 + epsilon * std::log10(teff / (beta * beta)));
    dangles[0] = theta_sigma * rng.gauss1(kTargetMuscNsig);
    dangles[1] = theta_sigma * rng.gauss1(kTargetMuscNsig);
}

// Fortran: target.f:extreme_target_musc (entry, l.581-588)
void HallCCryoTarget::extremeTargetMusc(double p, double beta, double teff, double& dangle, double& r) const {
    constexpr double Es = 13.6;
    constexpr double epsilon = 0.088;
    const double theta_sigma = Es / p / beta * std::sqrt(teff) * (1 + epsilon * std::log10(teff / (beta * beta)));
    dangle = theta_sigma * kTargetMuscNsig;
    r = kTargetMuscNsig;
}

// Fortran: init.f:target_init (l.1-88)
void HallCCryoTarget::targetInit(bool using_Eloss, bool using_Coulomb, double the_cent, double thp_cent,
                                 double Pp_cent, double Mh2, double ebeam, double Pe_cent) {
    NoRandom none;
    TargetInfo& targ = targ_;
    constexpr double zero = 0.0e0;
    constexpr double alpha = constants::alpha;

    // The radiation length of the target
    targ.L1 = std::log(184.15) - std::log(targ.Z) / 3.0;
    targ.L2 = std::log(1194.) - 2. * std::log(targ.Z) / 3.0;
    if (targ.Z == 1) {
        targ.L1 = 5.31;
        targ.L2 = 6.144;
    }
    const double za2 = (targ.Z * alpha) * (targ.Z * alpha);
    const double fc = za2 * (1.202 + za2 * (-1.0369 + za2 * 1.008 / (za2 + 1)));

    // PDB values (1999) for H, D, 4He; Tsai formula otherwise
    const long nA = std::lround(targ.A);
    if (nA == 1) {
        targ.X0 = 61.28;
    } else if (nA == 2) {
        targ.X0 = 122.4;
    } else if (nA == 4) {
        targ.X0 = 94.32;
    } else if (nA == 3) {
        if (log_) *log_ << " Using Tsai formula for He-3 radiation length!!! (is there a better value?)\n";
        targ.X0 = 716.405 * targ.A / targ.Z / (targ.Z * (targ.L1 - fc) + targ.L2);
    } else {
        targ.X0 = 716.405 * targ.A / targ.Z / (targ.Z * (targ.L1 - fc) + targ.L2);
    }
    targ.X0_cm = targ.X0 / targ.rho;

    // 'Average' ionization losses (MOST PROBABLE, REALLY)
    double energy = ebeam;
    trip(1, zero, energy, zero, targ.Eloss[0].ave, targ.teff[0].ave, constants::Me, 4, none);
    energy = Pe_cent;
    trip(2, zero, energy, the_cent, targ.Eloss[1].ave, targ.teff[1].ave, constants::Me, 4, none);
    energy = std::sqrt(Pp_cent * Pp_cent + Mh2);
    trip(3, zero, energy, thp_cent, targ.Eloss[2].ave, targ.teff[2].ave, std::sqrt(Mh2), 4, none);
    if (!using_Eloss) {
        for (int i = 0; i < 3; ++i) targ.Eloss[i].ave = 0.0;
    }

    // Coulomb potential energy (Aste et al., Eur. Phys. J. A 26, 167 (2005), 5/15/06 version)
    if (using_Coulomb) {
        targ.Coulomb.ave = 0.75 * 1.5 * (targ.Z - 1.) * alpha * constants::hbarc /
                           (1.1 * std::pow(targ.A, 1. / 3.) + 0.86 * std::pow(targ.A, -1. / 3.));
        targ.Coulomb_constant = targ.Coulomb.ave;
        targ.Coulomb.min = targ.Coulomb_constant;
        targ.Coulomb.max = targ.Coulomb_constant;
    } else {
        targ.Coulomb.ave = 0.0;
        targ.Coulomb_constant = 0.0;
        targ.Coulomb.min = 0.0;
        targ.Coulomb.max = 0.0;
    }
}

// ---- ITargetLossProvider: the set-up's view (task P7.4; was setup/HallCTargetLosses, P2.2) ------------------------
namespace {
void splitLines(const std::string& text, std::vector<std::string>& out) {
    std::istringstream is(text);
    std::string l;
    while (std::getline(is, l)) {
        const auto a = l.find_first_not_of(' ');
        if (a != std::string::npos) out.push_back(l.substr(a));
    }
}
}  // namespace

// Fortran: init.f:target_init (l.1-88), called at dbase.f:551, on the set-up's target description. The messages that
// target_init prints are returned in `messages`. Overwrites this object's target state (SimcRun uses a dedicated
// instance for the set-up).
void HallCCryoTarget::targetInit(const TargetInitArgs& a, simc::TargetInfo& t, std::vector<std::string>& messages) {
    std::ostream* const saved = log_;
    std::ostringstream log;
    setArms(a.electron_arm, a.hadron_arm);
    setLog(&log);
    targ_ = t;   // the set-up target description (same type since P7.4)
    targetInit(a.using_Eloss, a.using_Coulomb, a.the_cent, a.thp_cent, a.Pp_cent, a.Mh2, a.Ebeam, a.Pe_cent);
    t.L1 = targ_.L1;
    t.L2 = targ_.L2;
    t.X0 = targ_.X0;
    t.X0_cm = targ_.X0_cm;
    for (int i = 0; i < 3; ++i) {
        t.Eloss[i].ave = targ_.Eloss[i].ave;
        t.teff[i].ave = targ_.teff[i].ave;
    }
    t.Coulomb = {targ_.Coulomb.ave, targ_.Coulomb.min, targ_.Coulomb.max};
    t.Coulomb_constant = targ_.Coulomb_constant;
    setLog(saved);
    splitLines(log.str(), messages);
}

// Fortran: target.f:extreme_trip_thru_target (+ extreme_target_musc), called at init.f:245. Silent (as before P7.4).
void HallCCryoTarget::extremeLosses(const ExtremeTripArgs& x, simc::TargetInfo& t) {
    std::ostream* const saved = log_;
    setArms(x.electron_arm, x.hadron_arm);
    setLog(nullptr);
    targ_ = t;   // the set-up target description (same type since P7.4)
    extremeTripThruTarget(x.Ebeam, Limits{x.the.min, x.the.max}, Limits{x.thp.min, x.thp.max},
                          Limits{x.pe.min, x.pe.max}, Limits{x.pp.min, x.pp.max}, Limits{x.z.min, x.z.max}, x.m);
    for (int i = 0; i < 3; ++i) {
        t.Eloss[i].min = targ_.Eloss[i].min;
        t.Eloss[i].max = targ_.Eloss[i].max;
        t.teff[i].min = targ_.teff[i].min;
        t.teff[i].max = targ_.teff[i].max;
        t.musc_max[i] = targ_.musc_max[i];
    }
    t.musc_nsig_max = targ_.musc_nsig_max;
    setLog(saved);
}

}  // namespace material
}  // namespace simc

SIMC_REGISTER(simc::ITarget, "hallc_cryo", simc::material::HallCCryoTarget)
