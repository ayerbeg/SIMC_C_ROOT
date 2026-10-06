// include/simc/material/SimcMaterial.h
//
// IMaterialModel "simc_material": the SIMC energy-loss and multiple-scattering routines.
//   energyLoss : Fortran enerloss_new.f:enerloss_new (l.1-85)
//   musc       : Fortran shared/musc.f:musc (l.1-58)
//   muscExt    : Fortran shared/musc_ext.f:musc_ext (l.1-54)
// Units as in the Fortran: len cm, dens g/cm3, epart = TOTAL energy MeV, mpart MeV, Eloss MeV;
// m2 MeV^2, p MeV/c, rad_len radiation lengths, angles rad (the musc.f header says mrad: wrong,
// D-012 F5), x/y cm.
// Gaussian truncations exactly as in the Fortran (D-012 F5): enerloss 10 sigma, musc/musc_ext 99 sigma.
//
// Static state (Fortran -fno-automatic / DATA): enerloss_new keeps `x` and `numerr` between calls.
// They are members of this object, so one SimcMaterial instance = the Fortran process state.
//   FORTRAN-QUIRK(D-012 F3): a typeflag outside 1..4 is not rejected; x keeps the value of the
//   previous call (initially 0 -> lambda = 1e5).
//   FORTRAN-QUIRK(D-012 F3): Eloss > KE is clipped to KE - 1e-7 with a warning printed only the
//   first 10 times (numerr is a DATA variable counting over the whole process).
//   FORTRAN-QUIRK(D-012 F2): x == 0 (gauss1 returned 0) gives lambda = 1e5 instead of -2 ln 0.
// Fortran STOP (musc_ext with x_len <= 0 or rad_len < 0) -> throws simc::FortranStop.
#ifndef SIMC_MATERIAL_SIMCMATERIAL_H
#define SIMC_MATERIAL_SIMCMATERIAL_H

#include <iosfwd>
#include <stdexcept>
#include <string>

#include "simc/stages/IMaterialModel.h"
#include "simc/core/Errors.h"

namespace simc {
namespace material {


/// Gaussian truncations (Fortran PARAMETERs / literals).
constexpr double kEnerlossNsig = 10.0e0;   ///< enerloss_new.f: gauss1(10.0e0)
constexpr double kMuscNsig = 99.0e0;       ///< musc.f / musc_ext.f: nsig_max = 99.0e0

class SimcMaterial : public IMaterialModel {
public:
    SimcMaterial();

    std::string name() const override { return "simc_material"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}

    double energyLoss(double len_cm, double density, double zeff, double aeff, double epart, double mpart,
                      int typeflag, IRandom& rng) override;
    void musc(double m2, double p, double rad_len, double& dth, double& dph, IRandom& rng) override;
    void muscExt(double m2, double p, double rad_len, double x_len, double& dph, double& dth, double& y,
                 double& x, IRandom& rng) override;

    /// Where the Fortran write(6,*) warnings go (default std::cout; nullptr = silent).
    void setLog(std::ostream* os) { log_ = os; }
    std::ostream* log() const { return log_; }

    /// Static state of enerloss_new (inspection / reset for tests).
    int numerr() const { return numerr_; }
    double lastX() const { return x_; }
    void resetStaticState() { numerr_ = 0; x_ = 0.0; }
    /// True if the last energyLoss call clamped the loss to the kinetic energy (enerloss_new.f:75-77: the particle
    /// ranged out in the layer). Used by the corrected collimator stepping (D-037 F4).
    bool lastEnergyLossClamped() const { return lastClamped_; }

private:
    std::ostream* log_;
    int numerr_{0};   ///< enerloss_new.f: integer numerr; data numerr /0/
    double x_{0.0};   ///< enerloss_new.f: local x (static under -fno-automatic, zero-initialised)
    bool lastClamped_{false};   ///< diagnostic, not a Fortran variable
};

}  // namespace material
}  // namespace simc

#endif  // SIMC_MATERIAL_SIMCMATERIAL_H
