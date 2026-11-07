// include/simc/spectrometers/calo/CaloSpectrometer.h
//
// ISpectrometer "calo": the calorimeter arm of SIMC (arm codes 7 = HMS side and 8 = SHMS side of the beam; NPS-like),
// task P14.3. It adapts the kernel calo::McCalo (kernel/McCalo.h: calo/mc_calo.f) to the stage interface, as
// simc.f:montecarlo calls mc_calo (l.1489-1563 hadron arm, l.1756-1765 electron arm). One particle per call: the pi0
// two-photon logic (rotation to the calorimeter frame, two calls, ok_gamma, pizero_ngamma) is driver::PizeroCalorimeter
// (driver/PizeroCalorimeter.h), used by the driver for a pi0 hadron arm.
//
// JSON config (all optional): "drift_to_cal" (cm, default params.drift_to_cal: the distance target -> calorimeter face;
// FORTRAN-DIFF-CANDIDATE D-062 G5: 0 by default, unchecked, as dbase.f:1104), "ctau" (cm), "verbose".
// Outputs: the focal plane = the position and slopes at the calorimeter face (on success); the "recon" values are the
// input values unchanged (no reconstruction, mc_calo.f:167); no resolution factor (resmult is not written, D-062 G7).
// Stop counters: trials, slit_hor, slit_vert, successes (COMMON /spectrometer_calo/; .geni 'Calo Trials').
#ifndef SIMC_SPECTROMETERS_CALO_CALOSPECTROMETER_H
#define SIMC_SPECTROMETERS_CALO_CALOSPECTROMETER_H

#include <map>
#include <memory>
#include <string>

#include "simc/stages/ISpectrometer.h"

namespace simc {
namespace calo {

class McCalo;

class CaloSpectrometer : public ISpectrometer {
public:
    CaloSpectrometer();
    ~CaloSpectrometer() override;
    CaloSpectrometer(const CaloSpectrometer&) = delete;
    CaloSpectrometer& operator=(const CaloSpectrometer&) = delete;

    std::string name() const override { return "calo"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    SpecResult transport(const SpecInput& in, IRandom& rng) override;
    const std::map<std::string, long>& stopCounters() const override;

    McCalo& kernel() { return *k_; }
    double driftToCal() const { return drift_; }

private:
    std::unique_ptr<McCalo> k_;
    double drift_{0};
    double ctau_{0};
    int side_{0};
    mutable std::map<std::string, long> cntMap_;
};

}  // namespace calo
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_CALO_CALOSPECTROMETER_H
