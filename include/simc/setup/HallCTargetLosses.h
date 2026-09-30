// include/simc/setup/HallCTargetLosses.h
//
// Default ITargetLossProvider (task P2.2): a thin adapter over material::HallCCryoTarget (P3.1), the single
// C++ implementation of init.f:target_init and target.f:extreme_trip_thru_target (+ extreme_target_musc).
// It copies the RunSetup target (simc::TargetInfo) into the HallCCryoTarget state, calls the Fortran-order
// routine and copies the results back. The hadron mass comes from the caller (RunSetup::Mh via the
// TargetInitArgs::Mh2 / ExtremeTripArgs::m arguments), so it is reaction dependent as in the Fortran.
#ifndef SIMC_SETUP_HALLCTARGETLOSSES_H
#define SIMC_SETUP_HALLCTARGETLOSSES_H

#include <memory>

#include "simc/setup/TargetLossProvider.h"

namespace simc {
namespace material {
class HallCCryoTarget;
}

class HallCTargetLosses : public ITargetLossProvider {
public:
    HallCTargetLosses(int electron_arm, int hadron_arm);
    ~HallCTargetLosses() override;

    void targetInit(const TargetInitArgs& args, TargetInfo& targ, std::vector<std::string>& messages) override;
    void extremeLosses(const ExtremeTripArgs& args, TargetInfo& targ) override;

    /// The underlying target (state after the last call).
    material::HallCCryoTarget& target() { return *target_; }

private:
    std::unique_ptr<material::HallCCryoTarget> target_;
};

}  // namespace simc

#endif  // SIMC_SETUP_HALLCTARGETLOSSES_H
