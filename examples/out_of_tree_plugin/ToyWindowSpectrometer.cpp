// examples/out_of_tree_plugin/ToyWindowSpectrometer.cpp — a minimal ISpectrometer written outside the SIMC++ tree
// (task P15.3, docs/EXTENDING.md). "toy_window": accepts a track when |delta| <= delta_max (%) and returns it
// unchanged (perfect reconstruction, unit optics). JSON config: {"type": "toy_window", "delta_max": 8.0}.
// Stop counters (printed in the .geni and the summary): trials, delta, successes.
#include <cmath>
#include <map>
#include <string>

#include "simc/core/Registry.h"
#include "simc/stages/ISpectrometer.h"

namespace mylab {

class ToyWindowSpectrometer : public simc::ISpectrometer {
public:
    std::string name() const override { return "toy_window"; }
    void configure(const simc::Json& config, const simc::SimcParameters&) override {
        if (config.is_object()) deltaMax_ = config.value("delta_max", deltaMax_);
    }
    void init(const simc::StageContext&) override { counts_ = {{"trials", 0}, {"delta", 0}, {"successes", 0}}; }

    simc::SpecResult transport(const simc::SpecInput& in, simc::IRandom&) override {
        simc::SpecResult r;
        r.m2_final = in.m2;
        r.Mh2_final = in.m2;
        ++counts_["trials"];
        if (std::fabs(in.delta) > deltaMax_) {
            ++counts_["delta"];
            r.stop_reason = "delta";
            return r;
        }
        ++counts_["successes"];
        r.ok = true;
        r.delta_recon = in.delta;
        r.dxdz_recon = in.dxdz;
        r.dydz_recon = in.dydz;
        r.y_recon = in.y;
        r.x_fp = in.x;
        r.dx_fp = in.dxdz;
        r.y_fp = in.y;
        r.dy_fp = in.dydz;
        r.resolution_factor = 1.0;
        return r;
    }
    const std::map<std::string, long>& stopCounters() const override { return counts_; }

private:
    double deltaMax_{10.0};
    std::map<std::string, long> counts_;
};

}  // namespace mylab

SIMC_REGISTER(simc::ISpectrometer, "toy_window", mylab::ToyWindowSpectrometer)
