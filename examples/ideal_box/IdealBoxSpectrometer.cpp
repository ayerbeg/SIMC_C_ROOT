// examples/ideal_box/IdealBoxSpectrometer.cpp — an example ISpectrometer plug-in (task P7.4, docs/EXTENDING.md).
//
// "ideal_box": a spectrometer with PERFECT resolution and a rectangular acceptance in (delta, xptar, yptar) at the
// target. A track is accepted if
//     delta_min <= delta <= delta_max   (%)        xptar_min <= dx/dz <= xptar_max   (rad)
//     yptar_min <= dy/dz <= yptar_max   (rad)
// and then "reconstructed" exactly: delta, xptar, yptar and ytar (cm, the y at the plane z = 0 through the target
// centre) are returned unchanged. The focal-plane quantities are the target ones (unit optics), the path length is a
// constant, the resolution factor is 1. Nothing is random: transport() never draws a random number.
//
// It is selected for either arm purely through the JSON pipeline, e.g.
//     "pipeline": {"spectrometers": {"hadron": {"type": "ideal_box", "delta_min": -10, "delta_max": 10}}}
// JSON config (all optional; defaults = a large box):
//     delta_min/delta_max  [%]    (default -15 / +15)
//     xptar_min/xptar_max  [rad]  (default -0.1 / +0.1)
//     yptar_min/yptar_max  [rad]  (default -0.1 / +0.1)
//     path_length          [cm]   (default 0)
// Stop counters: trials, delta, xptar, yptar, successes (a track is counted at the first cut it fails, in that order).
//
// The whole plug-in is this file: one class derived from simc::ISpectrometer and one SIMC_REGISTER line. It is
// compiled into an OBJECT library linked by simc++ (examples/CMakeLists.txt) so that the registration is kept.
#include <map>
#include <string>

#include "simc/core/Registry.h"
#include "simc/stages/ISpectrometer.h"

namespace simc {
namespace examples {

class IdealBoxSpectrometer : public ISpectrometer {
public:
    std::string name() const override { return "ideal_box"; }

    void configure(const Json& config, const SimcParameters&) override {
        if (!config.is_object()) return;
        delta_ = {config.value("delta_min", delta_.lo), config.value("delta_max", delta_.hi)};
        xptar_ = {config.value("xptar_min", xptar_.lo), config.value("xptar_max", xptar_.hi)};
        yptar_ = {config.value("yptar_min", yptar_.lo), config.value("yptar_max", yptar_.hi)};
        pathLength_ = config.value("path_length", pathLength_);
    }

    void init(const StageContext&) override {
        counts_ = {{"trials", 0}, {"delta", 0}, {"xptar", 0}, {"yptar", 0}, {"successes", 0}};
    }

    SpecResult transport(const SpecInput& in, IRandom&) override {
        SpecResult r;
        r.m2_final = in.m2;
        r.Mh2_final = in.m2;
        ++counts_["trials"];
        if (!delta_.contains(in.delta)) return stop(r, "delta");
        if (!xptar_.contains(in.dxdz)) return stop(r, "xptar");
        if (!yptar_.contains(in.dydz)) return stop(r, "yptar");
        ++counts_["successes"];
        r.ok = true;
        r.delta_recon = in.delta;   // perfect resolution
        r.dxdz_recon = in.dxdz;
        r.dydz_recon = in.dydz;
        r.y_recon = in.y;
        r.x_fp = in.x;              // unit optics
        r.dx_fp = in.dxdz;
        r.y_fp = in.y;
        r.dy_fp = in.dydz;
        r.pathlen = pathLength_;
        r.resolution_factor = 1.0;
        return r;
    }

    const std::map<std::string, long>& stopCounters() const override { return counts_; }

private:
    struct Window {
        double lo, hi;
        bool contains(double x) const { return x >= lo && x <= hi; }
    };
    SpecResult& stop(SpecResult& r, const char* cut) {
        ++counts_[cut];
        r.stop_reason = cut;
        return r;
    }

    Window delta_{-15.0, 15.0}, xptar_{-0.1, 0.1}, yptar_{-0.1, 0.1};
    double pathLength_{0.0};
    std::map<std::string, long> counts_{{"trials", 0}, {"delta", 0}, {"xptar", 0}, {"yptar", 0}, {"successes", 0}};
};

}  // namespace examples
}  // namespace simc

SIMC_REGISTER(simc::ISpectrometer, "ideal_box", simc::examples::IdealBoxSpectrometer)
