// include/simc/core/random/ScriptedRandom.h
//
// IRandom that replays scripted numbers, for L1 fixed-RNG tests. It mirrors the Fortran
// harness stubs (fortran_harness/common/rng_queue.f, rng_stub_gauss.f, rng_stub_grnd.f):
//   * two FIFO queues, uniforms (grnd) and gaussians (gauss1);
//   * every call is counted (nUniform(), nGauss()); the last gauss1 argument is recorded;
//   * asking for more numbers than were supplied sets error() and returns NaN.
// Gaussian mode:
//   * Scripted  (default): gauss1(nsig) pops the gaussian queue; nsig is recorded, NOT applied
//                          (like rng_stub_gauss.f).
//   * Algorithm           : gauss1(nsig) runs the real gauss1.f algorithm (gauss1Polar) on the
//                          uniform queue (like the gauss1_fixedrng table: real gauss1.o, stubbed grnd).
#ifndef SIMC_CORE_RANDOM_SCRIPTEDRANDOM_H
#define SIMC_CORE_RANDOM_SCRIPTEDRANDOM_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "simc/stages/IRandom.h"

namespace simc {

class ScriptedRandom : public IRandom {
public:
    enum class GaussMode { Scripted, Algorithm };

    explicit ScriptedRandom(GaussMode mode = GaussMode::Scripted) : mode_(mode) {}

    std::string name() const override { return "scripted"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}

    /// Fortran harness: hrng_load(ng, g, nu, u). Replaces both queues and resets all counters.
    void load(std::vector<double> gaussians, std::vector<double> uniforms = {});
    void setMode(GaussMode m) { mode_ = m; }

    double uniform() override;
    double gauss1(double nsig_max) override;
    void seed(std::uint64_t) override {}

    std::size_t nUniform() const { return ncu_; }   ///< grnd calls (including those inside gauss1 in Algorithm mode)
    std::size_t nGauss() const { return ncg_; }     ///< gauss1 calls
    double lastNsig() const { return nsig_last_; }  ///< argument of the last gauss1 call (0 if none)
    bool error() const { return err_; }             ///< a queue ran dry
    std::size_t gaussiansLeft() const { return g_.size() - ig_; }
    std::size_t uniformsLeft() const { return u_.size() - iu_; }

private:
    GaussMode mode_;
    std::vector<double> g_, u_;
    std::size_t ig_{0}, iu_{0}, ncg_{0}, ncu_{0};
    double nsig_last_{0};
    bool err_{false};
};

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_SCRIPTEDRANDOM_H
