/**
 * @file
 * @ingroup validation
 * @brief simc::ScriptedRandom: a random source that replays scripted numbers, for routine-level (L1) tests.
 *
 * It mirrors the stubs of the reference test harness (fortran_harness/common/rng_queue.f, rng_stub_gauss.f,
 * rng_stub_grnd.f):
 *  - two FIFO queues, uniforms (grnd) and gaussians (gauss1);
 *  - every call is counted (nUniform(), nGauss()); the last gauss1 argument is recorded;
 *  - asking for more numbers than were supplied sets error() and returns NaN.
 *
 * Gaussian modes:
 *  - Scripted (default): gauss1(nsig) pops the gaussian queue; nsig is recorded, NOT applied (like rng_stub_gauss.f).
 *  - Algorithm: gauss1(nsig) runs the real gauss1 algorithm (simc::gauss1Polar) on the uniform queue (like the
 *    gauss1_fixedrng reference table: real gauss1, stubbed grnd).
 */
#ifndef SIMC_CORE_RANDOM_SCRIPTEDRANDOM_H
#define SIMC_CORE_RANDOM_SCRIPTEDRANDOM_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "simc/stages/IRandom.h"

namespace simc {

/** @addtogroup validation
 *  @{ */

/**
 * @brief Random source that replays fixed numbers and counts the calls (see the file description).
 *
 * Not registered in the plug-in registry: tests create it directly. Not thread-safe.
 */
class ScriptedRandom : public IRandom {
public:
    /// How gauss1() produces its numbers.
    enum class GaussMode {
        Scripted,    ///< Pop the gaussian queue (the argument is only recorded).
        Algorithm    ///< Run simc::gauss1Polar on the uniform queue.
    };

    /**
     * @brief Creates a source with empty queues.
     * @param mode The Gaussian mode.
     */
    explicit ScriptedRandom(GaussMode mode = GaussMode::Scripted) : mode_(mode) {}

    std::string name() const override { return "scripted"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}

    /**
     * @brief Replaces both queues and resets all counters and the error flag.
     * @param gaussians The numbers gauss1() returns in Scripted mode, in order.
     * @param uniforms The numbers uniform() returns, in order.
     * @note Origin: reference test harness hrng_load(ng, g, nu, u).
     */
    void load(std::vector<double> gaussians, std::vector<double> uniforms = {});
    /**
     * @brief Changes the Gaussian mode.
     * @param m The new mode (the queues and counters are kept).
     */
    void setMode(GaussMode m) { mode_ = m; }

    /**
     * @brief The next scripted uniform.
     * @return The next number of the uniform queue, or NaN (and error() set) if the queue is empty.
     */
    double uniform() override;
    /**
     * @brief The next Gaussian, by the current mode.
     * @param nsig_max The truncation (recorded as lastNsig(); applied only in Algorithm mode).
     * @return The next Gaussian, or NaN (and error() set) if a queue ran dry.
     */
    double gauss1(double nsig_max) override;
    /// Does nothing (a scripted source has no seed).
    void seed(std::uint64_t) override {}

    /// @return The number of uniform() (grnd) calls, including those inside gauss1 in Algorithm mode.
    std::size_t nUniform() const { return ncu_; }
    /// @return The number of gauss1() calls.
    std::size_t nGauss() const { return ncg_; }
    /// @return The argument of the last gauss1() call (0 if none).
    double lastNsig() const { return nsig_last_; }
    /// @return True if a queue ran dry since the last load().
    bool error() const { return err_; }
    /// @return The number of scripted gaussians not consumed yet.
    std::size_t gaussiansLeft() const { return g_.size() - ig_; }
    /// @return The number of scripted uniforms not consumed yet.
    std::size_t uniformsLeft() const { return u_.size() - iu_; }

private:
    GaussMode mode_;
    std::vector<double> g_, u_;
    std::size_t ig_{0}, iu_{0}, ncg_{0}, ncu_{0};
    double nsig_last_{0};
    bool err_{false};
};

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_SCRIPTEDRANDOM_H
