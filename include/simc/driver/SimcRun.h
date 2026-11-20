/**
 * @file
 * @ingroup driver
 * @brief simc::driver::SimcRun: one complete SIMC run, the counterpart of the simc.f main program.
 *
 * SimcRun: one complete SIMC run (task P7.2) — the C++ counterpart of the simc.f main program (l.1-623).
 *
 *   SimcRun run(params, pipeline, options);   // nothing happens yet
 *   run.initialise();     // simc.f l.53-160: RunSetup (dbase_read, init.f, radc_init), the stages, calculate_central
 *   while (!run.done()) run.step();           // simc.f l.169-352, one try per step
 *   run.finish();         // simc.f l.354-620: normalisation, resolutions, sinks (.root, summary, .gen.root, .hist)
 * or simply run.run().
 *
 * Ownership and data flow (docs/manual/markdown/reference/DRIVER.md):
 *   * the run owns every stage (unique_ptr, created by name through the Registry from the PipelineConfig) and hands
 *     the collaborators to them through a StageContext at init(), as INTERFACES (beam -> target, reaction -> target +
 *     radiator, montecarlo -> target + arms); the driver never casts a stage to a concrete class;
 *   * the set-up takes its start-up target losses from the pipeline's target model (ITargetLossProvider);
 *   * the per-event state is ONE EventState (main, vertex, orig, recon) handed to the stages on every try; the
 *     Fortran keeps stale values across tries in these structures and SIMC++ reproduces that (EventState.h);
 *   * the run-level accumulators are RunTotals; nothing is global.
 * The order of the stage calls and of the random-number draws is the Fortran one, so that a ScriptedRandom fed with
 * a captured RANLUX stream reproduces a Fortran run bit for bit (tests/l1 simcrun replay).
 */
#ifndef SIMC_DRIVER_SIMCRUN_H
#define SIMC_DRIVER_SIMCRUN_H

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/driver/Central.h"
#include "simc/driver/EventState.h"
#include "simc/driver/MonteCarlo.h"
#include "simc/driver/Pipeline.h"
#include "simc/driver/RunReport.h"
#include "simc/reactions/EventChain.h"
#include "simc/stages/IBeam.h"
#include "simc/stages/IOutputSink.h"
#include "simc/stages/IPhysicsModel.h"
#include "simc/stages/IRadiator.h"
#include "simc/stages/IRandom.h"
#include "simc/stages/IReaction.h"
#include "simc/stages/ITarget.h"

namespace simc {
namespace driver {

/// @addtogroup driver
/// @{

/// Options of a run that are not physics parameters (output names, overrides, logging, compatibility, threads).
struct RunOptions {
    std::string run_name{"simc"};   ///< Fortran `base`: output file stem
    std::string output_dir{"."};    ///< Fortran outfiles/ and worksim/
    std::string data_dir;           ///< Fortran working-directory data files (h2.theory, benharsf_*.dat, ...);
                                    ///< "" = the default data directory (simc::DataRoot(): SIMC_DATA_DIR, installed prefix, or source tree)
    bool override_ngen{false};  ///< true: use #ngen instead of the input ngen.
    long ngen{0};  ///< number of events when #override_ngen is set (> 0 successes, < 0 tries).
    bool override_seed{false};  ///< true: use #seed instead of the input random_seed.
    long seed{0};  ///< random seed when #override_seed is set.
    std::ostream* log{nullptr};     ///< the Fortran unit-6 messages of the driver (nullptr = std::cout)
    bool quiet_stages{false};       ///< silence the stage warnings (write(6,*) inside the physics routines)
    /// D-031: false (default) = clear Fortran bugs corrected; true = exact Fortran behaviour (--fortran-compat,
    /// JSON "run": {"fortran_compat": true}). Passed to the set-up (RunSetup::fortran_compat).
    bool fortran_compat{false};
    /// Task U.2 (D-082): layout of the tree h10. Upstream (default) = upstream SIMC 59ccdf6; Legacy = SIMC 294f662
    /// (--tree-layout legacy, JSON "output": {"layout": "legacy"}). Passed to RunSetup::tree_layout.
    TreeLayout tree_layout{TreeLayout::Upstream};
    /// Task U.R: upstream commit reproduced by fortran_compat (--compat-target 59ccdf6|294f662, JSON
    /// "run": {"compat_target": ...}); default 59ccdf6. Passed to RunSetup::compat_target.
    CompatTarget compat_target{CompatTarget::U59ccdf6};
    /// P15.4: worker threads of the event loop (--threads N); 0 = from the pipeline "threads" (default 1). N > 1:
    /// runThreaded() (docs/manual/markdown/reference/PERFORMANCE.md); refused with fortran_compat and with an rng other than mt64.
    int threads{0};
};

/// Fortran: simc.f l.226-238, the SPedge test on the reconstructed event (`pass_cuts`). With setup.fortran_compat the
/// hadron delta upper bound is the ELECTRON-arm limit as in simc.f:237 (FORTRAN-BUG D-006); otherwise the hadron one.
/// @param setup the run set-up (SPedge limits).
/// @param recon the reconstructed event.
/// @return true if the event passes the cuts.
bool PassCuts(const RunSetup& setup, const Event& recon);

/// @brief One complete SIMC run: set-up, stages, event loop, normalisation and outputs (see the file description).
///
/// Owns all stages and sinks. A SimcRun is used by one thread; runThreaded() creates its own worker stages.
class SimcRun {
public:
    /// @brief Store the inputs; nothing is built until initialise().
    /// @param params the run parameters (from the .inp or JSON).
    /// @param pipeline the "pipeline" JSON section (may be empty).
    /// @param options run options.
    SimcRun(SimcParameters params, const Json& pipeline, RunOptions options);
    /// Destroys the stages and sinks.
    ~SimcRun();

    /// Use this random source instead of the pipeline "rng" stage (e.g. ScriptedRandom for a replay).
    /// @param rng the random source (ownership taken).
    void setRandom(std::unique_ptr<IRandom> rng);
    /// Add a sink created outside the registry (tests, user code). Before initialise(). The run takes ownership;
    /// the returned reference stays valid for the lifetime of the run.
    /// @param sink the sink (ownership taken).
    /// @return a reference to the sink.
    IOutputSink& addSink(std::unique_ptr<IOutputSink> sink);

    /// @brief simc.f l.53-160: build the RunSetup, create and initialise the stages, compute the central event.
    /// @throws std::invalid_argument, std::runtime_error or FortranStop on an invalid configuration.
    void initialise();
    /// @return true when the requested number of events (or tries) is reached.
    bool done() const;
    /// One try of the event loop (simc.f l.169-352).
    void step();
    /// simc.f l.354-620: normalisation, resolutions and the sinks' outputs.
    void finish();
    /// initialise(), then the event loop on threadCount() threads (step() loop or runThreaded()), then finish().
    void run();

    // ---- parallel event loop (task P15.4, docs/manual/markdown/reference/PERFORMANCE.md) ------------------------------------------------
    /// The number of event-loop threads of this run: RunOptions::threads if > 0, else the pipeline "threads", else 1
    /// (after initialise()).
    /// @return the thread count.
    int threadCount() const;
    /// Why `n` > 1 threads cannot be used for this run ("" if they can): fortran_compat, an rng other than mt64 or a
    /// random source set with setRandom(), a sink without IParallelSink. After initialise().
    /// @param n requested number of threads.
    /// @return the reason, or "" if @p n threads can be used.
    std::string threadsRefusal(int n) const;
    /// The event loop on n worker threads (n > 1), each with its own stages (built from the same set-up) and its own
    /// mt64 stream Mt64Random::seedStream(seed, k); worker k makes ngen/n (+1 for k < ngen%n) of the successes (or of
    /// the tries for ngen < 0). The totals, stop counters and sinks are then merged in the worker order 0..n-1 and
    /// finish() writes the outputs. A given (seed, n) gives the same outputs on every run. Throws
    /// std::invalid_argument if threadsRefusal(n) is not empty. Calls initialise() and finish().
    /// @param n number of worker threads (> 1).
    void runThreaded(int n);

    // ---- access -------------------------------------------------------------------------------------------
    /// @return the run set-up (after initialise()).
    const RunSetup& setup() const { return setup_; }
    /// @return the pipeline configuration.
    const PipelineConfig& pipeline() const { return pipe_; }
    /// @return the per-event state of the last try.
    const EventState& event() const { return ev_; }
    /// @return the run accumulators.
    const RunTotals& totals() const { return report_.totals; }
    /// @return the run report.
    const RunReport& report() const { return report_; }
    /// @return the central event.
    const CentralResult& central() const { return central_; }
    /// The run's physics model and reaction (after initialise; L1 tests read their tables / Fortran state).
    /// @return the physics model (owned by the run).
    const IPhysicsModel& model() const { return *model_; }
    /// @return the reaction (owned by the run).
    const IReaction& reaction() const { return *reaction_; }
    /// @return the random source of the run (owned by the run).
    IRandom& rng() { return *rngPtr_; }
    /// @return the radiator (owned by the run).
    IRadiator& radiator() { return *radiator_; }
    /// @return the arm transport (owned by the run).
    MonteCarlo& montecarlo() { return *mc_; }
    /// @return the sinks created from the pipeline.
    const std::vector<std::unique_ptr<IOutputSink>>& sinks() const { return sinks_; }

private:
    void resetForTry();
    std::vector<StopCounters> collectStops() const;
    void log(const std::string& line);

    SimcParameters params_;
    Json pipelineJson_;
    RunOptions opt_;
    RunSetup setup_;
    PipelineConfig pipe_;
    bool initialised_{false}, finished_{false};

    std::unique_ptr<IRandom> rng_;
    IRandom* rngPtr_{nullptr};
    std::unique_ptr<ITarget> target_;
    std::unique_ptr<IRadiator> radiator_;
    std::unique_ptr<IBeam> beam_;
    std::unique_ptr<IReaction> reaction_;
    std::unique_ptr<IPhysicsModel> model_;
    std::unique_ptr<MonteCarlo> mc_;
    std::vector<std::unique_ptr<IOutputSink>> sinks_, extraSinks_;

    CentralResult central_;
    EventState ev_;
    RunReport report_;
    long ninform_{1000};
    // P15.4: set by runThreaded() for finish(): the merged stop counters and the reaction of the last worker
    bool mergedStopsSet_{false};
    std::vector<StopCounters> mergedStops_;
    const IReaction* finalStateReaction_{nullptr};
};

/// @}

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_SIMCRUN_H
