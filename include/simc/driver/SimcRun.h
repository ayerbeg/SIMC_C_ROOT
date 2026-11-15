// include/simc/driver/SimcRun.h
//
// SimcRun: one complete SIMC run (task P7.2) — the C++ counterpart of the simc.f main program (l.1-623).
//
//   SimcRun run(params, pipeline, options);   // nothing happens yet
//   run.initialise();     // simc.f l.53-160: RunSetup (dbase_read, init.f, radc_init), the stages, calculate_central
//   while (!run.done()) run.step();           // simc.f l.169-352, one try per step
//   run.finish();         // simc.f l.354-620: normalisation, resolutions, sinks (.root, summary, .gen.root, .hist)
// or simply run.run().
//
// Ownership and data flow (docs/DRIVER.md):
//   * the run owns every stage (unique_ptr, created by name through the Registry from the PipelineConfig) and hands
//     the collaborators to them through a StageContext at init(), as INTERFACES (beam -> target, reaction -> target +
//     radiator, montecarlo -> target + arms); the driver never casts a stage to a concrete class;
//   * the set-up takes its start-up target losses from the pipeline's target model (ITargetLossProvider);
//   * the per-event state is ONE EventState (main, vertex, orig, recon) handed to the stages on every try; the
//     Fortran keeps stale values across tries in these structures and SIMC++ reproduces that (EventState.h);
//   * the run-level accumulators are RunTotals; nothing is global.
// The order of the stage calls and of the random-number draws is the Fortran one, so that a ScriptedRandom fed with
// a captured RANLUX stream reproduces a Fortran run bit for bit (tests/l1 simcrun replay).
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

struct RunOptions {
    std::string run_name{"simc"};   ///< Fortran `base`: output file stem
    std::string output_dir{"."};    ///< Fortran outfiles/ and worksim/
    std::string data_dir;           ///< Fortran working-directory data files (h2.theory, benharsf_*.dat, ...);
                                    ///< "" = data/spectral of the source tree (SIMC_SPECTRAL_DATA_DIR)
    bool override_ngen{false};
    long ngen{0};
    bool override_seed{false};
    long seed{0};
    std::ostream* log{nullptr};     ///< the Fortran unit-6 messages of the driver (nullptr = std::cout)
    bool quiet_stages{false};       ///< silence the stage warnings (write(6,*) inside the physics routines)
    /// D-031: false (default) = clear Fortran bugs corrected; true = exact Fortran behaviour (--fortran-compat,
    /// JSON "run": {"fortran_compat": true}). Passed to the set-up (RunSetup::fortran_compat).
    bool fortran_compat{false};
    /// P15.4: worker threads of the event loop (--threads N); 0 = from the pipeline "threads" (default 1). N > 1:
    /// runThreaded() (docs/PERFORMANCE.md); refused with fortran_compat and with an rng other than mt64.
    int threads{0};
};

/// Fortran: simc.f l.226-238, the SPedge test on the reconstructed event (`pass_cuts`). With setup.fortran_compat the
/// hadron delta upper bound is the ELECTRON-arm limit as in simc.f:237 (FORTRAN-BUG D-006); otherwise the hadron one.
bool PassCuts(const RunSetup& setup, const Event& recon);

class SimcRun {
public:
    SimcRun(SimcParameters params, const Json& pipeline, RunOptions options);
    ~SimcRun();

    /// Use this random source instead of the pipeline "rng" stage (e.g. ScriptedRandom for a replay).
    void setRandom(std::unique_ptr<IRandom> rng);
    /// Add a sink created outside the registry (tests, user code). Before initialise(). The run takes ownership;
    /// the returned reference stays valid for the lifetime of the run.
    IOutputSink& addSink(std::unique_ptr<IOutputSink> sink);

    void initialise();
    bool done() const;
    void step();
    void finish();
    /// initialise(), then the event loop on threadCount() threads (step() loop or runThreaded()), then finish().
    void run();

    // ---- parallel event loop (task P15.4, docs/PERFORMANCE.md) ------------------------------------------------
    /// The number of event-loop threads of this run: RunOptions::threads if > 0, else the pipeline "threads", else 1
    /// (after initialise()).
    int threadCount() const;
    /// Why `n` > 1 threads cannot be used for this run ("" if they can): fortran_compat, an rng other than mt64 or a
    /// random source set with setRandom(), a sink without IParallelSink. After initialise().
    std::string threadsRefusal(int n) const;
    /// The event loop on n worker threads (n > 1), each with its own stages (built from the same set-up) and its own
    /// mt64 stream Mt64Random::seedStream(seed, k); worker k makes ngen/n (+1 for k < ngen%n) of the successes (or of
    /// the tries for ngen < 0). The totals, stop counters and sinks are then merged in the worker order 0..n-1 and
    /// finish() writes the outputs. A given (seed, n) gives the same outputs on every run. Throws
    /// std::invalid_argument if threadsRefusal(n) is not empty. Calls initialise() and finish().
    void runThreaded(int n);

    // ---- access -------------------------------------------------------------------------------------------
    const RunSetup& setup() const { return setup_; }
    const PipelineConfig& pipeline() const { return pipe_; }
    const EventState& event() const { return ev_; }
    const RunTotals& totals() const { return report_.totals; }
    const RunReport& report() const { return report_; }
    const CentralResult& central() const { return central_; }
    /// The run's physics model and reaction (after initialise; L1 tests read their tables / Fortran state).
    const IPhysicsModel& model() const { return *model_; }
    const IReaction& reaction() const { return *reaction_; }
    IRandom& rng() { return *rngPtr_; }
    IRadiator& radiator() { return *radiator_; }
    MonteCarlo& montecarlo() { return *mc_; }
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

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_SIMCRUN_H
