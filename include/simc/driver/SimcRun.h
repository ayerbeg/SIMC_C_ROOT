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
    std::string data_dir;           ///< Fortran working directory files (h2.theory, ...); "" = current directory
    bool override_ngen{false};
    long ngen{0};
    bool override_seed{false};
    long seed{0};
    std::ostream* log{nullptr};     ///< the Fortran unit-6 messages of the driver (nullptr = std::cout)
    bool quiet_stages{false};       ///< silence the stage warnings (write(6,*) inside the physics routines)
};

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
    void run();

    // ---- access -------------------------------------------------------------------------------------------
    const RunSetup& setup() const { return setup_; }
    const PipelineConfig& pipeline() const { return pipe_; }
    const EventState& event() const { return ev_; }
    const RunTotals& totals() const { return report_.totals; }
    const RunReport& report() const { return report_; }
    const CentralResult& central() const { return central_; }
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

    reactions::MainContext ctx_;
    CentralResult central_;
    EventState ev_;
    RunReport report_;
    long ninform_{1000};
};

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_SIMCRUN_H
