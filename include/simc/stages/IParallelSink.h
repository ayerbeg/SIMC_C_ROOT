/**
 * @file
 * @ingroup stages
 * @brief simc::IParallelSink: optional extension of simc::IOutputSink for the parallel event loop (`--threads N`).
 *
 * With N > 1 worker threads every worker runs its own copy of the event-loop stages and fills its OWN sinks; the run's
 * (master) sinks then take the workers' content in the fixed worker order 0..N-1 and write the outputs once. A worker
 * sink sees only the fill() arguments (main, vertex, orig, recon, success): the simc::StageContext collaborators it got
 * in init() belong to the master run. A sink that does not implement IParallelSink makes `--threads N > 1` fail at
 * start-up with a clear message (the single-thread run is unaffected). Task P15.4, docs/PERFORMANCE.md.
 */
#ifndef SIMC_STAGES_IPARALLELSINK_H
#define SIMC_STAGES_IPARALLELSINK_H

#include <memory>

#include "simc/stages/IOutputSink.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Mix-in interface of an output sink that supports the parallel event loop.
 *
 * Implemented in addition to simc::IOutputSink (multiple inheritance); the driver finds it with dynamic_cast.
 * makeWorkerSink() is called on the master sink before the workers start; mergeWorkerSink() is called on the master
 * sink, from one thread, after every worker has finished.
 */
class IParallelSink {
public:
    /// Virtual destructor (deleted through the simc::IOutputSink pointer of the implementation).
    virtual ~IParallelSink() = default;
    /**
     * @brief Creates the sink of one worker.
     *
     * The worker sink is a new, configured sink of the same kind; it gets beginRun(setup) and fill() only, never
     * endRun(): it must not write the run's output files.
     * @param worker Index of the worker, 0..nworkers-1.
     * @param nworkers Number of workers.
     * @return The worker sink (owned by the caller), or nullptr if fill() does nothing for this sink (a sink that only
     *         works at endRun, such as run_summary).
     */
    virtual std::unique_ptr<IOutputSink> makeWorkerSink(int worker, int nworkers) const = 0;
    /**
     * @brief Adds the content of a worker sink to this (master) sink.
     *
     * Called once per worker, in the order 0..n-1, after every worker has finished and before endRun().
     * @param worker A sink made by makeWorkerSink() of THIS sink.
     */
    virtual void mergeWorkerSink(IOutputSink& worker) = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_IPARALLELSINK_H
