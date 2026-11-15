// include/simc/stages/IParallelSink.h
//
// Optional extension of IOutputSink for the parallel event loop (task P15.4, `simc++ --threads N`, docs/PERFORMANCE.md).
// With N > 1 worker threads every worker runs its own copy of the event-loop stages and fills its OWN sinks; the run's
// (master) sinks then take the workers' content in the fixed worker order 0..N-1 and write the outputs once.
//
//   makeWorkerSink(k, n)  a new, configured sink of the same kind for worker k (it gets beginRun(setup) and fill() only,
//                         never endRun(): it must not write the run's output files); nullptr if fill() does nothing
//                         (a sink that only works at endRun, such as run_summary)
//   mergeWorkerSink(w)    add the content of worker sink `w` (made by makeWorkerSink of THIS sink) to this sink; called
//                         on the master sink after every worker has finished, once per worker in the order 0..n-1, before
//                         endRun()
//
// A worker sink sees only the fill() arguments (main, vertex, orig, recon, success): the StageContext collaborators
// it got in init() belong to the master run. A sink that does not implement IParallelSink makes `--threads N > 1`
// fail at start-up with a clear message (the single-thread run is unaffected).
#ifndef SIMC_STAGES_IPARALLELSINK_H
#define SIMC_STAGES_IPARALLELSINK_H

#include <memory>

#include "simc/stages/IOutputSink.h"

namespace simc {

class IParallelSink {
public:
    virtual ~IParallelSink() = default;
    virtual std::unique_ptr<IOutputSink> makeWorkerSink(int worker, int nworkers) const = 0;
    virtual void mergeWorkerSink(IOutputSink& worker) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_IPARALLELSINK_H
