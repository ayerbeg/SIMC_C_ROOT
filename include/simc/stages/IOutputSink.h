// include/simc/stages/IOutputSink.h
// Output: ntuple/tree, .hist-like report, generated histograms.
// Fortran: results_write.f, NtupleInit.f, simc.f:report (l.644), util/root_tree/make_root_tree.
#ifndef SIMC_STAGES_IOUTPUTSINK_H
#define SIMC_STAGES_IOUTPUTSINK_H

#include "simc/stages/IStage.h"

namespace simc {

class IOutputSink : public IStage {
public:
    virtual void beginRun(const RunSetup& setup) = 0;
    /// Called for every tried event (success tells whether it contributes).
    virtual void fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                      bool success) = 0;
    virtual void endRun(const RunSummary& summary) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_IOUTPUTSINK_H
