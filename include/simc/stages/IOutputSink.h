/**
 * @file
 * @ingroup stages
 * @brief simc::IOutputSink: an output of the run (ROOT tree, summary JSON, histograms, reports).
 * @note Origin: Fortran SIMC results_write.f, NtupleInit.f, simc.f:report (l.644), util/root_tree/make_root_tree.
 */
#ifndef SIMC_STAGES_IOUTPUTSINK_H
#define SIMC_STAGES_IOUTPUTSINK_H

#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief An output of the run: receives every try and the end-of-run summary.
 *
 * Call sequence (after configure()/init()): beginRun() once, fill() for every try, endRun() once. Sinks write their
 * files at endRun() (or open them at beginRun()). For `--threads N > 1` a sink must also implement
 * simc::IParallelSink. Implementations: simc::RootTreeSink, simc::RunSummarySink, simc::GenHistSink and the report
 * sinks of the driver.
 */
class IOutputSink : public IStage {
public:
    /**
     * @brief Starts the run (open files, book histograms).
     * @param setup The run set-up; valid until endRun() returns.
     * @throws simc::Error if an output cannot be created (implementation-specific).
     */
    virtual void beginRun(const RunSetup& setup) = 0;
    /**
     * @brief Called for every tried event.
     * @param main The weights, target and spectrometer quantities of the try.
     * @param vertex The vertex kinematics.
     * @param orig The kinematics after radiation and energy loss.
     * @param recon The reconstructed kinematics.
     * @param success True if the event passed every aperture and cut and contributes (weighted) to the results.
     */
    virtual void fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                      bool success) = 0;
    /**
     * @brief Ends the run: normalise and write the outputs.
     * @param summary The end-of-run numbers (counts, normalisation factor, stop counters).
     */
    virtual void endRun(const RunSummary& summary) = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_IOUTPUTSINK_H
