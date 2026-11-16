/**
 * @file
 * @ingroup output
 * @brief simc::RootTreeSink: sink "root_tree", the Fortran-compatible ntuple (TTree "h10", one Float_t branch per column).
 *
 * Writes the Fortran-compatible ntuple tree: TTree "h10" with one Float_t branch per
 * column, named and ordered exactly as NtupleInit.f / make_root_tree produce it
 * (46 columns for the (e,e'p) family; see NtupleLayout.h for the mapping).
 * Only successful events are written (results_write.f:34, "if (.not.success) return").
 *
 * JSON config (all optional):
 *   {"file": "<path>.root"}   default: <output_dir>/<run_name>.root
 *   {"tree": "h10"}
 */
#ifndef SIMC_OUTPUT_ROOTTREESINK_H
#define SIMC_OUTPUT_ROOTTREESINK_H

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "simc/core/RunSetup.h"
#include "simc/output/NtupleLayout.h"
#include "simc/stages/IOutputSink.h"
#include "simc/stages/IParallelSink.h"

class TFile;
class TTree;

namespace simc {

/// @addtogroup output
/// @{

/// @brief Sink "root_tree": writes one tree entry per successful event in the Fortran column layout (NtupleLayout.h).
///
/// Owns its TFile. Not thread-safe; a threaded run uses one worker copy per thread (makeWorkerSink).
class RootTreeSink : public IOutputSink, public IParallelSink {
public:
    /// Creates an unconfigured sink.
    RootTreeSink();
    /// Closes the file (a worker sink also removes its temporary file).
    ~RootTreeSink() override;
    std::string name() const override { return "root_tree"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    void beginRun(const RunSetup& setup) override;
    void fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
              bool success) override;
    void endRun(const RunSummary& summary) override;
    /// P15.4 (--threads N): a worker writes its tree to the temporary file <file>.w<k>.tmp; mergeWorkerSink() appends
    /// its entries to this tree (worker order) and removes the file.
    /// @param worker worker index (0-based).
    /// @param nworkers number of workers.
    /// @return a new sink for that worker.
    std::unique_ptr<IOutputSink> makeWorkerSink(int worker, int nworkers) const override;
    void mergeWorkerSink(IOutputSink& worker) override;

    /// @return the output file path.
    const std::string& fileName() const { return file_; }
    /// @return the number of entries written.
    long entries() const { return entries_; }

private:
    std::string file_, tree_name_{"h10"};
    RunSetup setup_;
    std::unique_ptr<TFile> tfile_;
    TTree* tree_{nullptr};                                  // owned by tfile_
    std::vector<double> ntu_;                               // static like the Fortran array (layout size)
    std::vector<float> out_;                                // make_root_tree: real*4 ntup_out
    long entries_{0};
    bool worker_{false};   ///< a worker sink of the parallel loop (temporary file, removed at destruction)
    void closeFile();
};

/// @}

}  // namespace simc

#endif
