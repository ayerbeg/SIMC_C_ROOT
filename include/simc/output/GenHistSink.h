/**
 * @file
 * @ingroup output
 * @brief simc::GenHistSink: sink "gen_hist", the Fortran diagnostic histograms H%geni, H%gen, H%RECON as TH1D.
 *
 * The Fortran diagnostic histograms H%geni, H%gen, H%RECON (histograms.inc, nHbins=50),
 * filled exactly like simc.f l.254-286 with inc() (simc.f:624-640), written as TH1D to
 * <run>.gen.root. Axes come from RunSetup::hist (init.f:limits_init l.515-570).
 *
 *   geni_*  every tried event, weight 1: vertex e/p delta, yptar, -xptar, Em, Pm
 *   gen_*   successful events, weight 1: vertex e/p delta, yptar, -xptar, Em   (no Pm, as the Fortran)
 *   recon_* successful events: main.RECON e/p delta, yptar, xptar (weight main.weight; NOTE: xptar
 *           NOT negated although the axis starts at -gen.xptar.max — FORTRAN-DIFF-CANDIDATE simc.f:274,277),
 *           recon Em, Pm (weight 1)
 * Values outside the 50 bins are dropped (inc() has no under/overflow).
 *
 * JSON config (optional): {"file": "<path>.gen.root"}
 */
#ifndef SIMC_OUTPUT_GENHISTSINK_H
#define SIMC_OUTPUT_GENHISTSINK_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "simc/core/RunSetup.h"
#include "simc/stages/IOutputSink.h"
#include "simc/stages/IParallelSink.h"

class TH1D;

namespace simc {

/// @addtogroup output
/// @{

/// @brief Sink "gen_hist": fills the 50-bin histograms like simc.f and writes them to the .gen.root file.
///
/// Owns its histograms. Not thread-safe; a threaded run uses one worker copy per thread (makeWorkerSink).
class GenHistSink : public IOutputSink, public IParallelSink {
public:
    static constexpr int kNHbins = 50;   ///< number of bins (histograms.inc nHbins).

    /// Creates an unconfigured sink.
    GenHistSink();
    /// Destroys the histograms.
    ~GenHistSink() override;
    std::string name() const override { return "gen_hist"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    void beginRun(const RunSetup& setup) override;
    void fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
              bool success) override;
    void endRun(const RunSummary& summary) override;
    /// P15.4 (--threads N): worker histograms in memory; merge = bin contents summed, errors in quadrature, entries summed.
    /// @param worker worker index (0-based).
    /// @param nworkers number of workers.
    /// @return a new in-memory sink for that worker.
    std::unique_ptr<IOutputSink> makeWorkerSink(int worker, int nworkers) const override;
    void mergeWorkerSink(IOutputSink& worker) override;

    /// Fortran inc(): ibin = nint(0.5+(val-min)/bin); counted only for 1 <= ibin <= nHbins.
    /// @param val the value.
    /// @param ax the histogram axis.
    /// @return the 1-based bin (may be outside 1..kNHbins: then the value is not counted).
    static int FortranBin(double val, const HistAxis& ax);

    /// @return the output file path.
    const std::string& fileName() const { return file_; }
    /// @param name histogram name (e.g. "geni_e_delta").
    /// @return the histogram (owned by the sink), or nullptr if there is none of that name.
    const TH1D* hist(const std::string& name) const;

private:
    bool rhoGenFromOrig_{false};   ///< D-057 F18
    void book(const std::string& name, const HistAxis& ax);
    /// P15.4: the histograms are resolved once at beginRun (slot = set*8 + variable, in booking order) instead of two
    /// std::map<string> lookups per inc() call; same TH1D operations, same results.
    struct Slot {
        TH1D* h{nullptr};
        HistAxis ax;
    };
    void inc(int slot, double val, double weight);

    std::string file_;
    RunSetup setup_;
    std::map<std::string, std::unique_ptr<TH1D>> h_;
    std::map<std::string, HistAxis> ax_;
    std::vector<Slot> slots_;
};

/// @}

}  // namespace simc

#endif
