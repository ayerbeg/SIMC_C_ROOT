// include/simc/output/GenHistSink.h
//
// The Fortran diagnostic histograms H%geni, H%gen, H%RECON (histograms.inc, nHbins=50),
// filled exactly like simc.f l.254-286 with inc() (simc.f:624-640), written as TH1D to
// <run>.gen.root. Axes come from RunSetup::hist (init.f:limits_init l.515-570).
//
//   geni_*  every tried event, weight 1: vertex e/p delta, yptar, -xptar, Em, Pm
//   gen_*   successful events, weight 1: vertex e/p delta, yptar, -xptar, Em   (no Pm, as the Fortran)
//   recon_* successful events: main.RECON e/p delta, yptar, xptar (weight main.weight; NOTE: xptar
//           NOT negated although the axis starts at -gen.xptar.max — FORTRAN-QUIRK simc.f:274,277),
//           recon Em, Pm (weight 1)
// Values outside the 50 bins are dropped (inc() has no under/overflow).
//
// JSON config (optional): {"file": "<path>.gen.root"}
#ifndef SIMC_OUTPUT_GENHISTSINK_H
#define SIMC_OUTPUT_GENHISTSINK_H

#include <map>
#include <memory>
#include <string>

#include "simc/core/RunSetup.h"
#include "simc/stages/IOutputSink.h"

class TH1D;

namespace simc {

class GenHistSink : public IOutputSink {
public:
    static constexpr int kNHbins = 50;   // histograms.inc nHbins

    GenHistSink();
    ~GenHistSink() override;
    std::string name() const override { return "gen_hist"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    void beginRun(const RunSetup& setup) override;
    void fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
              bool success) override;
    void endRun(const RunSummary& summary) override;

    /// Fortran inc(): ibin = nint(0.5+(val-min)/bin); counted only for 1 <= ibin <= nHbins.
    static int FortranBin(double val, const HistAxis& ax);

    const std::string& fileName() const { return file_; }
    const TH1D* hist(const std::string& name) const;

private:
    void book(const std::string& name, const HistAxis& ax);
    void inc(const std::string& name, double val, double weight);

    std::string file_;
    RunSetup setup_;
    std::map<std::string, std::unique_ptr<TH1D>> h_;
    std::map<std::string, HistAxis> ax_;
};

}  // namespace simc

#endif
