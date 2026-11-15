// src/output/GenHistSink.cpp — see include/simc/output/GenHistSink.h.
// The "gen_hist" output sink: the generated (geni: every try; gen: successes, unweighted) and reconstructed (recon:
// successes, weighted) distributions of delta, yptar, xptar of both arms and Em, Pm, as the Fortran H%geni/H%gen/
// H%RECON histograms that it prints to the .gen file. Stored as TH1D in <run>.gen.root (the .gen text is written by
// the report, driver/RunReport.cpp, from these histograms). Binning and the bin errors follow simc.f:inc exactly.
// Fortran: simc.f:inc (l.624-640), fills at simc.f l.254-286, axes init.f l.515-570,
// text output .gen at simc.f l.534-605.
#include "simc/output/GenHistSink.h"

#include <cmath>
#include <stdexcept>

#include <TFile.h>
#include <TROOT.h>
#include <TH1D.h>

#include "simc/core/Event.h"
#include "simc/core/Registry.h"

namespace simc {

namespace {
// slots_ indices (beginRun books geni, gen, recon x e_delta, e_yptar, e_xptar, p_delta, p_yptar, p_xptar, Em, Pm)
constexpr int kGeni_e_delta = 0;
constexpr int kGeni_e_yptar = 1;
constexpr int kGeni_e_xptar = 2;
constexpr int kGeni_p_delta = 3;
constexpr int kGeni_p_yptar = 4;
constexpr int kGeni_p_xptar = 5;
constexpr int kGeni_Em = 6;
constexpr int kGeni_Pm = 7;
constexpr int kGen_e_delta = 8;
constexpr int kGen_e_yptar = 9;
constexpr int kGen_e_xptar = 10;
constexpr int kGen_p_delta = 11;
constexpr int kGen_p_yptar = 12;
constexpr int kGen_p_xptar = 13;
constexpr int kGen_Em = 14;
constexpr int kGen_Pm = 15;
constexpr int kRecon_e_delta = 16;
constexpr int kRecon_e_yptar = 17;
constexpr int kRecon_e_xptar = 18;
constexpr int kRecon_p_delta = 19;
constexpr int kRecon_p_yptar = 20;
constexpr int kRecon_p_xptar = 21;
constexpr int kRecon_Em = 22;
constexpr int kRecon_Pm = 23;
}  // namespace

GenHistSink::GenHistSink() = default;
GenHistSink::~GenHistSink() = default;

void GenHistSink::configure(const Json& config, const SimcParameters&) {
    if (config.is_object()) file_ = config.value("file", file_);
}

void GenHistSink::init(const StageContext&) {}

int GenHistSink::FortranBin(double val, const HistAxis& ax) {
    // ibin = nint(0.5+(val-hist%min)/hist%bin); Fortran nint rounds half away from zero (= std::lround).
    // The Fortran bins are not ROOT's bins: the TH1D is used as storage with the same bin index (1..kNHbins).
    const double x = 0.5 + (val - ax.min) / ax.bin;
    if (!(std::abs(x) < 1.0e9)) return 0;   // integer*4 range guard (also NaN)
    return static_cast<int>(std::lround(x));
}

void GenHistSink::book(const std::string& name, const HistAxis& ax) {
    auto h = std::make_unique<TH1D>(name.c_str(), name.c_str(), kNHbins, ax.min, ax.min + kNHbins * ax.bin);
    h->SetDirectory(nullptr);
    h->Sumw2();
    slots_.push_back(Slot{h.get(), ax});
    h_[name] = std::move(h);
    ax_[name] = ax;
}

void GenHistSink::beginRun(const RunSetup& setup) {
    setup_ = setup;
    rhoGenFromOrig_ = setup.flags.doing_rho && !setup.fortran_compat;   // D-057 F18
    if (file_.empty()) file_ = setup.output_dir + "/" + setup.run_name + ".gen.root";
    const auto& H = setup.hist;
    slots_.clear();
    for (const char* set : {"geni", "gen", "recon"}) {
        const std::string s = set;
        book(s + "_e_delta", H.e.delta);
        book(s + "_e_yptar", H.e.yptar);
        book(s + "_e_xptar", H.e.xptar);
        book(s + "_p_delta", H.p.delta);
        book(s + "_p_yptar", H.p.yptar);
        book(s + "_p_xptar", H.p.xptar);
        book(s + "_Em", H.Em);
        book(s + "_Pm", H.Pm);   // gen_Pm stays empty (the Fortran never fills H%gen%Pm)
    }
}

// Fortran: simc.f:inc. Bin content += weight, error = sqrt(error^2 + weight^2), entries + 1 per call in range.
void GenHistSink::inc(int slot, double val, double weight) {
    const Slot& sl = slots_[static_cast<std::size_t>(slot)];
    const int ibin = FortranBin(val, sl.ax);
    if (ibin < 1 || ibin > kNHbins) return;   // no under/overflow in the Fortran
    TH1D* h = sl.h;
    const double entries = h->GetEntries();
    const double err = h->GetBinError(ibin);   // x*x, not std::pow(x,2): gcc 16 calls libm pow, not bitwise x*x
    const double err2 = err * err + weight * weight;
    h->SetBinContent(ibin, h->GetBinContent(ibin) + weight);
    h->SetBinError(ibin, std::sqrt(err2));
    h->SetEntries(entries + 1);   // SetBinContent changes fEntries; keep "number of inc() calls counted"
}

void GenHistSink::fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                       bool success) {
    // FORTRAN-BUG(D-057 F18): original Fortran simc.f:259-261, 283-285 fill the P-arm histograms with vertex%p, which for
    // doing_rho is the rho (its lab theta/phi in the slope fields); corrected: the decay pion (orig%p).
    const ArmFull& hp = GeneratedHadronArm(rhoGenFromOrig_, vertex, orig);
    // simc.f:253-262 — every tried event
    inc(kGeni_e_delta, vertex.e.delta, 1.0);
    inc(kGeni_e_yptar, vertex.e.yptar, 1.0);
    inc(kGeni_e_xptar, -vertex.e.xptar, 1.0);
    inc(kGeni_p_delta, hp.delta, 1.0);
    inc(kGeni_p_yptar, hp.yptar, 1.0);
    inc(kGeni_p_xptar, -hp.xptar, 1.0);
    inc(kGeni_Em, vertex.Em, 1.0);
    inc(kGeni_Pm, vertex.Pm, 1.0);
    if (!success) return;
    // simc.f:272-286
    inc(kRecon_e_delta, main.RECON.e.delta, main.weight);
    inc(kRecon_e_yptar, main.RECON.e.yptar, main.weight);
    // FORTRAN-DIFF-CANDIDATE(D-011): simc.f:274/277 fill H%RECON xptar NOT negated, while H%gen/H%geni use -vertex%xptar
    // (simc.f:256/262), so the .gen recon and generated xptar distributions are mirror images. Intended convention?
    inc(kRecon_e_xptar, main.RECON.e.xptar, main.weight);
    inc(kRecon_p_delta, main.RECON.p.delta, main.weight);
    inc(kRecon_p_yptar, main.RECON.p.yptar, main.weight);
    inc(kRecon_p_xptar, main.RECON.p.xptar, main.weight);   // FORTRAN-DIFF-CANDIDATE(D-011): not negated, see above
    inc(kRecon_Em, recon.Em, 1.0);
    inc(kRecon_Pm, recon.Pm, 1.0);
    inc(kGen_e_delta, vertex.e.delta, 1.0);
    inc(kGen_e_yptar, vertex.e.yptar, 1.0);
    inc(kGen_e_xptar, -vertex.e.xptar, 1.0);
    inc(kGen_p_delta, hp.delta, 1.0);
    inc(kGen_p_yptar, hp.yptar, 1.0);
    inc(kGen_p_xptar, -hp.xptar, 1.0);
    inc(kGen_Em, vertex.Em, 1.0);
}

// Parallel event loop (IParallelSink): each worker fills its own histograms; the merge adds contents and errors in
// quadrature bin by bin (the same bins as the sequential run, the sums in another order).
std::unique_ptr<IOutputSink> GenHistSink::makeWorkerSink(int, int) const {
    ROOT::EnableThreadSafety();   // every worker fills its own TH1Ds (no directory) in its own thread
    return std::make_unique<GenHistSink>();
}

void GenHistSink::mergeWorkerSink(IOutputSink& sink) {
    auto& w = dynamic_cast<GenHistSink&>(sink);
    if (w.slots_.size() != slots_.size()) throw std::logic_error("GenHistSink::mergeWorkerSink: different bookings");
    for (std::size_t i = 0; i < slots_.size(); ++i) {
        TH1D* h = slots_[i].h;
        const TH1D* g = w.slots_[i].h;
        const double entries = h->GetEntries() + g->GetEntries();
        for (int b = 1; b <= kNHbins; ++b) {
            const double e1 = h->GetBinError(b), e2 = g->GetBinError(b);
            h->SetBinContent(b, h->GetBinContent(b) + g->GetBinContent(b));
            h->SetBinError(b, std::sqrt(e1 * e1 + e2 * e2));
        }
        h->SetEntries(entries);
    }
}

const TH1D* GenHistSink::hist(const std::string& name) const {
    auto it = h_.find(name);
    return it == h_.end() ? nullptr : it->second.get();
}

void GenHistSink::endRun(const RunSummary&) {
    std::unique_ptr<TFile> f(TFile::Open(file_.c_str(), "RECREATE"));
    if (!f || f->IsZombie()) throw std::runtime_error("GenHistSink: cannot create " + file_);
    for (auto& [n, h] : h_) h->Write(n.c_str());
    f->Close();
}

}  // namespace simc

SIMC_REGISTER(simc::IOutputSink, "gen_hist", simc::GenHistSink)
