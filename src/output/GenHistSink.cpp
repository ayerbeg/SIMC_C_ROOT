// src/output/GenHistSink.cpp
// Fortran: simc.f:inc (l.624-640), fills at simc.f l.254-286, axes init.f l.515-570,
// text output .gen at simc.f l.534-605.
#include "simc/output/GenHistSink.h"

#include <cmath>
#include <stdexcept>

#include <TFile.h>
#include <TH1D.h>

#include "simc/core/Event.h"
#include "simc/core/Registry.h"

namespace simc {

GenHistSink::GenHistSink() = default;
GenHistSink::~GenHistSink() = default;

void GenHistSink::configure(const Json& config, const SimcParameters&) {
    if (config.is_object()) file_ = config.value("file", file_);
}

void GenHistSink::init(const StageContext&) {}

int GenHistSink::FortranBin(double val, const HistAxis& ax) {
    // ibin = nint(0.5+(val-hist%min)/hist%bin); Fortran nint rounds half away from zero.
    const double x = 0.5 + (val - ax.min) / ax.bin;
    if (!(std::abs(x) < 1.0e9)) return 0;   // integer*4 range guard (also NaN)
    return static_cast<int>(std::lround(x));
}

void GenHistSink::book(const std::string& name, const HistAxis& ax) {
    auto h = std::make_unique<TH1D>(name.c_str(), name.c_str(), kNHbins, ax.min, ax.min + kNHbins * ax.bin);
    h->SetDirectory(nullptr);
    h->Sumw2();
    h_[name] = std::move(h);
    ax_[name] = ax;
}

void GenHistSink::beginRun(const RunSetup& setup) {
    setup_ = setup;
    if (file_.empty()) file_ = setup.output_dir + "/" + setup.run_name + ".gen.root";
    const auto& H = setup.hist;
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

void GenHistSink::inc(const std::string& name, double val, double weight) {
    const int ibin = FortranBin(val, ax_.at(name));
    if (ibin < 1 || ibin > kNHbins) return;   // no under/overflow in the Fortran
    TH1D* h = h_.at(name).get();
    const double entries = h->GetEntries();
    const double err2 = std::pow(h->GetBinError(ibin), 2) + weight * weight;
    h->SetBinContent(ibin, h->GetBinContent(ibin) + weight);
    h->SetBinError(ibin, std::sqrt(err2));
    h->SetEntries(entries + 1);   // SetBinContent changes fEntries; keep "number of inc() calls counted"
}

void GenHistSink::fill(const MainEvent& main, const Event& vertex, const Event& /*orig*/, const Event& recon,
                       bool success) {
    // simc.f:253-262 — every tried event
    inc("geni_e_delta", vertex.e.delta, 1.0);
    inc("geni_e_yptar", vertex.e.yptar, 1.0);
    inc("geni_e_xptar", -vertex.e.xptar, 1.0);
    inc("geni_p_delta", vertex.p.delta, 1.0);
    inc("geni_p_yptar", vertex.p.yptar, 1.0);
    inc("geni_p_xptar", -vertex.p.xptar, 1.0);
    inc("geni_Em", vertex.Em, 1.0);
    inc("geni_Pm", vertex.Pm, 1.0);
    if (!success) return;
    // simc.f:272-286
    inc("recon_e_delta", main.RECON.e.delta, main.weight);
    inc("recon_e_yptar", main.RECON.e.yptar, main.weight);
    // FORTRAN-DIFF-CANDIDATE(D-011): simc.f:274/277 fill H%RECON xptar NOT negated, while H%gen/H%geni use -vertex%xptar
    // (simc.f:256/262), so the .gen recon and generated xptar distributions are mirror images. Intended convention?
    inc("recon_e_xptar", main.RECON.e.xptar, main.weight);
    inc("recon_p_delta", main.RECON.p.delta, main.weight);
    inc("recon_p_yptar", main.RECON.p.yptar, main.weight);
    inc("recon_p_xptar", main.RECON.p.xptar, main.weight);   // FORTRAN-DIFF-CANDIDATE(D-011): not negated, see above
    inc("recon_Em", recon.Em, 1.0);
    inc("recon_Pm", recon.Pm, 1.0);
    inc("gen_e_delta", vertex.e.delta, 1.0);
    inc("gen_e_yptar", vertex.e.yptar, 1.0);
    inc("gen_e_xptar", -vertex.e.xptar, 1.0);
    inc("gen_p_delta", vertex.p.delta, 1.0);
    inc("gen_p_yptar", vertex.p.yptar, 1.0);
    inc("gen_p_xptar", -vertex.p.xptar, 1.0);
    inc("gen_Em", vertex.Em, 1.0);
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
