// src/driver/ReportSink.cpp — see include/simc/driver/ReportSink.h.
#include "simc/driver/ReportSink.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/driver/FortranWriter.h"

namespace simc {
namespace driver {

void FortranHistSink::configure(const Json& config, const SimcParameters&) {
    if (config.is_object()) {
        hist_ = config.value("hist", hist_);
        geniFile_ = config.value("geni", geniFile_);
        genFile_ = config.value("gen", genFile_);
    }
}

void FortranHistSink::beginRun(const RunSetup& setup) {
    if (hist_.empty()) hist_ = setup.output_dir + "/" + setup.run_name + ".hist";
    if (geniFile_.empty()) geniFile_ = setup.output_dir + "/" + setup.run_name + ".geni";
    if (genFile_.empty()) genFile_ = setup.output_dir + "/" + setup.run_name + ".gen";
    // init.f:limits_init l.515-570: H%gen, H%geni and H%RECON share the axes
    const DiagHistAxes& H = setup.hist;
    for (HistSet* s : {&RECON_, &gen_, &geni_}) {
        *s = HistSet{};
        s->e.delta.ax = H.e.delta;
        s->e.yptar.ax = H.e.yptar;
        s->e.xptar.ax = H.e.xptar;
        s->p.delta.ax = H.p.delta;
        s->p.yptar.ax = H.p.yptar;
        s->p.xptar.ax = H.p.xptar;
        s->Em.ax = H.Em;
        s->Pm.ax = H.Pm;
    }
}

// Fortran: simc.f:inc (l.624-640). ibin = nint(0.5+(val-min)/bin) (integer*4; values far outside never land in 1..50).
void FortranHistSink::Inc(Hist& h, double val, double weight) {
    const double x = 0.5 + (val - h.ax.min) / h.ax.bin;
    if (!(std::fabs(x) < 1.0e9)) return;
    const long ibin = std::lround(x);   // nint: half away from zero
    if (ibin >= 1 && ibin <= kNHbins) h.buf[static_cast<std::size_t>(ibin)] = h.buf[static_cast<std::size_t>(ibin)] + weight;
}

// Fortran: simc.f l.250-287
void FortranHistSink::fill(const MainEvent& main, const Event& vertex, const Event&, const Event& recon, bool success) {
    const double one = 1.0;
    Inc(geni_.e.delta, vertex.e.delta, one);
    Inc(geni_.e.yptar, vertex.e.yptar, one);
    Inc(geni_.e.xptar, -vertex.e.xptar, one);
    Inc(geni_.p.delta, vertex.p.delta, one);
    Inc(geni_.p.yptar, vertex.p.yptar, one);
    Inc(geni_.p.xptar, -vertex.p.xptar, one);
    Inc(geni_.Em, vertex.Em, one);
    Inc(geni_.Pm, vertex.Pm, one);
    if (!success) return;
    Inc(RECON_.e.delta, main.RECON.e.delta, main.weight);
    Inc(RECON_.e.yptar, main.RECON.e.yptar, main.weight);
    Inc(RECON_.e.xptar, main.RECON.e.xptar, main.weight);   // FORTRAN-DIFF-CANDIDATE(D-011): not negated (simc.f:274)
    Inc(RECON_.p.delta, main.RECON.p.delta, main.weight);
    Inc(RECON_.p.yptar, main.RECON.p.yptar, main.weight);
    Inc(RECON_.p.xptar, main.RECON.p.xptar, main.weight);   // FORTRAN-DIFF-CANDIDATE(D-011): not negated (simc.f:277)
    Inc(RECON_.Em, recon.Em, one);
    Inc(RECON_.Pm, recon.Pm, one);
    Inc(gen_.e.delta, vertex.e.delta, one);
    Inc(gen_.e.yptar, vertex.e.yptar, one);
    Inc(gen_.e.xptar, -vertex.e.xptar, one);
    Inc(gen_.p.delta, vertex.p.delta, one);
    Inc(gen_.p.yptar, vertex.p.yptar, one);
    Inc(gen_.p.xptar, -vertex.p.xptar, one);
    Inc(gen_.Em, vertex.Em, one);
}

namespace {
double centre(const FortranHistSink::Hist& h, int i) { return h.ax.min + (i - 0.5) * h.ax.bin; }
}  // namespace

// Fortran: simc.f l.534-612
std::string FortranHistSink::formatGen() const {
    using fortran::Record;
    std::ostringstream os;
    auto header = [&](const char* title, std::initializer_list<const char*> cols) {
        os << fortran::List(title).str() << "\n";
        Record r;
        for (const char* c : cols) r.a(c, 12);
        os << r.str() << "\n";
    };
    // '(3(1x,2(e11.4,1x)))'
    auto rows3 = [&](const Hist& a, const Hist& ab, const Hist& b, const Hist& bb, const Hist& c, const Hist& cb) {
        for (int i = 1; i <= kNHbins; ++i) {
            Record r;
            r.x(1).e(centre(a, i), 11, 4).x(1).e(ab.buf[i], 11, 4).x(1);
            r.x(1).e(centre(b, i), 11, 4).x(1).e(bb.buf[i], 11, 4).x(1);
            r.x(1).e(centre(c, i), 11, 4).x(1).e(cb.buf[i], 11, 4).x(1);
            os << r.str() << "\n";
        }
    };
    header("E arm Experimental Target Distributions:", {"delta", "EXPERIM", "yptar", "EXPERIM", "xptar", "EXPERIM"});
    rows3(RECON_.e.delta, RECON_.e.delta, RECON_.e.yptar, RECON_.e.yptar, RECON_.e.xptar, RECON_.e.xptar);
    header("P arm Experimental Target Distributions:", {"delta", "EXPERIM", "yptar", "EXPERIM", "xptar", "EXPERIM"});
    rows3(RECON_.p.delta, RECON_.p.delta, RECON_.p.yptar, RECON_.p.yptar, RECON_.p.xptar, RECON_.p.xptar);
    header("Distributions of Contributing E arm Events:", {"delta", "CONTRIB", "yuptar", "CONTRIB", "xptar", "CONTRIB"});
    rows3(gen_.e.delta, gen_.e.delta, gen_.e.yptar, gen_.e.yptar, gen_.e.xptar, gen_.e.xptar);
    header("Distributions of Contributing P arm Events:", {"delta", "CONTRIB", "yptar", "CONTRIB", "xptar", "CONTRIB"});
    rows3(gen_.p.delta, gen_.p.delta, gen_.p.yptar, gen_.p.yptar, gen_.p.xptar, gen_.p.xptar);
    // FORTRAN-QUIRK simc.f:578-586: H%gen delta, then H%geni yptar/xptar, on the H%gen bin centres
    header("Original E arm Events:", {"delta", "ORIGIN", "yptar", "ORIGIN", "xptar", "ORIGIN"});
    rows3(gen_.e.delta, gen_.e.delta, gen_.e.yptar, geni_.e.yptar, gen_.e.xptar, geni_.e.xptar);
    header("Original P arm Events:", {"delta", "ORIGIN", "yptar", "ORIGIN", "xptar", "ORIGIN"});
    rows3(geni_.p.delta, geni_.p.delta, geni_.p.yptar, geni_.p.yptar, geni_.p.xptar, geni_.p.xptar);
    // '(4a12)' and '(3(1x,4(e11.4,1x)))' with 4 items
    header("Original Em/Pm distributions:", {"Em", "ORIGIN", "Pm", "ORIGIN"});
    for (int i = 1; i <= kNHbins; ++i) {
        Record r;
        r.x(1).e(centre(geni_.Em, i), 11, 4).x(1).e(geni_.Em.buf[i], 11, 4).x(1);
        r.e(centre(geni_.Pm, i), 11, 4).x(1).e(geni_.Pm.buf[i], 11, 4).x(1);
        os << r.str() << "\n";
    }
    return os.str();
}

void FortranHistSink::writeReport(const RunSetup& setup, const RunReport& report) {
    // simc.f order: .geni (l.462), .gen (l.534), .hist (l.616)
    const std::string geni = FormatGeni(setup, report);
    const std::string hist = FormatHist(setup, report);
    std::ofstream g(geniFile_);
    if (!g) throw std::runtime_error("fortran_hist: cannot write " + geniFile_);
    g << geni;
    std::ofstream gn(genFile_);
    if (!gn) throw std::runtime_error("fortran_hist: cannot write " + genFile_);
    gn << formatGen();
    std::ofstream h(hist_);
    if (!h) throw std::runtime_error("fortran_hist: cannot write " + hist_);
    h << hist;
}

}  // namespace driver
}  // namespace simc

SIMC_REGISTER(simc::IOutputSink, "fortran_hist", simc::driver::FortranHistSink)
