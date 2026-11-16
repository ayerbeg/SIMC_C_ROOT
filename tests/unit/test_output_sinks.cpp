// tests/unit/test_output_sinks.cpp
//   test_output_sinks <reference h10 .root> <docs/manual/markdown/reference/summary.schema.json>
// 1. The h10 branch list/order/types equal the Fortran make_root_tree tree.
// 2. Synthetic events through RootTreeSink / RunSummarySink / GenHistSink (created by name
//    from the Registry): column mapping (both arm sides), Float_t storage, summary JSON valid
//    against the schema, Fortran inc() binning.
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TTree.h>

#include "TestUtil.h"
#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/RunSummary.h"
#include "simc/output/GenHistSink.h"
#include "simc/output/NtupleLayout.h"
#include "simc/output/SummaryJson.h"
#include "simc/stages/IOutputSink.h"

using namespace simc;
using simctest::check;

namespace {

bool near(double a, double b, double rel = 1e-6) { return std::abs(a - b) <= rel * std::max(1.0, std::abs(b)); }

void synthetic(int i, MainEvent& m, Event& v, Event& o, Event& r) {
    const double x = 0.01 * i;
    v = Event{}; o = Event{}; r = Event{}; m = MainEvent{};
    v.e.delta = -5 + x; v.e.yptar = 0.01 + x / 100; v.e.xptar = 0.02 - x / 100; v.Em = 1 + x; v.Pm = 50 + x;
    v.p.delta = 3 - x; v.p.yptar = -0.01; v.p.xptar = 0.03;
    o = v;
    o.e.delta += 0.1;
    r = v;
    r.e.delta = 1.5 + x; r.e.yptar = 0.011; r.e.xptar = 0.021; r.e.z = 0.3; r.e.theta = 0.22; r.e.P = 8300;
    r.p.delta = -2.5; r.p.yptar = -0.012; r.p.xptar = 0.031; r.p.z = -0.4; r.p.theta = 0.64;
    r.q = 3000; r.nu = 2200; r.Q2 = 4.2e6; r.W = 940; r.epsilon = 0.95; r.Em = 12; r.Pm = 34;
    r.theta_pq = 0.01; r.phi_pq = 1.2; r.Pmx = 5; r.Pmy = 6; r.Pmz = 7; r.PmPar = 8; r.PmPer = 9; r.PmOop = 10;
    r.uq = {0.0, -0.6, 0.8};
    m.FP.e = {1.1, 0.002, 2.2, 0.003, 1800};
    m.FP.p = {-3.3, -0.004, 4.4, 0.005, 2000};
    m.target.z = 1.7; m.target.rastery = 0.05;
    m.ntup.radphot = 123.0; m.sigcc = 5.9e-4; m.weight = 3.2e-4 * (1 + x);
    m.RECON.e = {r.e.delta, r.e.yptar, r.e.xptar, r.e.z};
    m.RECON.p = {r.p.delta, r.p.yptar, r.p.xptar, r.p.z};
}

}  // namespace

int main(int argc, char** argv) {
    const std::string ref = argc > 1 ? argv[1] : "";
    const std::string schema_path = argc > 2 ? argv[2] : "";
    if (!simctest::exists(ref)) { std::printf("SKIP: reference tree %s missing\n", ref.c_str()); return simctest::kSkip; }

    // ---- 1. reference branch list
    std::vector<std::string> ref_names;
    {
        std::unique_ptr<TFile> f(TFile::Open(ref.c_str(), "READ"));
        TTree* t = f ? dynamic_cast<TTree*>(f->Get("h10")) : nullptr;
        check(t != nullptr, "reference file has tree h10");
        if (t) {
            for (auto* ob : *t->GetListOfBranches()) {
                auto* b = static_cast<TBranch*>(ob);
                ref_names.push_back(b->GetName());
                auto* leaf = static_cast<TLeaf*>(b->GetListOfLeaves()->At(0));
                check(std::string(leaf->GetTypeName()) == "Float_t", std::string("reference ") + b->GetName() + " is Float_t");
            }
        }
    }
    const auto& cols = EepNtupleColumns();
    check(ref_names.size() == cols.size() && cols.size() == 46, "46 columns in reference and layout");
    for (size_t i = 0; i < std::min(ref_names.size(), cols.size()); ++i)
        check(ref_names[i] == cols[i].name, "column " + std::to_string(i + 1) + ": layout " + cols[i].name + " vs Fortran " + ref_names[i]);

    // ---- 2. synthetic run
    std::system("mkdir -p output_sinks_tmp");
    SimcParameters par;
    par.electron_arm = 5; par.hadron_arm = 1; par.targ.A = 1; par.targ.Z = 1;
    par.Ebeam = 10600; par.spec.e.theta = 12.5; par.spec.p.theta = 36.0;
    par.SPedge.e.delta = {-10, 22}; par.SPedge.p.delta = {-10, 10};
    RunSetup setup = RunSetup::FromParameters(par, "synthetic");
    setup.output_dir = "output_sinks_tmp";
    setup.hist.e.delta = {-10, 0.5}; setup.hist.e.yptar = {-0.05, 0.002}; setup.hist.e.xptar = {-0.05, 0.002};
    setup.hist.p.delta = {-10, 0.5}; setup.hist.p.yptar = {-0.05, 0.002}; setup.hist.p.xptar = {-0.05, 0.002};
    setup.hist.Em = {0, 2.0}; setup.hist.Pm = {0, 2.0};
    check(setup.flags.doing_hyd_elast && UsesEepLayout(setup.flags), "H target gives the (e,e'p) layout");

    std::vector<std::unique_ptr<IOutputSink>> sinks;
    for (const char* n : {"root_tree", "run_summary", "gen_hist"}) {
        check(Registry<IOutputSink>::contains(n), std::string("sink registered: ") + n);
        sinks.push_back(Registry<IOutputSink>::create(n));
    }
    for (auto& s : sinks) { s->configure(Json::object(), par); s->init(setup); s->beginRun(setup); }

    const int ntried = 20;
    int nsucc = 0;
    MainEvent m; Event v, o, r;
    for (int i = 0; i < ntried; ++i) {
        synthetic(i, m, v, o, r);
        const bool ok = (i % 2 == 0);
        nsucc += ok;
        for (auto& s : sinks) s->fill(m, v, o, r, ok);
    }
    // GenHist content checks before endRun
    auto* gh = dynamic_cast<GenHistSink*>(sinks[2].get());
    check(gh && gh->hist("geni_e_delta")->GetEntries() == ntried, "geni filled for every try");
    check(gh && gh->hist("gen_e_delta")->GetEntries() == nsucc, "gen filled for successes");
    check(gh && gh->hist("gen_Pm")->GetEntries() == 0, "gen_Pm never filled (as the Fortran)");
    check(gh && gh->hist("geni_Pm")->GetEntries() == ntried, "geni_Pm filled with vertex Pm");
    check(gh && std::abs(gh->hist("recon_e_delta")->Integral() - [&] { double w = 0; for (int i = 0; i < ntried; i += 2) w += 3.2e-4 * (1 + 0.01 * i); return w; }()) < 1e-12,
          "recon histograms weighted by main.weight");
    {
        HistAxis ax{-10, 0.5};
        check(GenHistSink::FortranBin(-10.0, ax) == 1, "inc(): val=min -> bin 1");
        check(GenHistSink::FortranBin(-9.5000001, ax) == 1, "inc(): just below upper edge -> bin 1");
        check(GenHistSink::FortranBin(-9.5, ax) == 2, "inc(): upper edge -> bin 2");
        check(GenHistSink::FortranBin(-10.0000001, ax) == 0, "inc(): just below min -> dropped");
        check(GenHistSink::FortranBin(15.0, ax) == 51, "inc(): max -> 51 (dropped)");
    }

    RunSummary sum;
    sum.run_name = "synthetic"; sum.ngen = nsucc; sum.ntried = ntried; sum.ncontribute = nsucc; sum.nevent = nsucc;
    sum.wtcontribute = 42.0; sum.sum_sigcc = 1.0e-3; sum.central_sigcc = 5.8e-4; sum.charge = 1.0;
    sum.targetfac = 2.67954e-9; sum.luminosity = 3.73199e8; sum.genvol = 0.0224744; sum.normfac = 4.04142e5;
    StopCounters hms{"hms", "hadron", {{"trials", ntried}, {"slit_hor", 7}, {"successes", nsucc}}};
    StopCounters shms{"shms", "electron", {{"trials", nsucc}, {"successes", nsucc}}};
    sum.stops = {hms, shms};
    for (auto& s : sinks) s->endRun(sum);

    // Tree content
    {
        std::unique_ptr<TFile> f(TFile::Open("output_sinks_tmp/synthetic.root", "READ"));
        TTree* t = f ? dynamic_cast<TTree*>(f->Get("h10")) : nullptr;
        check(t && t->GetEntries() == nsucc, "h10 has one entry per success");
        if (t) {
            size_t k = 0;
            for (auto* ob : *t->GetListOfBranches()) {
                auto* b = static_cast<TBranch*>(ob);
                check(k < cols.size() && cols[k].name == std::string(b->GetName()), std::string("written branch order ") + b->GetName());
                check(std::string(static_cast<TLeaf*>(b->GetListOfLeaves()->At(0))->GetTypeName()) == "Float_t",
                      std::string("written ") + b->GetName() + " is Float_t");
                ++k;
            }
            float val[46];
            for (size_t c = 0; c < cols.size(); ++c) t->SetBranchAddress(cols[c].name, &val[c]);
            t->GetEntry(1);   // second success = synthetic(2)
            synthetic(2, m, v, o, r);
            // electron_arm=5 (left): hs* = hadron (HMS), ss* = electron (results_write.f:89-117)
            check(near(val[0], r.p.delta), "hsdelta = recon.p.delta");
            check(near(val[4], m.FP.p.x), "hsxfp = FP.p.x");
            check(near(val[8], o.p.delta), "hsdeltai = orig.p.delta");
            check(near(val[11], m.target.z * setup.p_sin_th), "hsytari = z*sin(theta_p)");
            check(near(val[12], r.e.delta), "ssdelta = recon.e.delta");
            check(near(val[20], o.e.delta), "ssdeltai = orig.e.delta");
            check(near(val[23], -m.target.z * setup.e_sin_th), "ssytari = -z*sin(theta_e)");
            check(near(val[24], r.q / 1000.), "q in GeV");
            check(near(val[26], r.Q2 / 1e6), "Q2 in GeV^2");
            const double poft = 938.27231 * 10600. / (2 * 10600. * (std::sin(r.e.theta / 2) * std::sin(r.e.theta / 2)) + 938.27231);
            check(near(val[33], (r.e.P - poft) / 1000., 1e-5), "corrsing");
            check(near(val[34], -r.Pmx / 1000.), "Pmx = -recon.Pmx");
            check(near(val[35], (r.Pmz * r.uq.y - r.Pmy * r.uq.z) / 1.0 / 1000.), "Pmy (Heep)");
            check(near(val[36], -(r.Pmy * r.uq.y + r.Pmz * r.uq.z) / 1.0 / 1000.), "Pmz (Heep)");
            check(near(val[40], -m.target.rastery), "fry = -rastery");
            check(near(val[41], m.ntup.radphot / 1000.), "radphot GeV");
            check(near(val[43], m.weight), "Weight");
            check(near(val[44], r.e.theta) && near(val[45], r.p.theta), "theta_e, theta_p");
        }
    }
    // Opposite arm assignment (electron_arm = 1: electron on the right)
    {
        RunSetup s2 = setup;
        s2.params.electron_arm = 1;
        std::array<double, kEepNtupleSize> ntu{};
        synthetic(3, m, v, o, r);
        FillEepNtuple(s2, m, v, o, r, ntu);
        check(near(ntu[0], r.e.delta) && near(ntu[12], r.p.delta), "electron_arm=1: hs* = electron, ss* = hadron");
        check(near(ntu[11], m.target.z * s2.e_sin_th) && near(ntu[23], -m.target.z * s2.p_sin_th), "ytari signs, electron right");
        s2.params.electron_arm = 9;   // unknown: columns 1-24 keep previous content (static array)
        std::array<double, kEepNtupleSize> keep = ntu;
        synthetic(4, m, v, o, r);
        check(!FillEepNtuple(s2, m, v, o, r, ntu), "unknown electron_arm reported");
        check(ntu[0] == keep[0] && ntu[23] == keep[23], "unknown arm leaves columns 1-24 untouched");
    }

    // Summary JSON
    {
        std::ifstream fs("output_sinks_tmp/synthetic.summary.json");
        check(fs.good(), "summary.json written");
        std::ifstream ft("output_sinks_tmp/synthetic.summary.txt");
        check(ft.good(), "summary.txt written");
        Json j = Json::parse(fs);
        Json schema;
        std::ifstream sc(schema_path);
        check(sc.good(), "schema file present");
        if (sc.good()) {
            schema = Json::parse(sc);
            auto errs = ValidateAgainstSchema(j, schema);
            for (const auto& e : errs) std::printf("  schema: %s\n", e.c_str());
            check(errs.empty(), "summary.json validates against docs/manual/markdown/reference/summary.schema.json");
            Json bad = j;
            bad.erase("normalization");
            bad["counters"]["ntried"] = "x";
            check(ValidateAgainstSchema(bad, schema).size() >= 2, "validator rejects a broken summary");
        }
        check(near(j["normalization"]["wtcontr"].get<double>(), 42.0 / nsucc), "wtcontr = wtcontribute/nevent");
        check(near(j["normalization"]["sigcc_average"].get<double>(), 1.0e-3 / nsucc), "average_sigcc = sum_sigcc/nevent");
        check(near(j["normalization"]["luminosity_GeV2"].get<double>(), 3.73199e8 * ((197.327053 / 1e5) * (197.327053 / 1e5))),
              "luminosity in GeV^2");
        check(j["stops"]["hms"]["slit_hor"] == 7, "stop counters by aperture name");
        check(j["stop_roles"]["hms"] == "hadron", "stop role");
        check(j["run"]["reaction"] == "H(e,e'p)" && j["run"]["producer"] == "simc++", "reaction title, producer");
    }
    // GenHist file
    {
        std::unique_ptr<TFile> f(TFile::Open("output_sinks_tmp/synthetic.gen.root", "READ"));
        check(f && f->Get("recon_p_xptar") && f->Get("geni_Em"), "gen.root has the histograms");
    }
    return simctest::finish("output_sinks");
}
