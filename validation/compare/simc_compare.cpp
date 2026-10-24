// simc_compare - decide whether two SIMC runs (Fortran or C++) are statistically
// the same. Task P1.5. Formulas and formats: validation/README.md.
//
// Usage:
//   simc_compare --a SPEC --b SPEC --out DIR [--config CFG.json]
//                [--label-a NAME] [--label-b NAME] [--perturb VAR:SHIFT]... [--empty-tree-b]
//                [--skip v1,v2] [--nbins N] [--prefer-hist] [--no-pdf]
//                [--expect PASS|FAIL|IDENTICAL] [--expect-flagged v1,v2]
//   simc_compare --dump-summary SPEC OUT.json
// SPEC: <run>.summary.json | <run>.hist | <run>.root | directory (see RunSummary.h)
// Exit: 0 = PASS (or expectation met), 1 = FAIL (or expectation not met), 2 = error.

#include "RunSummary.h"
#include "Stats.h"

#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TText.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace simcval;
namespace fs = std::filesystem;
using json = nlohmann::json;

namespace {

// ------------------------------------------------------------------ config
struct VarBinning {
    int nbins = 0;
    double lo = 0, hi = 0;
    bool fixed = false;
};

struct Config {
    double alpha = 0.01;          // per-test p threshold
    double globalAlpha = 0.001;   // binomial tail probability of the failure count
    double grossP = 1e-5;         // any counted test below this -> FAIL
    int nbins = 50;
    double qlo = 0.005, qhi = 0.995;
    double detRelTol = 1e-6;      // deterministic scalars (plus print precision)
    std::vector<std::string> variables;  // empty = all common branches
    std::set<std::string> skip;
    std::map<std::string, VarBinning> binning;
    std::string normA = "per_tried", normB = "per_tried";
    std::string treeName;         // override
    std::string weightBranch;     // override
    std::string sigccBranch = "sigcc";
    bool testChi2Abs = true, testChi2WW = true, testChi2UU = true, testKS = true;
    bool stopPerCounter = true;
    double ratioLo = 0.8, ratioHi = 1.2;
    bool pdf = true;
};

void loadConfig(const std::string& path, Config& c) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open config " + path);
    json j;
    in >> j;
    auto num = [&](const char* k, double& d) { if (j.contains(k)) d = j[k].get<double>(); };
    num("alpha", c.alpha);
    num("global_alpha", c.globalAlpha);
    num("gross_p", c.grossP);
    num("quantile_lo", c.qlo);
    num("quantile_hi", c.qhi);
    num("det_rel_tol", c.detRelTol);
    num("ratio_lo", c.ratioLo);
    num("ratio_hi", c.ratioHi);
    if (j.contains("nbins")) c.nbins = j["nbins"].get<int>();
    if (j.contains("variables")) c.variables = j["variables"].get<std::vector<std::string>>();
    if (j.contains("skip")) for (auto& s : j["skip"]) c.skip.insert(s.get<std::string>());
    if (j.contains("binning"))
        for (auto it = j["binning"].begin(); it != j["binning"].end(); ++it) {
            VarBinning b;
            b.nbins = it->value("nbins", c.nbins);
            if (it->contains("min") && it->contains("max")) {
                b.lo = (*it)["min"].get<double>();
                b.hi = (*it)["max"].get<double>();
                b.fixed = true;
            }
            c.binning[it.key()] = b;
        }
    auto normStr = [](const json& v) {
        return v.is_number() ? std::to_string(v.get<double>()) : v.get<std::string>();
    };
    if (j.contains("norm_a")) c.normA = normStr(j["norm_a"]);
    if (j.contains("norm_b")) c.normB = normStr(j["norm_b"]);
    if (j.contains("tree_name")) c.treeName = j["tree_name"];
    if (j.contains("weight_branch")) c.weightBranch = j["weight_branch"];
    if (j.contains("sigcc_branch")) c.sigccBranch = j["sigcc_branch"];
    if (j.contains("stop_per_counter")) c.stopPerCounter = j["stop_per_counter"];
    if (j.contains("pdf")) c.pdf = j["pdf"];
    if (j.contains("tests")) {
        auto& t = j["tests"];
        c.testChi2Abs = t.value("chi2_abs", c.testChi2Abs);
        c.testChi2WW = t.value("chi2_ww", c.testChi2WW);
        c.testChi2UU = t.value("chi2_uu", c.testChi2UU);
        c.testKS = t.value("ks_weighted", c.testKS);
    }
}

// ------------------------------------------------------------------ tree
struct TreeData {
    std::string file, name;
    std::vector<std::string> names;
    std::map<std::string, std::vector<double>> cols;
    size_t n = 0;
    const std::vector<double>* col(const std::string& k) const {
        auto it = cols.find(k);
        return it == cols.end() ? nullptr : &it->second;
    }
};

TreeData loadTree(const std::string& file, const std::string& name) {
    TreeData d;
    d.file = file;
    d.name = name;
    std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
    if (!f || f->IsZombie()) throw std::runtime_error("cannot open ROOT file " + file);
    TTree* t = dynamic_cast<TTree*>(f->Get(name.c_str()));
    if (!t) throw std::runtime_error("no tree '" + name + "' in " + file);
    d.n = static_cast<size_t>(t->GetEntries());
    std::vector<TLeaf*> leaves;
    TObjArray* ll = t->GetListOfLeaves();
    for (int i = 0; i < ll->GetEntriesFast(); ++i) {
        auto* lf = static_cast<TLeaf*>(ll->At(i));
        if (lf->GetLen() != 1) continue;  // scalar leaves only
        leaves.push_back(lf);
        d.names.push_back(lf->GetName());
        d.cols[lf->GetName()].resize(d.n);
    }
    std::vector<std::vector<double>*> dst;
    for (auto& nm : d.names) dst.push_back(&d.cols[nm]);
    for (size_t e = 0; e < d.n; ++e) {
        t->GetEntry(static_cast<Long64_t>(e));
        for (size_t k = 0; k < leaves.size(); ++k) (*dst[k])[e] = leaves[k]->GetValue(0);
    }
    return d;
}

// ------------------------------------------------------------------ results
struct Test {
    std::string category;  // "norm", "stop", "det", "info", "tree"
    std::string name;      // quantity / variable
    std::string method;    // binomial, contingency, pull, chi2_abs, chi2_ww, chi2_uu, ks_w, det
    double a = NAN, b = NAN, ea = NAN, eb = NAN;
    double stat = NAN;
    int ndf = -1;
    double pull = NAN;
    double p = 1;
    bool counted = true;
    bool pass = true;
    std::string note;
};

std::string fmt(double v, int prec = 6) {
    if (std::isnan(v)) return "-";
    std::ostringstream o;
    o << std::setprecision(prec) << v;
    return o.str();
}

// Per-run normalisation of the event weight for "yield per tried event" plots/tests.
struct Norm {
    double scale = 1;     // multiplies Weight
    bool perTried = true; // variance formula (S2 - S1^2/N)/N^2 applies
    double ntried = 0;
    std::string desc;
};

Norm makeNorm(const std::string& mode, const RunSummary& s, size_t entries,
              std::vector<std::string>& warn, const std::string& which) {
    Norm n;
    n.ntried = s.has_counters ? double(s.ntried) : double(entries);
    if (mode == "per_tried") {
        if (!s.has_counters) warn.push_back(which + ": no Ntried, using tree entries for per-tried normalisation");
        n.scale = 1.0 / n.ntried;
        n.desc = "Weight/Ntried";
    } else if (mode == "normfac") {
        const Scalar* nf = s.scalar("normfac");
        if (!nf || !s.has_counters) throw std::runtime_error(which + ": norm 'normfac' needs normfac and counters");
        n.scale = nf->value / double(s.nevent());
        n.desc = "Weight*normfac/nevent";
    } else {
        n.scale = std::stod(mode);
        n.perTried = false;
        n.desc = "Weight*" + mode;
    }
    return n;
}

void quantileRange(const std::vector<double>& a, const std::vector<double>& b, double qlo,
                   double qhi, double& lo, double& hi) {
    std::vector<double> v;
    v.reserve(a.size() + b.size());
    for (double x : a) if (std::isfinite(x)) v.push_back(x);
    for (double x : b) if (std::isfinite(x)) v.push_back(x);
    if (v.empty()) { lo = hi = 0; return; }
    auto q = [&](double f) {
        size_t k = static_cast<size_t>(std::floor(f * (v.size() - 1)));
        std::nth_element(v.begin(), v.begin() + k, v.end());
        return v[k];
    };
    lo = q(qlo);
    hi = q(qhi);
    if (!(hi > lo)) {
        auto mm = std::minmax_element(v.begin(), v.end());
        lo = *mm.first;
        hi = *mm.second;
    }
    if (hi > lo) {  // make the upper edge inclusive of the quantile value
        double pad = 1e-9 * (hi - lo);
        hi += pad;
    }
}

struct VarResult {
    std::string name;
    std::string status = "compared";  // compared | constant | skipped | missing | empty | nonfinite
    double lo = 0, hi = 0;
    int nbins = 0;
    std::vector<Test*> tests;
    double minP = 1;
    std::unique_ptr<TH1D> hwA, hwB, huA, huB;
};

// ------------------------------------------------------------------ text
void writeLines(TCanvas& c, const std::vector<std::string>& lines, const std::string& pdf,
                bool& first) {
    const size_t perPage = 50;
    for (size_t s = 0; s < lines.size(); s += perPage) {
        c.Clear();
        c.cd();
        TText t;
        t.SetTextFont(82);  // monospace
        t.SetTextSize(0.0155);
        double y = 0.975;
        for (size_t i = s; i < std::min(lines.size(), s + perPage); ++i) {
            t.DrawTextNDC(0.01, y, lines[i].c_str());
            y -= 0.0192;
        }
        c.Print((pdf + (first ? "(" : "")).c_str(), "pdf Landscape");
        first = false;
    }
}

void usage() {
    std::cerr << "usage: simc_compare --a SPEC --b SPEC --out DIR [--config CFG.json]\n"
                 "                    [--label-a NAME] [--label-b NAME] [--perturb VAR:SHIFT]... [--empty-tree-b]\n"
                 "                    [--skip v1,v2] [--nbins N] [--prefer-hist] [--no-pdf]\n"
                 "                    [--expect PASS|FAIL|IDENTICAL] [--expect-flagged v1,v2]\n"
                 "                    [--mode-a compat|corrected] [--mode-b compat|corrected]\n"
                 "       simc_compare --dump-summary SPEC OUT.json\n";
}

std::vector<std::string> splitComma(const std::string& s) {
    std::vector<std::string> v;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, ',')) if (!t.empty()) v.push_back(t);
    return v;
}

}  // namespace

// ============================================================================
int main(int argc, char** argv) {
    std::string specA, specB, outDir, cfgFile, labelA = "A", labelB = "B", expect, dumpSpec, dumpOut;
    std::vector<std::pair<std::string, double>> perturb;
    std::vector<std::string> expectFlagged, extraSkip;
    int nbinsCli = 0;
    bool preferHist = false, noPdf = false, emptyTreeB = false;
    int modeA = -2, modeB = -2;   // --mode-a/--mode-b: 1 compat, 0 corrected, -2 = from the run summary
    try {
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            auto next = [&]() -> std::string {
                if (i + 1 >= argc) throw std::runtime_error("missing value after " + a);
                return argv[++i];
            };
            if (a == "--a") specA = next();
            else if (a == "--b") specB = next();
            else if (a == "--out") outDir = next();
            else if (a == "--config") cfgFile = next();
            else if (a == "--label-a") labelA = next();
            else if (a == "--label-b") labelB = next();
            else if (a == "--expect") expect = next();
            else if (a == "--expect-flagged") expectFlagged = splitComma(next());
            else if (a == "--skip") extraSkip = splitComma(next());
            else if (a == "--nbins") nbinsCli = std::stoi(next());
            else if (a == "--prefer-hist") preferHist = true;
            else if (a == "--no-pdf") noPdf = true;
            else if (a == "--empty-tree-b") emptyTreeB = true;
            else if (a == "--mode-a" || a == "--mode-b") {
                const std::string v = next();
                if (v != "compat" && v != "corrected") throw std::runtime_error(a + " needs compat or corrected");
                (a == "--mode-a" ? modeA : modeB) = v == "compat" ? 1 : 0;
            }
            else if (a == "--perturb") {
                std::string v = next();
                size_t c = v.rfind(':');
                if (c == std::string::npos) throw std::runtime_error("--perturb needs VAR:SHIFT");
                perturb.emplace_back(v.substr(0, c), std::stod(v.substr(c + 1)));
            } else if (a == "--dump-summary") {
                dumpSpec = next();
                dumpOut = next();
            } else if (a == "-h" || a == "--help") {
                usage();
                return 0;
            } else throw std::runtime_error("unknown argument " + a);
        }
    } catch (const std::exception& e) {
        std::cerr << "simc_compare: " << e.what() << "\n";
        usage();
        return 2;
    }

    gErrorIgnoreLevel = kError;

    if (!dumpSpec.empty()) {
        try {
            RunSummary s = loadRun(dumpSpec, !preferHist);
            writeSummaryJson(s, dumpOut);
            for (auto& w : s.warnings) std::cerr << "warning: " << w << "\n";
            std::cout << "wrote " << dumpOut << "\n";
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "simc_compare: " << e.what() << "\n";
            return 2;
        }
    }
    if (specA.empty() || specB.empty() || outDir.empty()) {
        usage();
        return 2;
    }

    Config cfg;
    std::vector<std::string> warnings;
    RunSummary RA, RB;
    TreeData TA, TB;
    try {
        if (!cfgFile.empty()) loadConfig(cfgFile, cfg);
        if (nbinsCli > 0) cfg.nbins = nbinsCli;
        for (auto& s : extraSkip) cfg.skip.insert(s);
        if (noPdf) cfg.pdf = false;
        fs::create_directories(outDir);
        RA = loadRun(specA, !preferHist);
        RB = loadRun(specB, !preferHist);
        if (modeA != -2) RA.fortran_compat = modeA;
        if (modeB != -2) RB.fortran_compat = modeB;
        for (auto& w : RA.warnings) warnings.push_back("A: " + w);
        for (auto& w : RB.warnings) warnings.push_back("B: " + w);
        if (!cfg.treeName.empty()) RA.treeName = RB.treeName = cfg.treeName;
        if (!cfg.weightBranch.empty()) RA.weightBranch = RB.weightBranch = cfg.weightBranch;
        if (RA.treeFile.empty() || RB.treeFile.empty()) throw std::runtime_error("tree file not found for a run");
        TA = loadTree(RA.treeFile, RA.treeName);
        TB = loadTree(RB.treeFile, RB.treeName);
    } catch (const std::exception& e) {
        std::cerr << "simc_compare: " << e.what() << "\n";
        return 2;
    }
    // empty tree (test option): run B keeps its branches but loses every entry
    if (emptyTreeB) {
        for (auto& kv : TB.cols) kv.second.clear();
        TB.n = 0;
        warnings.push_back("B: TREE EMPTIED (test option --empty-tree-b)");
    }
    // perturbation (test option): applied to run B
    for (auto& [v, sh] : perturb) {
        auto it = TB.cols.find(v);
        if (it == TB.cols.end()) {
            std::cerr << "simc_compare: --perturb: no branch " << v << " in B\n";
            return 2;
        }
        for (auto& x : it->second) x += sh;
        warnings.push_back("B: PERTURBED " + v + " by " + fmt(sh) + " (test option)");
    }
    try {
        writeSummaryJson(RA, (fs::path(outDir) / "runA.summary.json").string());
        writeSummaryJson(RB, (fs::path(outDir) / "runB.summary.json").string());
    } catch (...) {}

    std::vector<std::unique_ptr<Test>> tests;
    auto add = [&](Test t) {
        tests.push_back(std::make_unique<Test>(std::move(t)));
        return tests.back().get();
    };

    // ----------------------------------------------------------- run level
    const bool haveCounters = RA.has_counters && RB.has_counters;
    Norm NA, NB;
    try {
        NA = makeNorm(cfg.normA, RA, TA.n, warnings, "A");
        NB = makeNorm(cfg.normB, RB, TB.n, warnings, "B");
    } catch (const std::exception& e) {
        std::cerr << "simc_compare: " << e.what() << "\n";
        return 2;
    }
    const std::vector<double>* WA = TA.col(RA.weightBranch);
    const std::vector<double>* WB = TB.col(RB.weightBranch);
    if (!WA || !WB) {
        std::cerr << "simc_compare: weight branch missing (" << RA.weightBranch << ")\n";
        return 2;
    }

    if (haveCounters) {
        if (RA.ncontribute != int64_t(TA.n)) warnings.push_back("A: tree entries " + std::to_string(TA.n) + " != Ncontribute " + std::to_string(RA.ncontribute));
        if (RB.ncontribute != int64_t(TB.n)) warnings.push_back("B: tree entries " + std::to_string(TB.n) + " != Ncontribute " + std::to_string(RB.ncontribute));
        Test t;
        t.category = "norm";
        t.name = "acceptance Ncontribute/Ntried";
        t.method = "binomial";
        binomialPull(RA.ncontribute, RA.ntried, RB.ncontribute, RB.ntried, t.a, t.b, t.ea, t.eb, t.pull);
        t.p = pullPValue(t.pull);
        add(t);
        Test u;
        u.category = "norm";
        u.name = "Nco_no_rad_prot/Ncontribute";
        u.method = "binomial";
        if (binomialPull(RA.nco_no_rad_prot, RA.ncontribute, RB.nco_no_rad_prot, RB.ncontribute, u.a, u.b, u.ea, u.eb, u.pull)) {
            u.p = pullPValue(u.pull);
            add(u);
        }
    } else {
        warnings.push_back("run-level counters missing in A or B: normalisation pulls skipped");
    }

    // Tree yield per tried event: Y = scale*sum W ; see README
    auto treeYield = [](const std::vector<double>& W, const Norm& n, double& Y, double& sY) {
        double S1 = 0, S2 = 0;
        for (double w : W) if (std::isfinite(w)) { S1 += w; S2 += w * w; }
        if (n.perTried) {
            perTriedMean(S1, S2, n.ntried, Y, sY);
            double f = n.scale * n.ntried;  // 1 for per_tried, lumi*genvol for normfac
            Y *= f;
            sY *= f;
        } else {
            Y = n.scale * S1;
            sY = n.scale * std::sqrt(S2);
        }
    };
    double YA, sYA, YB, sYB;
    treeYield(*WA, NA, YA, sYA);
    treeYield(*WB, NB, YB, sYB);
    {
        Test t;
        t.category = "norm";
        t.name = "yield sum(" + RA.weightBranch + ") [" + NA.desc + "]";
        t.method = "pull";
        t.a = YA; t.ea = sYA; t.b = YB; t.eb = sYB;
        t.pull = pull(YA, sYA, YB, sYB);
        t.p = pullPValue(t.pull);
        add(t);
    }
    // wtcontr (.hist): sum of main%weight over the pass_cuts events only (simc.f:229-241, 311-313), times normfac/nevent.
    // P12.3 (D-053): its relative error is that of the pass_cuts subset of the tree, sqrt(sum w^2)/sum w of the events in
    // the run's own pass_cuts window (perTriedMean over Ntried as the tree yield); without the window data (SPedge, slop,
    // electron_arm, the mode) the relative error of the whole tree yield is used, with a note.
    if (const Scalar *wa = RA.scalar("wtcontr"), *wb = RB.scalar("wtcontr"); wa && wb) {
        Test t;
        t.category = "norm";
        t.name = "wtcontr";
        t.method = "pull";
        t.a = wa->value; t.b = wb->value;
        std::string why[2];
        auto relerr = [&](int k) -> double {
            const RunSummary& R = k ? RB : RA;
            const TreeData& T = k ? TB : TA;
            const Norm& N = k ? NB : NA;
            const std::vector<double>& W = k ? *WB : *WA;
            const double Yall = k ? YB : YA, sYall = k ? sYB : sYA;
            auto fallback = [&](const std::string& reason) {
                why[k] = reason;
                return Yall > 0 ? sYall / Yall : 0.0;
            };
            if (!R.has_spedge || !R.has_slop) return fallback("no SPedge/slop");
            if (R.fortran_compat < 0) return fallback("pass_cuts mode unknown (--mode)");
            const bool right = R.electron_arm == 1 || R.electron_arm == 3 || R.electron_arm == 7;
            const bool left = R.electron_arm == 2 || R.electron_arm == 4 || R.electron_arm == 5 || R.electron_arm == 8;
            if (!right && !left) return fallback("electron_arm unknown");
            const char* ep[2] = {right ? "hs" : "ss", right ? "ss" : "hs"};   // e, p prefixes
            const char* var[3] = {"delta", "yptar", "xptar"};
            const std::vector<double>* c[6];
            for (int a = 0; a < 2; ++a)
                for (int v = 0; v < 3; ++v) {
                    c[3 * a + v] = T.col(std::string(ep[a]) + var[v]);
                    if (!c[3 * a + v]) return fallback(std::string("no branch ") + ep[a] + var[v]);
                }
            // the window: printed SPedge (already + slop) -/+ slop used; angles: slop in mr
            double lo[6], hi[6];
            for (int i = 0; i < 6; ++i) {
                const double sl = (i % 3 == 0) ? R.slop_used[i] : R.slop_used[i] / 1000.;
                lo[i] = R.spedge_min[i] + sl;
                hi[i] = R.spedge_max[i] - sl;
            }
            if (R.fortran_compat == 1) hi[3] = R.spedge_max[0] - R.slop_used[3];   // simc.f:237 (D-006)
            double S1 = 0, S2 = 0;
            for (size_t i = 0; i < W.size(); ++i) {
                bool pass = std::isfinite(W[i]);
                for (int k6 = 0; k6 < 6 && pass; ++k6) {
                    const double x = (*c[k6])[i];
                    if (x <= lo[k6] || x >= hi[k6]) pass = false;
                }
                if (pass) { S1 += W[i]; S2 += W[i] * W[i]; }
            }
            if (!(S1 > 0)) return fallback("empty pass_cuts subset");
            double Y = 0, sY = 0;
            if (N.perTried) perTriedMean(S1, S2, N.ntried, Y, sY);
            else { Y = S1; sY = std::sqrt(S2); }
            // consistency: S1*normfac/nevent should be wtcontr
            if (const Scalar* nf = R.scalar("normfac"); nf && R.has_counters && R.nevent() > 0) {
                const double pred = S1 * nf->value / double(R.nevent());
                const double wt = (k ? wb : wa)->value;
                if (std::fabs(pred - wt) > 2e-3 * std::fabs(wt))
                    warnings.push_back(std::string(k ? "B" : "A") + ": pass_cuts subset of the tree gives wtcontr " +
                                       fmt(pred) + " != " + fmt(wt) + " (window mode " +
                                       (R.fortran_compat ? "compat" : "corrected") + ")");
            }
            why[k] = std::string("pass_cuts subset, ") + (R.fortran_compat ? "compat" : "corrected") + " window";
            return Y > 0 ? sY / Y : 0.0;
        };
        const double ra = relerr(0), rb = relerr(1);
        t.ea = std::fabs(wa->value) * ra;
        t.eb = std::fabs(wb->value) * rb;
        t.pull = pull(t.a, t.ea, t.b, t.eb);
        t.p = pullPValue(t.pull);
        t.note = "sigma = wtcontr * relerr: A " + why[0] + "; B " + why[1];
        for (int k = 0; k < 2; ++k)
            if (why[k].rfind("pass_cuts subset", 0) != 0)
                warnings.push_back(std::string(k ? "B" : "A") + ": wtcontr sigma from the whole tree yield (" + why[k] + ")");
        add(t);
        // info: consistency of wtcontr with the tree (pass_cuts subset, simc.f:311-313)
        for (int k = 0; k < 2; ++k) {
            const RunSummary& R = k ? RB : RA;
            const std::vector<double>& W = k ? *WB : *WA;
            const Scalar* nf = R.scalar("normfac");
            if (!nf || !R.has_counters || R.nevent() == 0) continue;
            double S = 0;
            for (double w : W) if (std::isfinite(w)) S += w;
            double pred = S * nf->value / double(R.nevent());
            double wt = (k ? wb : wa)->value;
            if (std::fabs(pred - wt) > 1e-3 * std::fabs(wt))
                warnings.push_back(std::string(k ? "B" : "A") + ": sum(Weight)*normfac/nevent = " + fmt(pred) +
                                   " differs from wtcontr = " + fmt(wt) + " (tree holds all successes, wtcontr only pass_cuts)");
        }
    }
    // AVERAGE.sigcc = sum(sigcc over successes)/nevent (simc.f:248,959). The pion/kaon/delta ntuple calls the same
    // quantity "siglab" (results_write.f:146, NtupleInit.f:129; P9.2): used when the configured branch is absent.
    std::string sigcc_b = cfg.sigccBranch;
    if ((!TA.col(sigcc_b) || !TB.col(sigcc_b)) && TA.col("siglab") && TB.col("siglab")) sigcc_b = "siglab";
    if (const Scalar *sa = RA.scalar("sigcc_average"), *sb = RB.scalar("sigcc_average");
        sa && sb && haveCounters && TA.col(sigcc_b) && TB.col(sigcc_b)) {
        auto err = [&](const RunSummary& R, const std::vector<double>& x) {
            double S1 = 0, S2 = 0, n = 0;
            for (double v : x) if (std::isfinite(v)) { S1 += v; S2 += v * v; n += 1; }
            if (R.ngen < 0) {  // nevent = Ntried
                double Y, sY;
                perTriedMean(S1, S2, double(R.ntried), Y, sY);
                return sY;
            }
            if (n < 2) return 0.0;
            double m = S1 / n;
            return std::sqrt(std::max(0.0, S2 / n - m * m) / n);
        };
        Test t;
        t.category = "norm";
        t.name = "AVERAGE.sigcc";
        t.method = "pull";
        t.a = sa->value; t.b = sb->value;
        t.ea = std::hypot(err(RA, *TA.col(sigcc_b)), sa->halfulp);
        t.eb = std::hypot(err(RB, *TB.col(sigcc_b)), sb->halfulp);
        t.pull = pull(t.a, t.ea, t.b, t.eb);
        t.p = pullPValue(t.pull);
        add(t);
    }
    // normfac: random through Ntried when ngen>0 (info only: same content as acceptance)
    if (const Scalar *na = RA.scalar("normfac"), *nb = RB.scalar("normfac"); na && nb && haveCounters) {
        Test t;
        t.category = (RA.ngen > 0 && RB.ngen > 0) ? "info" : "det";
        t.name = "normfac";
        t.a = na->value; t.b = nb->value;
        if (t.category == "info") {
            t.method = "pull";
            auto rel = [](const RunSummary& R) {
                double a = double(R.ncontribute) / double(R.ntried);
                return std::sqrt(std::max(0.0, 1 - a) / double(R.ncontribute));
            };
            t.ea = t.a * rel(RA);
            t.eb = t.b * rel(RB);
            t.pull = pull(t.a, t.ea, t.b, t.eb);
            t.p = pullPValue(t.pull);
            t.counted = false;
            t.note = "Var(Ntried)/Ntried^2=(1-acc)/Ncontribute; not counted (duplicates acceptance)";
            add(t);
        }
    }
    // deterministic scalars
    {
        std::vector<std::string> det = {"sigcc_central", "charge_mC", "targetfac", "luminosity_ub",
                                        "luminosity_GeV2", "genvol"};
        if (haveCounters && !(RA.ngen > 0 && RB.ngen > 0)) det.push_back("normfac");
        for (const auto& k : det) {
            const Scalar *a = RA.scalar(k), *b = RB.scalar(k);
            if (!a || !b) continue;
            Test t;
            t.category = "det";
            t.name = k;
            t.method = "det";
            t.a = a->value; t.b = b->value;
            double tol = cfg.detRelTol * std::max(std::fabs(t.a), std::fabs(t.b)) + a->halfulp + b->halfulp;
            t.stat = std::fabs(t.a - t.b);
            t.ea = tol;
            t.pass = t.stat <= tol;
            t.p = t.pass ? 1.0 : 0.0;
            t.note = "tol " + fmt(tol, 3);
            add(t);
        }
        if (haveCounters) {
            Test t;
            t.category = "det";
            t.name = "Ngen (request)";
            t.method = "det";
            t.a = double(RA.ngen); t.b = double(RB.ngen);
            t.stat = std::fabs(t.a - t.b);
            t.pass = RA.ngen == RB.ngen;
            t.p = t.pass ? 1.0 : 0.0;
            add(t);
        }
    }
    // stop counters
    if (haveCounters) {
        static const std::set<std::string> nonCategory = {"trials", "hut", "q2_mid_dup"};
        for (const auto& [spec, mA] : RA.stops) {
            auto itB = RB.stops.find(spec);
            if (itB == RB.stops.end()) {
                warnings.push_back("stop counters for " + spec + " only in A");
                continue;
            }
            const auto& mB = itB->second;
            // Only the counters reported by BOTH runs are compared (P9.2): the Fortran .geni does not print every
            // counter (SHMS q3_mid, D1_men, targ_*, cer1/2, BP_in/out; OUTPUT_FORMAT.md), while the SIMC++ summary
            // carries all of them; a counter missing on one side is not a measured 0. Such keys are listed as warnings
            // when they are non-zero; their counts end up in the 'unaccounted' cell of the partition test.
            std::vector<std::string> keys;
            const auto& oB = RB.stopOrder.at(spec);
            for (const auto& k : RA.stopOrder.at(spec))
                if (std::find(oB.begin(), oB.end(), k) != oB.end()) keys.push_back(k);
            for (int side = 0; side < 2; ++side) {
                const auto& order = side ? oB : RA.stopOrder.at(spec);
                const auto& mm = side ? itB->second : mA;
                for (const auto& k : order) {
                    if (std::find(keys.begin(), keys.end(), k) != keys.end()) continue;
                    auto it = mm.find(k);
                    if (it != mm.end() && it->second != 0)
                        warnings.push_back(spec + "." + k + " = " + std::to_string(it->second) + " only in " +
                                           (side ? "B" : "A") + " (not reported by the other run; not compared)");
                }
            }
            auto get = [](const std::map<std::string, int64_t>& m, const std::string& k) {
                auto it = m.find(k);
                return it == m.end() ? int64_t(0) : it->second;
            };
            if (cfg.stopPerCounter) {
                for (const auto& k : keys) {
                    Test t;
                    t.category = "stop";
                    t.name = spec + "." + k + "/Ntried";
                    t.method = "binomial";
                    if (!binomialPull(get(mA, k), RA.ntried, get(mB, k), RB.ntried, t.a, t.b, t.ea, t.eb, t.pull)) continue;
                    t.p = pullPValue(t.pull);
                    add(t);
                }
            }
            // homogeneity of the partition of the spectrometer trials
            std::vector<int64_t> ca, cb;
            int64_t sa = 0, sb = 0;
            for (const auto& k : keys) {
                if (nonCategory.count(k)) continue;
                ca.push_back(get(mA, k));
                cb.push_back(get(mB, k));
                sa += ca.back();
                sb += cb.back();
            }
            int64_t ua = get(mA, "trials") - sa, ub = get(mB, "trials") - sb;
            if (ua < 0 || ub < 0) {
                warnings.push_back(spec + ": stop counters exceed trials; partition test skipped");
                continue;
            }
            ca.push_back(ua);  // 'unaccounted' (e.g. SHMS q3_mid, never printed: simc.f:520)
            cb.push_back(ub);
            Test t;
            t.category = "stop";
            t.name = spec + " partition of trials";
            t.method = "contingency";
            t.a = double(get(mA, "trials"));
            t.b = double(get(mB, "trials"));
            if (contingencyChi2(ca, cb, t.stat, t.ndf, t.p)) add(t);
        }
        for (const auto& [spec, m] : RB.stops)
            if (!RA.stops.count(spec)) warnings.push_back("stop counters for " + spec + " only in B");
    }

    // ----------------------------------------------------------- tree level
    std::vector<std::string> vars;
    std::vector<std::string> onlyA, onlyB;
    for (auto& n : TA.names) {
        if (TB.cols.count(n)) vars.push_back(n);
        else onlyA.push_back(n);
    }
    for (auto& n : TB.names) if (!TA.cols.count(n)) onlyB.push_back(n);
    if (!cfg.variables.empty()) {
        std::vector<std::string> v2;
        for (auto& v : cfg.variables) if (TA.cols.count(v) && TB.cols.count(v)) v2.push_back(v);
        vars = v2;
    }
    if (!onlyA.empty()) { std::string s; for (auto& x : onlyA) s += " " + x; warnings.push_back("branches only in A:" + s); }
    if (!onlyB.empty()) { std::string s; for (auto& x : onlyB) s += " " + x; warnings.push_back("branches only in B:" + s); }

    // Non-finite tree entries (P7.4): counted per run at the EVENT level (an event with any NaN/inf branch among the
    // compared variables or the weight) and reported as an info test; the per-variable tests below exclude them and
    // warn per variable. The Fortran can write NaN events (e.g. K2 norad, D-028 F2).
    {
        auto scan = [&](const TreeData& T, const std::vector<double>* W, std::string& brs) {
            std::vector<char> bad(T.n, 0);
            std::vector<std::string> names(vars);
            if (std::find(names.begin(), names.end(), RA.weightBranch) == names.end()) names.push_back(RA.weightBranch);
            for (const auto& v : names) {
                const std::vector<double>* c = (v == RA.weightBranch) ? W : T.col(v);
                if (!c) continue;
                size_t nv = 0;
                for (size_t i = 0; i < c->size(); ++i)
                    if (!std::isfinite((*c)[i])) { bad[i] = 1; ++nv; }
                if (nv) brs += (brs.empty() ? "" : ",") + v + "(" + std::to_string(nv) + ")";
            }
            size_t n = 0;
            for (char b : bad) n += b ? 1 : 0;
            return n;
        };
        std::string brA, brB;
        const size_t nA = scan(TA, WA, brA), nB = scan(TB, WB, brB);
        Test t;
        t.category = "info";
        t.name = "tree events with non-finite values";
        t.method = "count";
        t.a = double(nA);
        t.b = double(nB);
        t.counted = false;
        t.note = "excluded from the per-variable tests";
        if (nA) t.note += "; A: " + brA;
        if (nB) t.note += "; B: " + brB;
        add(t);
        if (nA || nB)
            warnings.push_back("non-finite tree events (excluded from the per-variable tests): A=" + std::to_string(nA) +
                               (nA ? " [" + brA + "]" : "") + " B=" + std::to_string(nB) + (nB ? " [" + brB + "]" : ""));
    }

    // Empty trees (P7.4): both empty (e.g. Ncontribute = 0 in both runs) -> nothing to compare per variable, the
    // run-level tests decide; exactly one empty -> a real difference (deterministic FAIL).
    const bool emptyA = TA.n == 0, emptyB = TB.n == 0;
    if (emptyA && emptyB) {
        warnings.push_back("both trees are empty (0 entries): the per-variable tests are skipped");
    } else if (emptyA || emptyB) {
        Test t;
        t.category = "det";
        t.name = "tree entries (one tree is empty)";
        t.method = "det";
        t.a = double(TA.n);
        t.b = double(TB.n);
        t.pass = false;
        t.p = 0.0;
        t.note = std::string("tree ") + (emptyA ? "A" : "B") + " has no entries, the other has " +
                 std::to_string(emptyA ? TB.n : TA.n) + "; per-variable tests skipped";
        add(t);
    }

    std::vector<VarResult> vres;
    int hid = 0;
    for (const auto& v : vars) {
        VarResult r;
        r.name = v;
        if (cfg.skip.count(v)) {
            r.status = "skipped";
            vres.push_back(std::move(r));
            continue;
        }
        if (emptyA || emptyB) {
            r.status = "empty";
            vres.push_back(std::move(r));
            continue;
        }
        const auto& xa = *TA.col(v);
        const auto& xb = *TB.col(v);
        // finite events only
        std::vector<double> fa, fwa, fb, fwb;
        size_t nfA = 0, nfB = 0;
        for (size_t i = 0; i < xa.size(); ++i) {
            if (std::isfinite(xa[i]) && std::isfinite((*WA)[i])) { fa.push_back(xa[i]); fwa.push_back((*WA)[i]); }
            else ++nfA;
        }
        for (size_t i = 0; i < xb.size(); ++i) {
            if (std::isfinite(xb[i]) && std::isfinite((*WB)[i])) { fb.push_back(xb[i]); fwb.push_back((*WB)[i]); }
            else ++nfB;
        }
        if (nfA || nfB) warnings.push_back(v + ": non-finite entries A=" + std::to_string(nfA) + " B=" + std::to_string(nfB));

        auto binIt = cfg.binning.find(v);
        r.nbins = (binIt != cfg.binning.end() && binIt->second.nbins > 0) ? binIt->second.nbins : cfg.nbins;
        if (binIt != cfg.binning.end() && binIt->second.fixed) {
            r.lo = binIt->second.lo;
            r.hi = binIt->second.hi;
        } else {
            quantileRange(fa, fb, cfg.qlo, cfg.qhi, r.lo, r.hi);
        }
        if (fa.empty() && fb.empty()) {
            // every entry non-finite in both runs: nothing to compare (reported above and in the warnings)
            r.status = "nonfinite";
            vres.push_back(std::move(r));
            continue;
        }
        if (!(r.hi > r.lo)) {
            // constant branch in both runs: compare the constant value
            r.status = "constant";
            double ca = fa.empty() ? NAN : fa[0], cb = fb.empty() ? NAN : fb[0];
            Test t;
            t.category = "tree";
            t.name = v;
            t.method = "constant";
            t.a = ca; t.b = cb;
            t.pass = (ca == cb);
            t.p = t.pass ? 1.0 : 0.0;
            t.note = "constant branch";
            r.tests.push_back(add(t));
            r.minP = t.p;
            vres.push_back(std::move(r));
            continue;
        }
        const std::string sfx = std::to_string(hid++);
        r.hwA = std::make_unique<TH1D>(("hwA" + sfx).c_str(), v.c_str(), r.nbins, r.lo, r.hi);
        r.hwB = std::make_unique<TH1D>(("hwB" + sfx).c_str(), v.c_str(), r.nbins, r.lo, r.hi);
        r.huA = std::make_unique<TH1D>(("huA" + sfx).c_str(), v.c_str(), r.nbins, r.lo, r.hi);
        r.huB = std::make_unique<TH1D>(("huB" + sfx).c_str(), v.c_str(), r.nbins, r.lo, r.hi);
        for (auto* h : {r.hwA.get(), r.hwB.get(), r.huA.get(), r.huB.get()}) { h->SetDirectory(nullptr); h->Sumw2(); }
        for (size_t i = 0; i < fa.size(); ++i) { r.hwA->Fill(fa[i], fwa[i]); r.huA->Fill(fa[i]); }
        for (size_t i = 0; i < fb.size(); ++i) { r.hwB->Fill(fb[i], fwb[i]); r.huB->Fill(fb[i]); }

        if (cfg.testChi2Abs) {
            // absolute (yield per tried event) chi2 over all bins incl. under/overflow
            Test t;
            t.category = "tree";
            t.name = v;
            t.method = "chi2_abs";
            double chi2 = 0;
            int ndf = 0;
            for (int b = 0; b <= r.nbins + 1; ++b) {
                double a = NA.scale * r.hwA->GetBinContent(b), bb = NB.scale * r.hwB->GetBinContent(b);
                double ea = NA.scale * r.hwA->GetBinError(b), eb = NB.scale * r.hwB->GetBinError(b);
                double s2 = ea * ea + eb * eb;
                if (s2 <= 0) continue;
                chi2 += (a - bb) * (a - bb) / s2;
                ++ndf;
            }
            t.stat = chi2;
            t.ndf = ndf;
            t.p = ndf > 0 ? TMath::Prob(chi2, ndf) : 1.0;
            if (chi2 == 0) t.p = 1.0;
            r.tests.push_back(add(t));
        }
        auto rootChi2 = [&](TH1D* h1, TH1D* h2, const char* opt, const char* meth) {
            Test t;
            t.category = "tree";
            t.name = v;
            t.method = meth;
            double chi2 = 0;
            Int_t ndf = 0, igood = 0;
            std::vector<double> res(r.nbins + 2);
            double p = h1->Chi2TestX(h2, chi2, ndf, igood, opt, res.data());
            t.stat = chi2;
            t.ndf = ndf;
            t.p = (chi2 == 0) ? 1.0 : p;
            if (igood != 0) t.note = "Chi2Test igood=" + std::to_string(igood);
            r.tests.push_back(add(t));
        };
        if (cfg.testChi2WW) rootChi2(r.hwA.get(), r.hwB.get(), "WW UF OF", "chi2_ww");
        if (cfg.testChi2UU) rootChi2(r.huA.get(), r.huB.get(), "UU UF OF", "chi2_uu");
        if (cfg.testKS) {
            Test t;
            t.category = "tree";
            t.name = v;
            t.method = "ks_w";
            weightedKS(fa, fwa, fb, fwb, t.stat, t.p);
            r.tests.push_back(add(t));
        }
        for (auto* t : r.tests) r.minP = std::min(r.minP, t->p);
        vres.push_back(std::move(r));
    }

    // ----------------------------------------------------------- verdict
    int N = 0, F = 0, gross = 0;
    double minP = 1;
    std::string minPName;
    for (auto& tp : tests) {
        Test& t = *tp;
        if (t.method != "det" && t.method != "constant") t.pass = !(t.p < cfg.alpha);
        if (!t.counted) continue;
        ++N;
        if (!t.pass) ++F;
        if (t.method == "det" || t.method == "constant") { if (!t.pass) ++gross; continue; }
        if (t.p < cfg.grossP) ++gross;
        if (t.p < minP) { minP = t.p; minPName = t.name + " [" + t.method + "]"; }
    }
    double expected = N * cfg.alpha;
    double pBinom = (F == 0) ? 1.0 : TMath::BinomialI(cfg.alpha, N, F);
    bool pass = (pBinom >= cfg.globalAlpha) && gross == 0;
    std::string verdict = pass ? "PASS" : "FAIL";

    double maxPull = 0;
    std::string maxPullName;
    for (auto& tp : tests)
        if (tp->counted && !std::isnan(tp->pull) && std::fabs(tp->pull) > std::fabs(maxPull)) {
            maxPull = tp->pull;
            maxPullName = tp->name;
        }

    // ----------------------------------------------------------- outputs
    std::vector<std::string> L;
    auto line = [&](const std::string& s) { L.push_back(s); };
    {
        std::ostringstream o;
        o << "simc_compare  A=" << labelA << "  B=" << labelB;
        line(o.str());
        line("A: " + RA.source + "  [" + RA.format + "]  tree " + RA.treeFile + " (" + std::to_string(TA.n) + " entries)");
        line("B: " + RB.source + "  [" + RB.format + "]  tree " + RB.treeFile + " (" + std::to_string(TB.n) + " entries)");
        line("normalisation: A " + NA.desc + ", B " + NB.desc);
        std::ostringstream v;
        v << "VERDICT: " << verdict << "   tests counted N=" << N << "  failed (p<" << cfg.alpha << ") F=" << F
          << "  expected " << std::fixed << std::setprecision(2) << expected << "  fraction "
          << std::setprecision(4) << (N ? double(F) / N : 0.0) << "  P(X>=F|N,alpha)=" << std::scientific
          << std::setprecision(3) << pBinom << " (need >= " << cfg.globalAlpha << ")  gross failures="
          << gross << " (det mismatch or p<" << cfg.grossP << ")";
        line(v.str());
        std::ostringstream w;
        w << "worst p = " << std::scientific << std::setprecision(3) << minP << "  (" << minPName << ")   largest |pull| = "
          << std::fixed << std::setprecision(3) << maxPull << "  (" << maxPullName << ")";
        line(w.str());
        line("");
        line("RUN-LEVEL TESTS");
        std::ostringstream h;
        h << std::left << std::setw(46) << "quantity" << std::setw(12) << "method" << std::right << std::setw(14) << "A"
          << std::setw(11) << "errA" << std::setw(14) << "B" << std::setw(11) << "errB" << std::setw(9) << "pull"
          << std::setw(11) << "p" << "  status";
        line(h.str());
        for (auto& tp : tests) {
            const Test& t = *tp;
            if (t.category == "tree") continue;
            std::ostringstream o;
            o << std::left << std::setw(46) << t.name.substr(0, 45) << std::setw(12) << t.method << std::right
              << std::setw(14) << fmt(t.a) << std::setw(11) << fmt(t.ea, 3) << std::setw(14) << fmt(t.b)
              << std::setw(11) << fmt(t.eb, 3) << std::setw(9) << fmt(t.pull, 3) << std::setw(11) << fmt(t.p, 3)
              << "  " << (t.pass ? "ok" : "FAIL") << (t.counted ? "" : " (info)")
              << (t.method == "contingency" ? "  chi2/ndf=" + fmt(t.stat, 4) + "/" + std::to_string(t.ndf) : "")
              << (t.note.empty() ? "" : "  " + t.note);
            line(o.str());
        }
        line("");
        line("TREE VARIABLES (p-values; chi2_abs: yield/tried incl. normalisation, chi2_ww/uu: TH1::Chi2Test shape, ks_w: unbinned weighted KS)");
        std::ostringstream h2;
        h2 << std::left << std::setw(11) << "variable" << std::setw(10) << "status" << std::right << std::setw(12) << "lo"
           << std::setw(12) << "hi" << std::setw(12) << "chi2_abs" << std::setw(12) << "chi2_ww" << std::setw(12)
           << "chi2_uu" << std::setw(12) << "ks_w" << "  result";
        line(h2.str());
        for (auto& r : vres) {
            std::map<std::string, double> pm;
            bool ok = true;
            for (auto* t : r.tests) { pm[t->method] = t->p; ok = ok && t->pass; }
            auto P = [&](const char* m) { return pm.count(m) ? fmt(pm[m], 3) : std::string("-"); };
            std::ostringstream o;
            o << std::left << std::setw(11) << r.name << std::setw(10) << r.status << std::right << std::setw(12)
              << fmt(r.lo, 4) << std::setw(12) << fmt(r.hi, 4) << std::setw(12) << P("chi2_abs") << std::setw(12)
              << P("chi2_ww") << std::setw(12) << P("chi2_uu") << std::setw(12) << P("ks_w") << "  "
              << ((r.status == "skipped" || r.status == "empty" || r.status == "nonfinite") ? "-" : (ok ? "ok" : "FAIL"));
            line(o.str());
        }
        if (!warnings.empty()) {
            line("");
            line("WARNINGS");
            for (auto& w : warnings) line("  " + w);
        }
    }
    {
        std::ofstream o(fs::path(outDir) / "summary.txt");
        for (auto& s : L) o << s << "\n";
    }
    {
        json j;
        j["verdict"] = verdict;
        j["labels"] = {{"A", labelA}, {"B", labelB}};
        j["sources"] = {{"A", RA.source}, {"B", RB.source}};
        j["trees"] = {{"A", RA.treeFile}, {"B", RB.treeFile}};
        j["tree_entries"] = {{"A", TA.n}, {"B", TB.n}};
        j["config"] = {{"alpha", cfg.alpha}, {"global_alpha", cfg.globalAlpha}, {"gross_p", cfg.grossP},
                       {"nbins", cfg.nbins}, {"quantile_lo", cfg.qlo}, {"quantile_hi", cfg.qhi},
                       {"det_rel_tol", cfg.detRelTol}, {"norm_a", cfg.normA}, {"norm_b", cfg.normB}};
        j["counts"] = {{"n_tests", N}, {"n_fail", F}, {"expected_fail", expected}, {"fail_fraction", N ? double(F) / N : 0.0},
                       {"p_binomial", pBinom}, {"n_gross", gross}};
        j["worst"] = {{"p", minP}, {"p_test", minPName}, {"pull", maxPull}, {"pull_test", maxPullName}};
        json arr = json::array();
        for (auto& tp : tests) {
            const Test& t = *tp;
            auto num = [](double x) { return std::isfinite(x) ? json(x) : json(nullptr); };
            arr.push_back({{"category", t.category}, {"name", t.name}, {"method", t.method}, {"a", num(t.a)},
                           {"err_a", num(t.ea)}, {"b", num(t.b)}, {"err_b", num(t.eb)}, {"stat", num(t.stat)},
                           {"ndf", t.ndf}, {"pull", std::isnan(t.pull) ? json(nullptr) : json(std::isfinite(t.pull) ? t.pull : (t.pull > 0 ? 1e300 : -1e300))},
                           {"p", t.p}, {"pass", t.pass}, {"counted", t.counted}, {"note", t.note}});
        }
        j["tests"] = arr;
        json va = json::array();
        for (auto& r : vres) {
            json o = {{"name", r.name}, {"status", r.status}, {"lo", r.lo}, {"hi", r.hi}, {"nbins", r.nbins}, {"min_p", r.minP}};
            for (auto* t : r.tests) o[t->method] = t->p;
            va.push_back(o);
        }
        j["variables"] = va;
        j["warnings"] = warnings;
        std::ofstream o(fs::path(outDir) / "summary.json");
        o << j.dump(2) << "\n";
    }

    // ----------------------------------------------------------- PDF
    if (cfg.pdf) {
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetPaperSize(28.0, 18.0);  // landscape page matching the 1400x900 canvas
        std::string pdf = (fs::path(outDir) / "compare.pdf").string();
        TCanvas c("c", "c", 1400, 900);
        bool first = true;
        writeLines(c, L, pdf, first);
        for (auto& r : vres) {
            if (r.status != "compared") continue;
            c.Clear();
            std::vector<std::unique_ptr<TObject>> keep;
            auto panel = [&](double x0, double x1, TH1D* ha, TH1D* hb, double sa, double sb, const std::string& ytit, bool shape) {
                c.cd();
                auto* top = new TPad(("t" + ytit).c_str(), "", x0, 0.30, x1, 0.94);
                auto* bot = new TPad(("b" + ytit).c_str(), "", x0, 0.0, x1, 0.30);
                keep.emplace_back(top);
                keep.emplace_back(bot);
                top->SetBottomMargin(0.02);
                bot->SetTopMargin(0.02);
                bot->SetBottomMargin(0.3);
                top->Draw();
                bot->Draw();
                auto* a = static_cast<TH1D*>(ha->Clone());
                auto* b = static_cast<TH1D*>(hb->Clone());
                a->SetDirectory(nullptr);
                b->SetDirectory(nullptr);
                keep.emplace_back(a);
                keep.emplace_back(b);
                if (shape) {
                    double ia = a->Integral(0, a->GetNbinsX() + 1), ib = b->Integral(0, b->GetNbinsX() + 1);
                    if (ia > 0) a->Scale(1.0 / ia);
                    if (ib > 0) b->Scale(1.0 / ib);
                } else {
                    a->Scale(sa);
                    b->Scale(sb);
                }
                a->SetLineColor(kBlack);
                a->SetMarkerColor(kBlack);
                a->SetMarkerStyle(20);
                a->SetMarkerSize(0.5);
                b->SetLineColor(kRed);
                b->SetMarkerColor(kRed);
                b->SetMarkerStyle(24);
                b->SetMarkerSize(0.5);
                top->cd();
                a->SetMaximum(1.15 * std::max(a->GetMaximum(), b->GetMaximum()));
                a->SetMinimum(0);
                a->GetYaxis()->SetTitle(ytit.c_str());
                a->GetXaxis()->SetLabelSize(0);
                a->Draw("E1");
                b->Draw("E1 SAME");
                auto* leg = new TLegend(0.65, 0.78, 0.89, 0.89);
                keep.emplace_back(leg);
                leg->AddEntry(a, labelA.c_str(), "lep");
                leg->AddEntry(b, labelB.c_str(), "lep");
                leg->Draw();
                bot->cd();
                auto* rt = static_cast<TH1D*>(b->Clone());
                rt->SetDirectory(nullptr);
                keep.emplace_back(rt);
                rt->Divide(a);
                rt->SetMinimum(cfg.ratioLo);
                rt->SetMaximum(cfg.ratioHi);
                rt->GetYaxis()->SetTitle("B/A");
                rt->GetYaxis()->SetNdivisions(505);
                rt->GetYaxis()->SetTitleSize(0.10);
                rt->GetYaxis()->SetTitleOffset(0.45);
                rt->GetYaxis()->SetLabelSize(0.08);
                rt->GetXaxis()->SetLabelSize(0.09);
                rt->GetXaxis()->SetTitle(r.name.c_str());
                rt->GetXaxis()->SetTitleSize(0.10);
                rt->Draw("E1");
                auto* ln = new TLine(r.lo, 1, r.hi, 1);
                keep.emplace_back(ln);
                ln->SetLineStyle(2);
                ln->Draw();
            };
            panel(0.0, 0.5, r.hwA.get(), r.hwB.get(), NA.scale, NB.scale, "weighted yield / tried event", false);
            panel(0.5, 1.0, r.huA.get(), r.huB.get(), 1, 1, "unweighted, unit area", true);
            c.cd();
            std::ostringstream hd;
            hd << r.name << "   ";
            for (auto* t : r.tests) {
                hd << t->method << " p=" << std::setprecision(3) << t->p;
                if (t->ndf >= 0) hd << " (" << fmt(t->stat, 4) << "/" << t->ndf << ")";
                hd << "   ";
            }
            hd << (std::all_of(r.tests.begin(), r.tests.end(), [](Test* t) { return t->pass; }) ? "ok" : "FAIL");
            TText tl;
            tl.SetTextSize(0.024);
            tl.DrawTextNDC(0.02, 0.962, hd.str().c_str());
            c.Print(pdf.c_str(), "pdf Landscape");
        }
        c.Clear();
        TText tl;
        tl.SetTextSize(0.03);
        tl.DrawTextNDC(0.3, 0.5, ("end of comparison: " + verdict).c_str());
        c.Print((pdf + ")").c_str(), "pdf Landscape");
    }

    // ----------------------------------------------------------- console + exit
    for (size_t i = 0; i < std::min<size_t>(L.size(), 6); ++i) std::cout << L[i] << "\n";
    std::cout << "outputs: " << outDir << "/{summary.txt,summary.json" << (cfg.pdf ? ",compare.pdf" : "") << "}\n";

    bool ok = true;
    if (expect.empty()) ok = pass;
    else if (expect == "PASS") ok = pass;
    else if (expect == "FAIL") ok = !pass;
    else if (expect == "IDENTICAL") {
        for (auto& tp : tests) {
            if (!tp->counted) continue;
            if (tp->p < 1.0 - 1e-12 || (!std::isnan(tp->pull) && tp->pull != 0.0) || !tp->pass) {
                std::cout << "not identical: " << tp->name << " [" << tp->method << "] p=" << tp->p << " pull=" << tp->pull << "\n";
                ok = false;
            }
        }
    } else {
        std::cerr << "simc_compare: unknown --expect " << expect << "\n";
        return 2;
    }
    for (auto& v : expectFlagged) {
        bool flagged = false;
        for (auto& r : vres)
            if (r.name == v)
                for (auto* t : r.tests) flagged = flagged || !t->pass;
        std::cout << "expect-flagged " << v << ": " << (flagged ? "flagged" : "NOT flagged") << "\n";
        ok = ok && flagged;
    }
    if (!expect.empty() || !expectFlagged.empty())
        std::cout << "expectation " << (ok ? "MET" : "NOT MET") << "\n";
    return ok ? 0 : 1;
}
