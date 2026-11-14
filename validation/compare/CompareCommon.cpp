// validation/compare/CompareCommon.cpp — see CompareCommon.h (task P15.3 split of simc_compare.cpp).
#include "CompareCommon.h"

namespace simccmp {

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

std::string fmt(double v, int prec) {
    if (std::isnan(v)) return "-";
    std::ostringstream o;
    o << std::setprecision(prec) << v;
    return o.str();
}

// Per-run normalisation of the event weight for "yield per tried event" plots/tests.
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


}  // namespace simccmp
