// validation/compare/CompareTreeLevel.cpp — the tree-level tests of simc_compare: per ntuple variable, the weighted
// histograms of the two runs (chi2 and Kolmogorov-Smirnov), plus the non-finite and empty-tree checks.
#include "CompareCommon.h"

namespace simccmp {

void TreeLevelTests(Session& st) {
    Config& cfg = st.cfg;
    std::vector<std::string>& warnings = st.warnings;
    RunSummary &RA = st.RA;
    TreeData &TA = st.TA, &TB = st.TB;
    Norm &NA = st.NA, &NB = st.NB;
    const std::vector<double>* WA = st.WA;
    const std::vector<double>* WB = st.WB;
    auto add = [&](Test t) { return st.add(std::move(t)); };
    std::vector<VarResult>& vres = st.vres;
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

    // Non-finite tree entries: counted per run at the EVENT level (an event with any NaN/inf branch among the
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

    // Empty trees: both empty (e.g. Ncontribute = 0 in both runs) -> nothing to compare per variable, the
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

}

}  // namespace simccmp
