// validation/compare/CompareReport.cpp — verdict, summary.txt/.json and the PDF (task P15.3 split of simc_compare.cpp, no logic change).
#include "CompareCommon.h"

namespace simccmp {

bool VerdictAndReports(Session& st, std::vector<std::string>& L) {
    Config& cfg = st.cfg;
    std::vector<std::string>& warnings = st.warnings;
    RunSummary &RA = st.RA, &RB = st.RB;
    TreeData &TA = st.TA, &TB = st.TB;
    std::vector<std::unique_ptr<Test>>& tests = st.tests;
    Norm &NA = st.NA, &NB = st.NB;
    std::vector<VarResult>& vres = st.vres;
    const std::string& labelA = st.labelA;
    const std::string& labelB = st.labelB;
    const std::string& outDir = st.outDir;
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

    return pass;
}

}  // namespace simccmp
