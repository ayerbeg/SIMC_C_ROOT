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
//
// Task P15.3: the stages are in CompareRunLevel.cpp, CompareTreeLevel.cpp, CompareReport.cpp (helpers: CompareCommon.*).
#include "CompareCommon.h"

using namespace simccmp;


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
        // D-054 (P15.1): the run summaries carry run.fortran_compat (a Fortran .hist is compat); an explicit --mode-a/-b
        // overrides it, with a warning when the two disagree
        auto applyMode = [&](RunSummary& R, int mode, const char* which) {
            if (mode == -2) return;
            if (R.fortran_compat >= 0 && R.fortran_compat != mode)
                warnings.push_back(std::string("--mode-") + which + " " + (mode ? "compat" : "corrected") +
                                   " overrides the run summary (" + (R.fortran_compat ? "compat" : "corrected") + ")");
            R.fortran_compat = mode;
        };
        applyMode(RA, modeA, "a");
        applyMode(RB, modeB, "b");
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
    Session st{cfg, warnings, RA, RB, TA, TB, tests, NA, NB, WA, WB, haveCounters, labelA, labelB, outDir, {}};
    RunLevelTests(st);
    TreeLevelTests(st);
    std::vector<std::string> L;
    const bool pass = VerdictAndReports(st, L);
    const std::vector<VarResult>& vres = st.vres;

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
