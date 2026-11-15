// validation/compare/CompareCommon.h — simc_compare internals:
// configuration, tree loading, test/result records, normalisation, text pages, the CLI helpers and the comparison
// session shared by the run-level (CompareRunLevel.cpp), tree-level (CompareTreeLevel.cpp) and report
// (CompareReport.cpp) stages. simc_compare.cpp keeps main (CLI, loading, console, exit status).
#ifndef SIMC_VALIDATION_COMPARECOMMON_H
#define SIMC_VALIDATION_COMPARECOMMON_H

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

namespace simccmp {

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

void loadConfig(const std::string& path, Config& c);
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

TreeData loadTree(const std::string& file, const std::string& name);
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

std::string fmt(double v, int prec = 6);
struct Norm {
    double scale = 1;     // multiplies Weight
    bool perTried = true; // variance formula (S2 - S1^2/N)/N^2 applies
    double ntried = 0;
    std::string desc;
};

Norm makeNorm(const std::string& mode, const RunSummary& s, size_t entries, std::vector<std::string>& warn, const std::string& which);
void quantileRange(const std::vector<double>& a, const std::vector<double>& b, double qlo, double qhi, double& lo, double& hi);
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
void writeLines(TCanvas& c, const std::vector<std::string>& lines, const std::string& pdf, bool& first);
void usage();
std::vector<std::string> splitComma(const std::string& s);

/// The state of one comparison (the locals of the former single main()).
struct Session {
    Config& cfg;
    std::vector<std::string>& warnings;
    RunSummary &RA, &RB;
    TreeData &TA, &TB;
    std::vector<std::unique_ptr<Test>>& tests;
    Norm &NA, &NB;
    const std::vector<double>* WA;
    const std::vector<double>* WB;
    bool haveCounters;
    std::string labelA, labelB, outDir;
    std::vector<VarResult> vres;
    Test* add(Test t) {
        tests.push_back(std::make_unique<Test>(std::move(t)));
        return tests.back().get();
    }
};

/// Run-level tests: normalisation pulls, wtcontr, sigcc, deterministic scalars, stop counters (CompareRunLevel.cpp).
void RunLevelTests(Session& st);
/// Tree-level tests per variable (CompareTreeLevel.cpp); fills st.vres.
void TreeLevelTests(Session& st);
/// Verdict, summary.txt, summary.json and compare.pdf (CompareReport.cpp). Returns the PASS flag; L = the text lines.
bool VerdictAndReports(Session& st, std::vector<std::string>& L);

}  // namespace simccmp

#endif  // SIMC_VALIDATION_COMPARECOMMON_H
