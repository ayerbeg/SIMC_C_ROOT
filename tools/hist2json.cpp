// tools/hist2json.cpp
//
// Convert a Fortran SIMC .hist (+ .geni) into the run-summary JSON
// (docs/summary.schema.json), so the comparison tool reads both codes in one format.
//
// Usage: hist2json run.hist [run.geni] [-o out.json] [--text out.txt] [--run NAME] [-q]
//   The .geni defaults to the .hist path with the extension replaced; it is optional.
//   Parser warnings (unrecognised lines) go to stderr unless -q. Exit 0 ok, 1 usage, 2 I/O.
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "simc/output/SummaryJson.h"

using namespace simc;

namespace {
bool slurp(const std::string& path, std::string& out) {
    std::ifstream f(path);
    if (!f) return false;
    std::stringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return true;
}
}  // namespace

int main(int argc, char** argv) {
    std::string hist, geni, out, text, run;
    bool quiet = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-o" && i + 1 < argc) out = argv[++i];
        else if (a == "--text" && i + 1 < argc) text = argv[++i];
        else if (a == "--run" && i + 1 < argc) run = argv[++i];
        else if (a == "-q") quiet = true;
        else if (a == "-h" || a == "--help") {
            std::cerr << "usage: hist2json run.hist [run.geni] [-o out.json] [--text out.txt] [--run NAME] [-q]\n";
            return 0;
        } else if (hist.empty()) hist = a;
        else if (geni.empty()) geni = a;
        else { std::cerr << "hist2json: too many arguments\n"; return 1; }
    }
    if (hist.empty()) { std::cerr << "usage: hist2json run.hist [run.geni] [-o out.json]\n"; return 1; }
    const auto slash = hist.find_last_of('/');
    std::string stem = hist.substr(slash == std::string::npos ? 0 : slash + 1);
    if (stem.size() > 5 && stem.compare(stem.size() - 5, 5, ".hist") == 0) stem.resize(stem.size() - 5);
    if (run.empty()) run = stem;
    if (geni.empty()) {
        geni = hist;
        if (geni.size() > 5 && geni.compare(geni.size() - 5, 5, ".hist") == 0) geni.replace(geni.size() - 5, 5, ".geni");
        else geni += ".geni";
    }
    std::string htxt, gtxt;
    if (!slurp(hist, htxt)) { std::cerr << "hist2json: cannot read " << hist << "\n"; return 2; }
    if (!slurp(geni, gtxt) && !quiet) std::cerr << "hist2json: no .geni (" << geni << "), stop counters empty\n";

    std::vector<std::string> warnings;
    Json j = ParseFortranHist(htxt, gtxt, run, &warnings);
    if (!quiet) for (const auto& w : warnings) std::cerr << "hist2json: " << w << "\n";
    // make_root_tree output next to the .hist, if present
    std::string root = hist.substr(0, hist.size() - (hist.size() > 5 && hist.compare(hist.size() - 5, 5, ".hist") == 0 ? 5 : 0)) + ".root";
    if (std::ifstream(root).good()) j["tree"] = {{"file", root}, {"name", "h10"}, {"weight_branch", "Weight"}};

    if (out.empty()) std::cout << j.dump(2) << "\n";
    else {
        std::ofstream f(out);
        if (!f) { std::cerr << "hist2json: cannot write " << out << "\n"; return 2; }
        f << j.dump(2) << "\n";
    }
    if (!text.empty()) {
        std::ofstream f(text);
        f << FormatSummaryText(j);
    }
    return 0;
}
