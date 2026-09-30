// RunSummary.h - run-level numbers of one SIMC run (Fortran or C++), as read
// from a summary JSON (preferred) or from the Fortran .hist/.geni text files
// (fallback parser). See validation/README.md for the formats.
#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace simcval {

// A scalar printed in a text file carries a finite precision: `halfulp` is half
// a unit in the last printed digit (0 for values read from JSON with full
// double precision). Deterministic comparisons use it as an extra tolerance.
struct Scalar {
    double value = 0.0;
    double halfulp = 0.0;
    bool present = false;
};

struct RunSummary {
    std::string source;       // file the summary was read from
    std::string format;       // "fortran-hist" | "summary-json"
    std::string tag;          // run tag (basename)
    std::string producer;     // "fortran" / "simc++" / ...

    // COUNTERS (simc.f:913-916)
    bool has_counters = false;
    int64_t ngen = 0;         // Ngen (request); <0 means |ngen| tries (simc.f:344)
    int64_t ntried = 0;
    int64_t ncontribute = 0;
    int64_t nco_no_rad_prot = 0;
    int64_t random_seed = 0;
    bool has_seed = false;

    // nevent at the end of the event loop (simc.f:169,343-350)
    int64_t nevent() const { return ngen < 0 ? -ngen : ngen; }

    // Normalisation scalars (simc.f:920, 955-965), keys:
    //   wtcontr sigcc_central sigcc_average charge_mC targetfac
    //   luminosity_ub luminosity_GeV2 genvol normfac
    std::map<std::string, Scalar> scalars;

    // Selected flags (FLAGS block of the .hist), e.g. doing_phsp.
    std::map<std::string, bool> flags;

    // Stop counters per spectrometer ("hms","sos","hrsr","hrsl","shms","calo"),
    // keyed by the Fortran COMMON variable name without the prefix
    // (hSTOP_slit_hor -> "slit_hor"). Order of appearance is kept in stopOrder.
    std::map<std::string, std::map<std::string, int64_t>> stops;
    std::map<std::string, std::vector<std::string>> stopOrder;

    // Event tree
    std::string treeFile;
    std::string treeName = "h10";
    std::string weightBranch = "Weight";

    std::vector<std::string> warnings;

    const Scalar* scalar(const std::string& k) const {
        auto it = scalars.find(k);
        return (it == scalars.end() || !it->second.present) ? nullptr : &it->second;
    }
};

// Parse a Fortran .hist file (simc.f:report, l.650-1000). Throws on I/O error.
void parseFortranHist(const std::string& path, RunSummary& s);

// Parse a Fortran .geni file (stop counters, simc.f:462-535). Throws on I/O error.
void parseFortranGeni(const std::string& path, RunSummary& s);

// Parse a summary JSON (validation/README.md "Run summary JSON").
void parseSummaryJson(const std::string& path, RunSummary& s);

// Resolve a run specification:
//   <file>.summary.json | <file>.json  -> summary JSON (tree from its "tree" block,
//                                         else <base>.root next to it)
//   <file>.hist                        -> .hist + <base>.geni + <base>.root
//   <file>.root                        -> <base>.summary.json if present, else
//                                         <base>.hist/.geni if present
//   <directory>                        -> the single *.summary.json, else *.hist in it
// preferJson: if both a summary JSON and a .hist exist, use the JSON.
RunSummary loadRun(const std::string& spec, bool preferJson = true);

// Write a RunSummary as summary JSON (used to produce examples for the C++ writer).
void writeSummaryJson(const RunSummary& s, const std::string& path);

}  // namespace simcval
