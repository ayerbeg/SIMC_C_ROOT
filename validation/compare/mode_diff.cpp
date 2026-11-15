// validation/compare/mode_diff.cpp — the simc_mode_diff tool.
//
// Event-by-event difference of two SIMC++ runs made with the SAME seed and input but different modes (fortran_compat
// vs the corrected default, D-031). When a correction does not change the random-number consumption, the two runs
// generate the same events and the difference is exact; this tool measures it:
//
//   simc_mode_diff <a.root> <b.root> [Pm bin width GeV/c (default 0.05)]
//
// Prints: the entries of both h10 trees; per branch the number of entries that differ (bitwise, float); the sums of
// Weight and their ratio B/A; for the events whose Weight differs, the ratio B/A in bins of Pm (the recon missing
// momentum branch). Exit 0 always (it measures, it does not judge), 2 on I/O errors.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s a.root b.root [Pm bin width GeV/c]\n", argv[0]);
        return 2;
    }
    const double pmBin = argc > 3 ? std::atof(argv[3]) : 0.05;
    std::unique_ptr<TFile> fa(TFile::Open(argv[1])), fb(TFile::Open(argv[2]));
    TTree* ta = fa ? dynamic_cast<TTree*>(fa->Get("h10")) : nullptr;
    TTree* tb = fb ? dynamic_cast<TTree*>(fb->Get("h10")) : nullptr;
    if (!ta || !tb) {
        std::fprintf(stderr, "cannot read h10 from both files\n");
        return 2;
    }
    std::vector<std::string> names;
    for (auto* o : *ta->GetListOfBranches()) names.emplace_back(o->GetName());
    std::vector<float> va(names.size()), vb(names.size());
    int iW = -1, iPm = -1;
    for (std::size_t i = 0; i < names.size(); ++i) {
        ta->SetBranchAddress(names[i].c_str(), &va[i]);
        tb->SetBranchAddress(names[i].c_str(), &vb[i]);
        if (names[i] == "Weight") iW = static_cast<int>(i);
        if (names[i] == "Pm") iPm = static_cast<int>(i);
    }
    const Long64_t na = ta->GetEntries(), nb = tb->GetEntries(), n = std::min(na, nb);
    std::printf("entries: A %lld, B %lld\n", na, nb);
    std::vector<long> ndiff(names.size(), 0);
    double sa = 0, sb = 0;
    long nwdiff = 0, nfiniteA = 0;
    struct Bin { double wa{0}, wb{0}; long n{0}; };
    std::map<int, Bin> bins;
    for (Long64_t k = 0; k < n; ++k) {
        ta->GetEntry(k);
        tb->GetEntry(k);
        for (std::size_t i = 0; i < names.size(); ++i)
            if (std::memcmp(&va[i], &vb[i], 4) != 0 && !(std::isnan(va[i]) && std::isnan(vb[i]))) ++ndiff[i];
        if (iW >= 0) {
            const double wa = va[iW], wb = vb[iW];
            if (std::isfinite(wa)) {
                sa += wa;
                ++nfiniteA;
            }
            if (std::isfinite(wb)) sb += wb;
            if (std::memcmp(&va[iW], &vb[iW], 4) != 0) ++nwdiff;
            if (iPm >= 0) {
                Bin& b = bins[static_cast<int>(std::floor(va[iPm] / pmBin))];
                if (std::isfinite(wa)) b.wa += wa;
                if (std::isfinite(wb)) b.wb += wb;
                ++b.n;
            }
        }
    }
    std::printf("branches with differences (entries differing, bitwise float):\n");
    int nb_ = 0;
    for (std::size_t i = 0; i < names.size(); ++i)
        if (ndiff[i]) {
            std::printf("  %-10s %ld\n", names[i].c_str(), ndiff[i]);
            ++nb_;
        }
    if (!nb_) std::printf("  none: the two trees are identical\n");
    std::printf("sum Weight: A %.10g (finite entries %ld), B %.10g, B/A = %.6f; entries with a different Weight: %ld\n", sa,
                nfiniteA, sb, sa != 0 ? sb / sa : 0., nwdiff);
    if (nwdiff && iPm >= 0) {
        std::printf("Weight B/A per recon Pm bin (%.3f GeV/c):\n  %-16s %8s %14s %14s %9s\n", pmBin, "Pm [GeV/c]",
                    "events", "sumW A", "sumW B", "B/A");
        for (const auto& kv : bins)
            std::printf("  [%6.3f,%6.3f) %8ld %14.6g %14.6g %9.4f\n", kv.first * pmBin, (kv.first + 1) * pmBin, kv.second.n,
                        kv.second.wa, kv.second.wb, kv.second.wa != 0 ? kv.second.wb / kv.second.wa : 0.);
    }
    return 0;
}
