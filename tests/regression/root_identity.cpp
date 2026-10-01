// tests/regression/root_identity.cpp — bitwise identity of two ROOT output files (task P7.4).
//
//   simc_root_identity A.root B.root
//
// Every key of A must exist in B with the same class, and vice versa. Compared:
//   * TTree: entry count, the branch list (names, leaf types, lengths) and, for every entry, the raw bytes of every
//     leaf buffer (memcmp; NaNs compare equal when their bit patterns are equal);
//   * TH1 (any dimension): name, title, number of bins, axis limits, every bin content and error (including under- and
//     overflow), entries and the statistics array — all bitwise;
//   * other objects: only the class name.
// The TFile/TKey date stamps are NOT compared (they change with every run).
// Exit status: 0 identical, 1 different, 2 usage / cannot open.
#include <TClass.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TTree.h>

#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

std::map<std::string, std::string> keys(TFile& f) {
    std::map<std::string, std::string> m;
    for (TObject* o : *f.GetListOfKeys()) {
        auto* k = static_cast<TKey*>(o);
        m[k->GetName()] = k->GetClassName();
    }
    return m;
}

long compareTree(TTree& a, TTree& b, const std::string& name) {
    long bad = 0;
    if (a.GetEntries() != b.GetEntries()) {
        std::printf("  %s: entries %lld vs %lld\n", name.c_str(), a.GetEntries(), b.GetEntries());
        return 1;
    }
    TObjArray* la = a.GetListOfLeaves();
    TObjArray* lb = b.GetListOfLeaves();
    if (la->GetEntriesFast() != lb->GetEntriesFast()) {
        std::printf("  %s: %d leaves vs %d\n", name.c_str(), la->GetEntriesFast(), lb->GetEntriesFast());
        return 1;
    }
    std::vector<std::pair<TLeaf*, TLeaf*>> leaves;
    for (int i = 0; i < la->GetEntriesFast(); ++i) {
        auto* x = static_cast<TLeaf*>(la->At(i));
        auto* y = static_cast<TLeaf*>(lb->At(i));
        if (std::strcmp(x->GetName(), y->GetName()) != 0 || std::strcmp(x->GetTypeName(), y->GetTypeName()) != 0) {
            std::printf("  %s: leaf %d is %s/%s vs %s/%s\n", name.c_str(), i, x->GetName(), x->GetTypeName(),
                        y->GetName(), y->GetTypeName());
            return 1;
        }
        leaves.emplace_back(x, y);
    }
    std::map<std::string, long> diffs;
    for (Long64_t e = 0; e < a.GetEntries(); ++e) {
        a.GetEntry(e);
        b.GetEntry(e);
        for (auto& p : leaves) {
            const int na = p.first->GetLen() * p.first->GetLenType();
            const int nb = p.second->GetLen() * p.second->GetLenType();
            if (na != nb || std::memcmp(p.first->GetValuePointer(), p.second->GetValuePointer(), na) != 0) {
                if (diffs[p.first->GetName()]++ == 0)
                    std::printf("  %s: leaf %s differs first at entry %lld (%.9g vs %.9g)\n", name.c_str(),
                                p.first->GetName(), e, p.first->GetValue(), p.second->GetValue());
                ++bad;
            }
        }
    }
    for (const auto& d : diffs) std::printf("  %s: leaf %s differs in %ld entries\n", name.c_str(), d.first.c_str(), d.second);
    std::printf("  %-24s TTree  %lld entries x %zu leaves: %s\n", name.c_str(), a.GetEntries(), leaves.size(),
                bad ? "DIFFERENT" : "identical");
    return bad;
}

long compareHist(TH1& a, TH1& b, const std::string& name) {
    long bad = 0;
    auto fail = [&](const char* what) {
        std::printf("  %s: %s differs\n", name.c_str(), what);
        ++bad;
    };
    if (std::strcmp(a.GetTitle(), b.GetTitle()) != 0) fail("title");
    if (a.GetNcells() != b.GetNcells() || a.GetDimension() != b.GetDimension()) {
        fail("binning");
        return bad;
    }
    const TAxis* ax[3] = {a.GetXaxis(), a.GetYaxis(), a.GetZaxis()};
    const TAxis* bx[3] = {b.GetXaxis(), b.GetYaxis(), b.GetZaxis()};
    for (int i = 0; i < 3; ++i)
        if (!sameBits(ax[i]->GetXmin(), bx[i]->GetXmin()) || !sameBits(ax[i]->GetXmax(), bx[i]->GetXmax()))
            fail("axis limits");
    for (int c = 0; c < a.GetNcells(); ++c) {
        if (!sameBits(a.GetBinContent(c), b.GetBinContent(c))) { fail("bin content"); break; }
        if (!sameBits(a.GetBinError(c), b.GetBinError(c))) { fail("bin error"); break; }
    }
    if (!sameBits(a.GetEntries(), b.GetEntries())) fail("entries");
    double sa[TH1::kNstat] = {0}, sb[TH1::kNstat] = {0};
    a.GetStats(sa);
    b.GetStats(sb);
    for (int i = 0; i < TH1::kNstat; ++i)
        if (!sameBits(sa[i], sb[i])) { fail("statistics"); break; }
    if (bad) std::printf("  %-24s %s: DIFFERENT\n", name.c_str(), a.ClassName());
    return bad;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: simc_root_identity A.root B.root\n");
        return 2;
    }
    std::unique_ptr<TFile> fa(TFile::Open(argv[1], "READ")), fb(TFile::Open(argv[2], "READ"));
    if (!fa || fa->IsZombie() || !fb || fb->IsZombie()) {
        std::fprintf(stderr, "simc_root_identity: cannot open %s or %s\n", argv[1], argv[2]);
        return 2;
    }
    const auto ka = keys(*fa), kb = keys(*fb);
    long bad = 0;
    if (ka != kb) {
        std::printf("  the key lists differ (%zu vs %zu keys)\n", ka.size(), kb.size());
        ++bad;
    }
    long nh = 0;
    for (const auto& k : ka) {
        auto it = kb.find(k.first);
        if (it == kb.end() || it->second != k.second) continue;
        TObject* oa = fa->Get(k.first.c_str());
        TObject* ob = fb->Get(k.first.c_str());
        if (!oa || !ob) {
            std::printf("  %s: cannot read\n", k.first.c_str());
            ++bad;
            continue;
        }
        if (auto* ta = dynamic_cast<TTree*>(oa)) {
            bad += compareTree(*ta, *static_cast<TTree*>(ob), k.first);
        } else if (auto* ha = dynamic_cast<TH1*>(oa)) {
            bad += compareHist(*ha, *static_cast<TH1*>(ob), k.first);
            ++nh;
        }
    }
    if (nh) std::printf("  %ld histograms compared\n", nh);
    std::printf("%s vs %s: %s\n", argv[1], argv[2], bad ? "DIFFERENT" : "IDENTICAL");
    return bad ? 1 : 0;
}
