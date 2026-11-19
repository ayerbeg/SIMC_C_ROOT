// tests/l1/test_l1_nuc_elastic.cpp — L1 of upstream SIMC nuclear_elastic.f (task U.1) against the U.F1 tables
// reference/functions/nuclear_elastic/{sig_nuc_elastic.tsv, fgauss.tsv} (README.md there). Required: bitwise.
//   test_l1_nuc_elastic <dir>     exit 0 pass, 1 fail, 77 tables missing
// Every output column is compared: Fgauss, W1, W2 (NUC_FORM_FACTOR), sigMott(E', theta, Q2), recoil, sig_nuc_elastic.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/physics/nuclear/NuclearElastic.h"
#include "simc/physics/proton/ProtonCrossSections.h"

namespace {
struct Stat {
    std::string name;
    std::size_t n{0}, bitwise{0};
    double maxrel{0};
    void add(double cpp, double ref) {
        ++n;
        if (std::memcmp(&cpp, &ref, sizeof cpp) == 0) { ++bitwise; return; }
        maxrel = std::max(maxrel, simcl1::relDiff(cpp, ref));
    }
    bool print() const {
        const bool ok = n > 0 && bitwise == n;
        std::printf("  %-22s rows %4zu  bitwise %4zu  max|rel| %.3e  %s\n", name.c_str(), n, bitwise, maxrel,
                    ok ? "PASS" : "FAIL");
        return ok;
    }
};
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s <reference/functions/nuclear_elastic>\n", argv[0]); return 2; }
    const std::string dir = argv[1];
    simcl1::Table sig, fg;
    std::string err;
    if (!sig.load(dir + "/sig_nuc_elastic.tsv", err) || !fg.load(dir + "/fgauss.tsv", err)) {
        std::printf("SKIP: %s\n", err.c_str());
        return 77;
    }
    bool ok = true;
    // fgauss.tsv: Fgauss(T = Q2, A)
    {
        Stat s{"Fgauss (fgauss.tsv)"};
        for (std::size_t i = 0; i < fg.rows.size(); ++i) {
            simcl1::Row r(fg, i);
            s.add(simc::nuclear::Fgauss(r["Q2[MeV2]"], r["A"]), r["Fgauss"]);
        }
        ok = s.print() && ok;
    }
    // sig_nuc_elastic.tsv
    {
        Stat sF{"Fgauss"}, sW1{"W1"}, sW2{"W2"}, sM{"sigMott"}, sR{"recoil"}, sS{"sig_nuc_elastic"};
        for (std::size_t i = 0; i < sig.rows.size(); ++i) {
            simcl1::Row r(sig, i);
            const double A = r["A"], Z = r["Z"], M = r["targ_M[MeV]"], Ein = r["Ein[MeV]"], th = r["theta_e[rad]"];
            const double Ep = r["Eprime[MeV]"], Q2 = r["Q2[MeV2]"];
            double W1, W2, FF;
            simc::nuclear::NucFormFactor(Q2, W1, W2, FF, Z, A);
            sF.add(FF, r["Fgauss"]);
            sW1.add(W1, r["W1"]);
            sW2.add(W2, r["W2"]);
            sM.add(simc::proton::sigMott(Ep, th, Q2), r["sigMott[ub/sr]"]);
            sR.add(M / (M + Ein * (1. - std::cos(th))), r["recoil"]);
            sS.add(simc::nuclear::SigNucElastic(Ein, Ep, th, Q2, Z, A, M), r["sig_nuc_elastic[ub/sr]"]);
        }
        for (const Stat* s : {&sF, &sW1, &sW2, &sM, &sR, &sS}) ok = s->print() && ok;
    }
    std::printf("nuclear_elastic L1: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
