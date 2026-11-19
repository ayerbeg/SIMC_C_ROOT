// tests/l1/Upstream294.h — task U.1: the values in which the frozen 294f662 Fortran references differ from SIMC++
// because SIMC++ follows the upstream SIMC (59ccdf6) there. Used by the L1 / replay comparators against 294f662 tables,
// captures and runs; every accepted value is counted and the total is printed at exit. Narrow by construction:
//   * z and pt^2 (columns *zhad, *pt2, zi, pt2i): upstream event.f:990-991 / 1338-1342 (commit e4dd175) computes them for
//     every reaction; 294f662 assigned them only for SIDIS (and recon for rho), so the reference holds 0 -> accepted
//     only when the reference is exactly 0 (a non-zero reference is still compared bitwise);
//   * pfermi (column *pfermi): 294f662 results_write.f:208-209 wrote 0/0 = NaN for pfer = 0 (D-051 F3), upstream
//     results_write.f:216-218 writes 0 -> accepted only for reference NaN and C++ exactly 0;
//   * the .hist line `doing_nuc_elast= F` (upstream simc.f:894, 4df7351): dropped from the C++ .hist when the reference
//     .hist does not have it (u1StripHist).
// The upstream behaviour itself is proven bitwise against 59ccdf6 runs (l1_simcrun_ranlux_up_*).
#ifndef SIMC_TESTS_L1_UPSTREAM294_H
#define SIMC_TESTS_L1_UPSTREAM294_H

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

namespace simcl1 {

struct U1Count {
    long n{0};
    ~U1Count() {
        if (n)
            std::printf("  U.1: %ld reference value(s) of the 294f662 Fortran accepted as upstream-SIMC changes "
                        "(z/pt2 never assigned -> computed, e4dd175; pfermi NaN -> 0, D-051 F3)\n", n);
    }
};
inline U1Count& u1Count() {
    static U1Count c;
    return c;
}

inline bool u1EndsWith(const std::string& s, const char* t) {
    const std::string x(t);
    return s.size() >= x.size() && s.compare(s.size() - x.size(), x.size(), x) == 0;
}

/// True (and counted) if (ref, cpp) of column `col` is one of the U.1 upstream changes described above.
inline bool u1UpstreamChange(const std::string& col, double ref, double cpp) {
    const bool zpt = u1EndsWith(col, "zhad") || u1EndsWith(col, "pt2") || col == "zi" || col == "pt2i";
    if (zpt && ref == 0) { ++u1Count().n; return true; }   // any C++ value (also NaN for nu = 0, as upstream)
    if (u1EndsWith(col, "pfermi") && std::isnan(ref) && cpp == 0) { ++u1Count().n; return true; }
    return false;
}

/// The C++ .hist text without the `doing_nuc_elast= ...` line when the reference .hist (294f662) has none.
inline std::string u1StripHist(const std::string& cpp, const std::string& ref) {
    if (ref.find("doing_nuc_elast=") != std::string::npos) return cpp;
    std::istringstream is(cpp);
    std::string out, l;
    while (std::getline(is, l))
        if (l.find("doing_nuc_elast=") == std::string::npos) out += l + "\n";
    return out;
}

}  // namespace simcl1

#endif  // SIMC_TESTS_L1_UPSTREAM294_H
