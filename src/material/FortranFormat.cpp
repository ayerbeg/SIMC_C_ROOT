// src/material/FortranFormat.cpp — see include/simc/material/FortranFormat.h.
// gfortran list-directed output of one real*8 (write(6,*) x): the 25-character field that the Fortran messages of
// the target, material and driver code print, reproduced so that the logs can be compared with the Fortran ones.
#include "simc/material/FortranFormat.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace simc {
namespace fortran {

// gfortran chooses F20.d (0.1 <= |x| < 1e17, after rounding to 17 significant digits) or 1PE25.16E3 otherwise.
std::string listReal(double x) {
    char buf[64];
    if (std::isnan(x)) {
        std::snprintf(buf, sizeof buf, "%25s", "NaN");
        return buf;
    }
    if (std::isinf(x)) {
        std::snprintf(buf, sizeof buf, "%25s", x > 0 ? "Infinity" : "-Infinity");
        return buf;
    }
    // Decimal exponent after rounding to 17 significant digits.
    std::snprintf(buf, sizeof buf, "%.16e", x);
    int e10 = std::atoi(std::strchr(buf, 'e') + 1);
    if (x == 0.0) e10 = 0;
    if (e10 >= -1 && e10 <= 16) {
        // F20.d with d = 17 - (digits before the point), then 5 blanks (the unused exponent field)
        int d = (e10 >= 0) ? 16 - e10 : 17;
        char f[48];
        std::snprintf(f, sizeof f, "%20.*f", d, x);
        std::snprintf(buf, sizeof buf, "%s     ", f);
        return buf;
    }
    // 1PE25.16E3: mantissa d.dddddddddddddddd, exponent with sign and 3 digits
    char m[48];
    std::snprintf(m, sizeof m, "%.16E", x);
    char* e = std::strchr(m, 'E');
    int ex = std::atoi(e + 1);
    *e = '\0';
    char out[64];
    std::snprintf(out, sizeof out, "%sE%c%03d", m, ex < 0 ? '-' : '+', std::abs(ex));
    std::snprintf(buf, sizeof buf, "%25s", out);
    return buf;
}

}  // namespace fortran
}  // namespace simc
