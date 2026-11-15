// src/optics/Lfit.cpp — see include/simc/optics/Lfit.h.
// CERNLIB LFIT (unweighted straight-line least-squares fit), used by the detector huts to fit the track through the
// wire-chamber planes.
#include "simc/optics/Lfit.h"

namespace simc {
namespace cern {

// Fortran: cern/lfit.f:LFIT (l.11-60). Accumulators in double (implicit REAL under -fdefault-real-8),
// inputs/outputs float; operation order as in the Fortran (and as in lfit.o).
void lfit(const float* x, const float* y, int l, int key, float& a, float& b, float& e) {
    double count, sumx, sumy, sumxy, sumxx, sumyy, xmed, ymed;
    if (l - 2 < 0) goto l25;   // IF(L-2) 25,1,1
    count = 0.0;
    sumx = 0.0;
    sumy = 0.0;
    sumxy = 0.0;
    sumxx = 0.0;
    sumyy = 0.0;
    for (int j = 0; j < l; ++j) {
        if (y[j] == 0.0f && key == 0) continue;   // FORTRAN-QUIRK: points with Y = 0 are ignored
        sumx = sumx + static_cast<double>(x[j]);
        sumy = sumy + static_cast<double>(y[j]);
        count = count + 1.0;
    }
    if (count <= 1.) goto l25;
    ymed = sumy / count;
    xmed = sumx / count;
    for (int j = 0; j < l; ++j) {
        if (y[j] == 0.0f && key == 0) continue;
        const double scartx = static_cast<double>(x[j]) - xmed;
        const double scarty = static_cast<double>(y[j]) - ymed;
        sumxy = sumxy + scartx * scarty;
        sumxx = sumxx + scartx * scartx;
        sumyy = sumyy + scarty * scarty;
    }
    // FIT PARAMETERS
    if (sumxx == 0.) goto l25;
    a = static_cast<float>(sumxy / sumxx);
    b = static_cast<float>(ymed - static_cast<double>(a) * xmed);
    if (count < 3.) goto l101;
    e = static_cast<float>((sumyy - sumxy * static_cast<double>(a)) / (count - 2.0));
    return;
    // INSUFFICIENT POINTS
l25:
    a = 0.0f;
    b = 0.0f;
l101:
    e = 0.0f;
}

}  // namespace cern
}  // namespace simc
