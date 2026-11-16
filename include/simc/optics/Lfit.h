/**
 * @file
 * @ingroup optics
 * @brief simc::cern::lfit: the CERNLIB LFIT least-squares straight line, with the precision of the SIMC build.
 *
 * CERNLIB LFIT (Mathlib E250, cern/lfit.f): least-squares straight line y = a*x + b through l points.
 * Fortran: cern/lfit.f:LFIT(X,Y,L,KEY,A,B,E) (l.11-60), as compiled by SIMC.
 *
 * Precision (verified on cern/lfit.o, objdump): X, Y, A, B, E are declared REAL*4, but SIMC compiles with
 * -fdefault-real-8, so the IMPLICITLY typed locals COUNT, SUMX, SUMY, SUMXY, SUMXX, SUMYY, XMED, YMED, SCARTX,
 * SCARTY are REAL*8. The sums are therefore accumulated in double from float inputs (cvtss2sd), and only the
 * three results are rounded to float (cvtsd2ss). B = YMED - A*XMED and E use the ALREADY ROUNDED float A.
 *
 * FORTRAN-QUIRK (CERNLIB): with KEY = 0 every point whose Y is exactly 0.0 (float compare, -0.0 included) is
 * ignored. mc_hms_hut / mc_shms_hut (both call it) rely on this: the planes that measure the other coordinate are set to 0,
 * and a real hit that happens to round to exactly 0.0f is silently dropped as well.
 * Fewer than 2 used points, or all used X equal, gives A = B = E = 0.
 */
#ifndef SIMC_OPTICS_LFIT_H
#define SIMC_OPTICS_LFIT_H

namespace simc {
namespace cern {

/// @addtogroup optics
/// @{

/// Fortran: LFIT(X,Y,L,KEY,A,B,E). a = slope, b = intercept, e = residual variance estimate.
/// @param x abscissae (@p l values).
/// @param y ordinates (@p l values).
/// @param l number of points.
/// @param key 0: ignore points with y == 0.0f; otherwise use all points.
/// @param[out] a slope.
/// @param[out] b intercept.
/// @param[out] e residual variance estimate.
void lfit(const float* x, const float* y, int l, int key, float& a, float& b, float& e);

/// @}

}  // namespace cern
}  // namespace simc

#endif  // SIMC_OPTICS_LFIT_H
