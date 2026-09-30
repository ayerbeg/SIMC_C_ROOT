// include/simc/material/FortranFormat.h
//
// Helpers to print warnings with the same text as gfortran list-directed output
// (write(6,*) ...), so L1 tests can compare the WARNING lines of the Fortran driver logs.
// gfortran rules used (REAL*8, checked against the reference driver.log files):
//   * a record starts with one blank;
//   * items are separated by one blank, except between two character items;
//   * a real*8 item is a 25-character field: G-editing with 17 significant digits,
//     F20.d + 5 blanks when 0.1 <= |x| < 1e17 (or x == 0), else 1PE25.16E3.
#ifndef SIMC_MATERIAL_FORTRANFORMAT_H
#define SIMC_MATERIAL_FORTRANFORMAT_H

#include <string>

namespace simc {
namespace fortran {

/// The 25-character list-directed field of a real*8 value.
std::string listReal(double x);

}  // namespace fortran
}  // namespace simc

#endif  // SIMC_MATERIAL_FORTRANFORMAT_H
