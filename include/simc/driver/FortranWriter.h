// include/simc/driver/FortranWriter.h
//
// gfortran formatted-output emulation for the text reports of the driver (.hist, .geni), so that the C++ can
// write files that are text-identical to the Fortran ones (task P7.2, replay A).
//
// A Record is one output line built edit descriptor by edit descriptor, like a Fortran FORMAT:
//     Record r;  r.x(9).a("Ebeam", 12).lit(" = ").f(Ebeam, 15, 4).x(2).a("MeV", 10);
// Rules reproduced (gfortran 15, checked against the reference .hist/.geni files):
//   Aw    w > len: right-justified (leading blanks); w < len: the leftmost w characters; A: the whole string
//   Fw.d  glibc %w.df (same correctly rounded digits as libgfortran); field overflow -> w asterisks
//   Ew.d  0.ddddddE+xx; the leading 0 is dropped when the field is one character short (negative values)
//   Iw, Lw (w-1 blanks + T/F), nX (blanks, but trailing X at the end of a record writes nothing), Tc (tab to
//   column c, 1-based, may move left and overwrite), '/' = end of record (handled by the caller).
//   List-directed (write(u,*)): the record starts with one blank; a character item is copied as is; an
//   integer*4 item is right-justified in a 12-character field (no extra separator after a character item).
#ifndef SIMC_DRIVER_FORTRANWRITER_H
#define SIMC_DRIVER_FORTRANWRITER_H

#include <string>

namespace simc {
namespace fortran {

/// Fw.d
std::string F(double v, int w, int d);
/// Ew.d
std::string E(double v, int w, int d);
/// Iw
std::string I(long v, int w);
/// Lw
std::string L(bool v, int w);
/// Aw
std::string A(const std::string& s, int w);

/// One formatted record (line).
class Record {
public:
    Record& lit(const std::string& s) { return put(s); }                      ///< 'literal' or A without width
    Record& a(const std::string& s, int w) { return put(A(s, w)); }
    Record& f(double v, int w, int d) { return put(F(v, w, d)); }
    Record& e(double v, int w, int d) { return put(E(v, w, d)); }
    Record& i(long v, int w) { return put(I(v, w)); }
    Record& l(bool v, int w) { return put(L(v, w)); }
    Record& x(int n) { pos_ += static_cast<std::size_t>(n); return *this; }   ///< nX
    Record& t(int col) { pos_ = static_cast<std::size_t>(col - 1); return *this; }   ///< Tc
    /// List-directed integer*4 item.
    Record& listInt(long v) { return put(I(v, 12)); }
    std::string str() const { return buf_; }

private:
    Record& put(const std::string& s);
    std::string buf_;
    std::size_t pos_{0};
};

/// A list-directed record start: write(u,*) 'text' -> " text".
inline Record List(const std::string& text = "") {
    Record r;
    r.lit(" " + text);
    return r;
}

}  // namespace fortran
}  // namespace simc

#endif  // SIMC_DRIVER_FORTRANWRITER_H
