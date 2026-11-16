/**
 * @file
 * @ingroup driver
 * @brief Emulation of gfortran formatted output (edit descriptors), so that the text reports are identical to the Fortran ones.
 *
 * gfortran formatted-output emulation for the text reports of the driver (.hist, .geni), so that the C++ can
 * write files that are text-identical to the Fortran ones (task P7.2, replay A).
 *
 * A Record is one output line built edit descriptor by edit descriptor, like a Fortran FORMAT:
 *     Record r;  r.x(9).a("Ebeam", 12).lit(" = ").f(Ebeam, 15, 4).x(2).a("MeV", 10);
 * Rules reproduced (gfortran 15, checked against the reference .hist/.geni files):
 *   Aw    w > len: right-justified (leading blanks); w < len: the leftmost w characters; A: the whole string
 *   Fw.d  glibc %w.df (same correctly rounded digits as libgfortran); field overflow -> w asterisks
 *   Ew.d  0.ddddddE+xx; the leading 0 is dropped when the field is one character short (negative values)
 *   Iw, Lw (w-1 blanks + T/F), nX (blanks, but trailing X at the end of a record writes nothing), Tc (tab to
 *   column c, 1-based, may move left and overwrite), '/' = end of record (handled by the caller).
 *   List-directed (write(u,*)): the record starts with one blank; a character item is copied as is; an
 *   integer*4 item is right-justified in a 12-character field (no extra separator after a character item).
 */
#ifndef SIMC_DRIVER_FORTRANWRITER_H
#define SIMC_DRIVER_FORTRANWRITER_H

#include <string>

/// @addtogroup driver
/// @{

namespace simc {
namespace fortran {

/// Fw.d
/// @param v value to write.
/// @param w field width (characters).
/// @param d digits after the decimal point.
/// @return the field, exactly @p w characters (w asterisks on overflow).
std::string F(double v, int w, int d);
/// Ew.d
/// @param v value to write.
/// @param w field width (characters).
/// @param d digits of the mantissa.
/// @return the field, exactly @p w characters (w asterisks on overflow).
std::string E(double v, int w, int d);
/// Iw
/// @param v value to write.
/// @param w field width (characters).
/// @return the right-justified integer field (w asterisks on overflow).
std::string I(long v, int w);
/// Lw
/// @param v value to write.
/// @param w field width (characters).
/// @return w-1 blanks followed by `T` or `F`.
std::string L(bool v, int w);
/// Aw
/// @param s string to write.
/// @param w field width (characters).
/// @return @p s right-justified in @p w characters, or its leftmost @p w characters if it is longer.
std::string A(const std::string& s, int w);

/// One formatted record (line).
///
/// Built descriptor by descriptor with chained calls; every method returns `*this`. Not thread-safe (plain value
/// object; use one per thread).
class Record {
public:
    /// @param s text copied as is at the current position.
    /// @return this record.
    Record& lit(const std::string& s) { return put(s); }                      ///< 'literal' or A without width
    /// Aw.
    /// @param s string; @param w field width.
    /// @return this record.
    Record& a(const std::string& s, int w) { return put(A(s, w)); }
    /// Fw.d.
    /// @param v value; @param w field width; @param d digits after the decimal point.
    /// @return this record.
    Record& f(double v, int w, int d) { return put(F(v, w, d)); }
    /// Ew.d.
    /// @param v value; @param w field width; @param d mantissa digits.
    /// @return this record.
    Record& e(double v, int w, int d) { return put(E(v, w, d)); }
    /// Iw.
    /// @param v value; @param w field width.
    /// @return this record.
    Record& i(long v, int w) { return put(I(v, w)); }
    /// Lw.
    /// @param v value; @param w field width.
    /// @return this record.
    Record& l(bool v, int w) { return put(L(v, w)); }
    /// @param n number of blanks to skip (a trailing nX writes nothing).
    /// @return this record.
    Record& x(int n) { pos_ += static_cast<std::size_t>(n); return *this; }   ///< nX
    /// @param col 1-based column to move to (may move left and overwrite).
    /// @return this record.
    Record& t(int col) { pos_ = static_cast<std::size_t>(col - 1); return *this; }   ///< Tc
    /// List-directed integer*4 item.
    /// @param v the integer, right-justified in a 12-character field.
    /// @return this record.
    Record& listInt(long v) { return put(I(v, 12)); }
    /// @return the record text built so far.
    std::string str() const { return buf_; }

private:
    Record& put(const std::string& s);
    std::string buf_;
    std::size_t pos_{0};
};

/// A list-directed record start: write(u,*) 'text' -> " text".
/// @param text the first character item (copied after the leading blank).
/// @return a new record.
inline Record List(const std::string& text = "") {
    Record r;
    r.lit(" " + text);
    return r;
}

/// @}
}  // namespace fortran
}  // namespace simc

#endif  // SIMC_DRIVER_FORTRANWRITER_H
