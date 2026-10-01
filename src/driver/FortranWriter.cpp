// src/driver/FortranWriter.cpp — see include/simc/driver/FortranWriter.h.
#include "simc/driver/FortranWriter.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace simc {
namespace fortran {

namespace {
std::string stars(int w) { return std::string(static_cast<std::size_t>(w), '*'); }
std::string rjust(const std::string& s, int w) {
    if (static_cast<int>(s.size()) > w) return stars(w);
    return std::string(static_cast<std::size_t>(w) - s.size(), ' ') + s;
}
}  // namespace

std::string F(double v, int w, int d) {
    if (std::isnan(v)) return rjust("NaN", w);
    if (std::isinf(v)) return rjust(v > 0 ? "Infinity" : "-Infinity", w);
    char buf[512];
    std::snprintf(buf, sizeof buf, "%.*f", d, v);
    std::string s = buf;
    if (static_cast<int>(s.size()) > w) {
        // gfortran drops the optional leading zero of |v| < 1 before giving up
        if (s.rfind("0.", 0) == 0) s = s.substr(1);
        else if (s.rfind("-0.", 0) == 0) s = "-" + s.substr(2);
    }
    return rjust(s, w);
}

std::string E(double v, int w, int d) {
    if (std::isnan(v)) return rjust("NaN", w);
    if (std::isinf(v)) return rjust(v > 0 ? "Infinity" : "-Infinity", w);
    std::string digits;
    int exp10 = 0;
    const bool neg = std::signbit(v);
    if (v == 0.0) {
        digits.assign(static_cast<std::size_t>(d), '0');
        exp10 = 0;
    } else {
        // d significant digits, correctly rounded: d.ddd...e+XX -> 0.dddd...E+(XX+1)
        char buf[512];
        std::snprintf(buf, sizeof buf, "%.*e", d - 1, std::fabs(v));
        const std::string s = buf;
        const std::size_t epos = s.find('e');
        for (std::size_t k = 0; k < epos; ++k)
            if (s[k] != '.') digits += s[k];
        exp10 = std::atoi(s.c_str() + epos + 1) + 1;
    }
    char ebuf[16];
    if (std::abs(exp10) <= 99)
        std::snprintf(ebuf, sizeof ebuf, "E%c%02d", exp10 < 0 ? '-' : '+', std::abs(exp10));
    else
        std::snprintf(ebuf, sizeof ebuf, "%c%03d", exp10 < 0 ? '-' : '+', std::abs(exp10));
    std::string body = "." + digits + ebuf;
    std::string withZero = std::string(neg ? "-" : "") + "0" + body;
    if (static_cast<int>(withZero.size()) <= w) return rjust(withZero, w);
    std::string noZero = std::string(neg ? "-" : "") + body;
    return rjust(noZero, w);
}

std::string I(long v, int w) { return rjust(std::to_string(v), w); }

std::string L(bool v, int w) { return std::string(static_cast<std::size_t>(w - 1), ' ') + (v ? "T" : "F"); }

std::string A(const std::string& s, int w) {
    if (static_cast<int>(s.size()) >= w) return s.substr(0, static_cast<std::size_t>(w));
    return std::string(static_cast<std::size_t>(w) - s.size(), ' ') + s;
}

Record& Record::put(const std::string& s) {
    if (buf_.size() < pos_) buf_.resize(pos_, ' ');
    if (pos_ + s.size() > buf_.size()) buf_.resize(pos_ + s.size(), ' ');
    buf_.replace(pos_, s.size(), s);
    pos_ += s.size();
    return *this;
}

}  // namespace fortran
}  // namespace simc
