// src/driver/RunReportInternal.h — private to src/driver (the report writer is split over RunReport*.cpp): the gfortran
// formatted-record helpers shared by the .hist (RunReportHist.cpp), .geni (RunReportGeni.cpp) and summary writers.
#ifndef SIMC_SRC_DRIVER_RUNREPORTINTERNAL_H
#define SIMC_SRC_DRIVER_RUNREPORTINTERNAL_H

#include <cmath>
#include <initializer_list>
#include <sstream>
#include <string>
#include <utility>

#include "simc/driver/FortranWriter.h"
#include "simc/driver/RunReport.h"

namespace simc {
namespace driver {
namespace report_detail {

using fortran::List;
using fortran::Record;

// gfortran MAX(a, b) for reals: b if b > a or a is NaN
inline double fmaxF(double a, double b) { return (b > a || std::isnan(a)) ? b : a; }
inline double fminF(double a, double b) { return (b < a || std::isnan(a)) ? b : a; }

class Out {
public:
    void line(const Record& r) { os_ << r.str() << "\n"; }
    void line(const std::string& s) { os_ << s << "\n"; }
    void list(const std::string& s) { line(List(s)); }   // write(iun,*) 'text'
    std::string str() const { return os_.str(); }

private:
    std::ostringstream os_;
};

// RunReportRecords.cpp: simc.f:report l.667-788 reaction lines and the record formats of the report
void reactionLines(Out& o, const RunSetup& s);
Record kin1(const std::string& label, double v, const std::string& unit);
Record kin2(const std::string& label, double e, double p, const std::string& unit);
Record cen1(const std::string& label, double v, const std::string& unit);
Record targ2(const std::string& l1, double v1, const std::string& u1, const std::string& l2, double v2, const std::string& u2);
Record targ3(const std::string& l, const TargetFunction& tf, const std::string& u);
Record targ1(const std::string& l, double v, const std::string& u);
Record flags(std::initializer_list<std::pair<const char*, bool>> items);
Record lim(const std::string& l, double umin, double umax, double lo, double hi, const std::string& u);
Record limHeader1();
Record limHeader2();
Record slopLine(const std::string& a, const std::string& b, double used, double lo, double hi, const std::string& u);
Record rad4(const std::string& l, std::initializer_list<double> v);
Record misc(const std::string& l, double v, const std::string& u);
Record recon(const std::string& l, double a, double b, const std::string& u);
Record spedge(const std::string& l, const Cut& c, const std::string& u);

}  // namespace report_detail
}  // namespace driver
}  // namespace simc

#endif  // SIMC_SRC_DRIVER_RUNREPORTINTERNAL_H
