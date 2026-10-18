// include/simc/physics/sidis/Cteq5Pdf.h
//
// CTEQ5 parton distribution functions, a bit-for-bit port of cteq5/Ctq5Pdf.f (task P12.1):
//   Cteq5Pdf::load     Ctq5Pdf.f:SetCtq5 + ReadTbl (l.199-276)   table file -> COMMON /CtqPar1/ /CtqPar2/ /XQrange/
//                                                                 /QCDtable/ /Masstbl/
//   Cteq5Pdf::pdf      Ctq5Pdf.f:Ctq5Pdf (l.61-95)              x f(x, Q) of one parton, negative values -> 0
//   Cteq5Pdf::partonX  Ctq5Pdf.f:PartonX (l.97-196)             quadratic (M = 2) interpolation in x and ln(Q/Lambda)
//   polint             Ctq5Pdf.f:POLINT (l.297-334)             Numerical Recipes polynomial interpolation
//
// SIMC calls SetCtq5(iset = 1) (semi_physics.f:121, 300 and semi_dilution.f:79): only cteq5/cteq5m.tbl (CTEQ5M, MSbar)
// is read; a copy is in data/sidis/. The other sets are supported by load() (same file format) but never used.
//
// Numerics. Ctq5Pdf.f has its `Implicit Double Precision` lines commented out, so under -fdefault-real-8 every variable is
// REAL(8) (implicit typing, I-N integer); the literals 0D0/1D0 are REAL(16) but are only compared with or assigned to
// REAL(8) values (exact). Ctq5Pdf.o contains no binary128 arithmetic (only __extenddftf2/__trunctfdf2 of the -ff2c
// return convention, exact), so the port is plain double. The table values are read list-directed (correctly rounded,
// strtod).
//
// State. The SAVEd locals `Warn` (Ctq5Pdf) and `First` (PartonX) and the COMMON blocks are members: one Cteq5Pdf object
// is "one Fortran process" (shared by peepiX and semi_dilution, as the COMMON is).
#ifndef SIMC_PHYSICS_SIDIS_CTEQ5PDF_H
#define SIMC_PHYSICS_SIDIS_CTEQ5PDF_H

#include <array>
#include <ostream>
#include <string>
#include <vector>

namespace simc {
namespace sidis {

/// Fortran: Ctq5Pdf.f:POLINT(XA, YA, N, X, Y, DY) (l.297-334), N <= 10. xa/ya point to N values (XA(1) = xa[0]).
/// Returns Y; DY in `dy`. DEN = 0 (two equal abscissae) is a Fortran PAUSE: throws simc::FortranStop.
double polint(const double* xa, const double* ya, int n, double x, double& dy);

class Cteq5Pdf {
public:
    static constexpr int MXX = 105, MXQ = 25, MXF = 6;
    static constexpr int MXPQX = (MXF * 2 + 2) * MXQ * MXX;
    static constexpr int M = 2, M1 = M + 1;
    /// Ctq5Pdf.f:207-212 Flnm(1..9) without the directory 'cteq5/'
    static const std::array<const char*, 9> kFiles;

    /// SetCtq5(iset) + ReadTbl: reads <dir>/<kFiles[iset-1]>. Returns false and sets `err` if iset is out of range
    /// (Fortran: 'Invalid Iset number in SetCtq5' STOP), if the file is missing (Fortran: STOP l.226) or malformed.
    /// The table read stops silently at an end of file (IOSTAT=IRET, l.272): `points_read()` tells how many were read.
    bool load(const std::string& dir, int iset, std::string* err = nullptr);
    bool loaded() const { return iset_ > 0; }
    int iset() const { return iset_; }

    /// Fortran: Ctq5Pdf(Iparton, X, Q) (l.61-95). Q is passed by reference as in the Fortran: Q < Alambda prints
    /// 'Q out of range' and SETS Q = Alambda (the caller's variable). X outside [0,1] is a Fortran STOP (FortranStop).
    /// Iparton outside +-NfMx returns 0 (a warning the first time, SAVEd Warn). Messages go to `log` (may be null).
    double pdf(int iparton, double x, double& q, std::ostream* log = nullptr);
    /// Fortran: PartonX(IPRTN, X, Q) (l.97-196). Updates the SAVEd `First` (the X < Xmin warning is printed once).
    double partonX(int iprtn, double x, double q, std::ostream* log = nullptr);

    // COMMON /CtqPar1/ Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX) as ONE contiguous storage in the Fortran order: PartonX
    // can index XV(-1) (= Al) and UPD(0) (= QL(MXQ)) when x <= XV(0) after its first call (Jx = -1, the clamp to 0 is
    // only done the first time, Ctq5Pdf.f:149-156); the reads are reproduced as the Fortran memory layout makes them.
    double& Al() { return c1_[0]; }
    double Al() const { return c1_[0]; }
    double& XV(int i) { return c1_[static_cast<std::size_t>(kXV + i)]; }     ///< XV(i), i = -1 .. MXX
    double XV(int i) const { return c1_[static_cast<std::size_t>(kXV + i)]; }
    double& QL(int i) { return c1_[static_cast<std::size_t>(kQL + i)]; }     ///< QL(i) = ln(Q_i/Al), i = 0 .. MXQ
    double QL(int i) const { return c1_[static_cast<std::size_t>(kQL + i)]; }
    double& UPD(int i) { return c1_[static_cast<std::size_t>(kUPD + i - 1)]; }   ///< UPD(i), i = 0 .. MXPQX (1-based)
    double UPD(int i) const { return c1_[static_cast<std::size_t>(kUPD + i - 1)]; }
    // /CtqPar2/ Nx, Nt, NfMx; /XQrange/ Qini, Qmax, Xmin; /QCDtable/ Alambda, Nfl, Iorder; /Masstbl/ Amass(6)
    int Nx{0}, Nt{0}, NfMx{0};
    double Qini{0}, Qmax{0}, Xmin{0};
    double Alambda{0};
    int Nfl{0}, Iorder{0};
    std::array<double, 6> Amass{};
    int points_read() const { return npts_read_; }
    int points_expected() const { return npts_; }

    // SAVEd locals
    bool warn{true};    ///< Ctq5Pdf: Data Warn /.true./
    bool first{true};   ///< PartonX: Data First /.true./

private:
    static constexpr int kXV = 1, kQL = kXV + MXX + 1, kUPD = kQL + MXQ + 1;
    std::vector<double> c1_ = std::vector<double>(static_cast<std::size_t>(kUPD + MXPQX), 0.0);
    int iset_{0};
    int npts_{0}, npts_read_{0};
};

}  // namespace sidis
}  // namespace simc

#endif  // SIMC_PHYSICS_SIDIS_CTEQ5PDF_H
