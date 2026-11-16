/**
 * @file
 * @ingroup physics
 * @brief CTEQ5 parton distribution functions, a bit-for-bit port of cteq5/Ctq5Pdf.f.
 *
 * CTEQ5 parton distribution functions, a bit-for-bit port of cteq5/Ctq5Pdf.f (task P12.1):
 *   Cteq5Pdf::load     Ctq5Pdf.f:SetCtq5 + ReadTbl (l.199-276)   table file -> COMMON /CtqPar1/ /CtqPar2/ /XQrange/
 *                                                                 /QCDtable/ /Masstbl/
 *   Cteq5Pdf::pdf      Ctq5Pdf.f:Ctq5Pdf (l.61-95)              x f(x, Q) of one parton, negative values -> 0
 *   Cteq5Pdf::partonX  Ctq5Pdf.f:PartonX (l.97-196)             quadratic (M = 2) interpolation in x and ln(Q/Lambda)
 *   polint             Ctq5Pdf.f:POLINT (l.297-334)             Numerical Recipes polynomial interpolation
 *
 * SIMC calls SetCtq5(iset = 1) (semi_physics.f:121, 300 and semi_dilution.f:79): only cteq5/cteq5m.tbl (CTEQ5M, MSbar)
 * is read; a copy is in data/sidis/. The other sets are supported by load() (same file format) but never used.
 *
 * Numerics. Ctq5Pdf.f has its `Implicit Double Precision` lines commented out, so under -fdefault-real-8 every variable is
 * REAL(8) (implicit typing, I-N integer); the literals 0D0/1D0 are REAL(16) but are only compared with or assigned to
 * REAL(8) values (exact). Ctq5Pdf.o contains no binary128 arithmetic (only __extenddftf2/__trunctfdf2 of the -ff2c
 * return convention, exact), so the port is plain double. The table values are read list-directed (correctly rounded,
 * strtod).
 *
 * State. The SAVEd locals `Warn` (Ctq5Pdf) and `First` (PartonX) and the COMMON blocks are members: one Cteq5Pdf object
 * is "one Fortran process" (shared by peepiX and semi_dilution, as the COMMON is).
 */
#ifndef SIMC_PHYSICS_SIDIS_CTEQ5PDF_H
#define SIMC_PHYSICS_SIDIS_CTEQ5PDF_H

#include <array>
#include <ostream>
#include <string>
#include <vector>

namespace simc {
namespace sidis {

/// @addtogroup physics
/// @{

/// Fortran: Ctq5Pdf.f:POLINT(XA, YA, N, X, Y, DY) (l.297-334), N <= 10. xa/ya point to N values (XA(1) = xa[0]).
/// Returns Y; DY in `dy`. DEN = 0 (two equal abscissae) is a Fortran PAUSE: throws simc::FortranStop.
/// @param xa abscissae.
/// @param ya ordinates.
/// @param n number of points (<= 10).
/// @param x the point.
/// @param[out] dy error estimate.
/// @return the interpolated value.
double polint(const double* xa, const double* ya, int n, double x, double& dy);

/// @brief One CTEQ5 set with the COMMON blocks and SAVEd locals of Ctq5Pdf.f. Not thread-safe (SAVEd warnings).
class Cteq5Pdf {
public:
    static constexpr int MXX = 105, MXQ = 25, MXF = 6;  ///< table dimensions (Ctq5Pdf.f PARAMETERs).
    static constexpr int MXPQX = (MXF * 2 + 2) * MXQ * MXX;  ///< size of UPD.
    static constexpr int M = 2, M1 = M + 1;  ///< interpolation order and points.
    /// Ctq5Pdf.f:207-212 Flnm(1..9) without the directory 'cteq5/'
    static const std::array<const char*, 9> kFiles;

    /// SetCtq5(iset) + ReadTbl: reads <dir>/<kFiles[iset-1]>. Returns false and sets `err` if iset is out of range
    /// (Fortran: 'Invalid Iset number in SetCtq5' STOP), if the file is missing (Fortran: STOP l.226) or malformed.
    /// The table read stops silently at an end of file (IOSTAT=IRET, l.272): `points_read()` tells how many were read.
    /// @param dir the data directory.
    /// @param iset the set number (1-9).
    /// @param[out] err error text (optional, may be nullptr).
    /// @return false on error.
    bool load(const std::string& dir, int iset, std::string* err = nullptr);
    /// @return true after a successful load().
    bool loaded() const { return iset_ > 0; }
    /// @return the loaded set number (0 if none).
    int iset() const { return iset_; }

    /// Fortran: Ctq5Pdf(Iparton, X, Q) (l.61-95). Q is passed by reference as in the Fortran: Q < Alambda prints
    /// 'Q out of range' and SETS Q = Alambda (the caller's variable). X outside [0,1] is a Fortran STOP (FortranStop).
    /// Iparton outside +-NfMx returns 0 (a warning the first time, SAVEd Warn). Messages go to `log` (may be null).
    /// @param iparton parton (-NfMx..NfMx, 0 = gluon).
    /// @param x Bjorken x.
    /// @param[in,out] q scale Q [GeV].
    /// @param log stream for the messages (not owned), or nullptr.
    /// @return the parton distribution.
    /// @throws FortranStop for x outside [0,1].
    double pdf(int iparton, double x, double& q, std::ostream* log = nullptr);
    /// Fortran: PartonX(IPRTN, X, Q) (l.97-196). Updates the SAVEd `First` (the X < Xmin warning is printed once).
    /// @param iprtn parton.
    /// @param x Bjorken x.
    /// @param q scale Q [GeV].
    /// @param log stream for the warning (not owned), or nullptr.
    /// @return x times the distribution as tabulated.
    double partonX(int iprtn, double x, double q, std::ostream* log = nullptr);

    // COMMON /CtqPar1/ Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX) as ONE contiguous storage in the Fortran order: PartonX
    // can index XV(-1) (= Al) and UPD(0) (= QL(MXQ)) when x <= XV(0) after its first call (Jx = -1, the clamp to 0 is
    // only done the first time, Ctq5Pdf.f:149-156); the reads are reproduced as the Fortran memory layout makes them.
    /// @return Al (modifiable).
    double& Al() { return c1_[0]; }
    /// @return Al.
    double Al() const { return c1_[0]; }
    /// @param i index (-1 .. MXX).
    /// @return XV(i) (modifiable).
    double& XV(int i) { return c1_[static_cast<std::size_t>(kXV + i)]; }     ///< XV(i), i = -1 .. MXX
    /// @param i index (-1 .. MXX).
    /// @return XV(i).
    double XV(int i) const { return c1_[static_cast<std::size_t>(kXV + i)]; }
    /// @param i index (0 .. MXQ).
    /// @return QL(i) (modifiable).
    double& QL(int i) { return c1_[static_cast<std::size_t>(kQL + i)]; }     ///< QL(i) = ln(Q_i/Al), i = 0 .. MXQ
    /// @param i index (0 .. MXQ).
    /// @return QL(i).
    double QL(int i) const { return c1_[static_cast<std::size_t>(kQL + i)]; }
    /// @param i index (0 .. MXPQX).
    /// @return UPD(i) (modifiable).
    double& UPD(int i) { return c1_[static_cast<std::size_t>(kUPD + i - 1)]; }   ///< UPD(i), i = 0 .. MXPQX (1-based)
    /// @param i index (0 .. MXPQX).
    /// @return UPD(i).
    double UPD(int i) const { return c1_[static_cast<std::size_t>(kUPD + i - 1)]; }
    // /CtqPar2/ Nx, Nt, NfMx; /XQrange/ Qini, Qmax, Xmin; /QCDtable/ Alambda, Nfl, Iorder; /Masstbl/ Amass(6)
    int Nx{0}, Nt{0}, NfMx{0};  ///< /CtqPar2/ number of x and Q points, number of flavours.
    double Qini{0}, Qmax{0}, Xmin{0};  ///< /XQrange/ Q range [GeV] and minimum x.
    double Alambda{0};  ///< /QCDtable/ Lambda_QCD [GeV].
    int Nfl{0}, Iorder{0};  ///< /QCDtable/ number of flavours and order.
    std::array<double, 6> Amass{};  ///< /Masstbl/ quark masses [GeV].
    /// @return the number of table points read.
    int points_read() const { return npts_read_; }
    /// @return the number of table points expected.
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

/// @}

}  // namespace sidis
}  // namespace simc

#endif  // SIMC_PHYSICS_SIDIS_CTEQ5PDF_H
