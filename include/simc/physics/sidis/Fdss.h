/**
 * @file
 * @ingroup physics
 * @brief fDSS fragmentation functions (de Florian, Sassot, Stratmann 2007), a bit-for-bit port of fdss.f.
 *
 * fDSS fragmentation functions (de Florian, Sassot, Stratmann, PRD 75 114010 and PRD 76 074033 (2007)), a bit-for-bit
 * port of fdss/fdss.f (task P12.1):
 *   FdssGrid::load   fdss.f:fDSS l.70-139  grid selection and reading (fdss/<set>.GRID, format 9(1PE10.3)), the
 *                                          interpolation arrays XUTOTF ... XSVALF, ARRF (SAVEd in the Fortran)
 *   FdssGrid::eval   fdss.f:fDSS l.58-69, 140-200   z D(z, Q2) of u, ubar, d, dbar, s, sbar, c, b, g for charge IC
 *   ffint            fdss.f:fFINT (l.206-247)       CERNLIB E104 multilinear interpolation, NARG <= 2 here
 *
 * SIMC calls fDSS only for doing_semika (semi_physics.f:496-506): IH = 2 (kaons), IO = 1 (NLO), IC = +1/-1 from
 * doing_hplus, i.e. only fdss/KANLO.GRID is read (copy in data/sidis/). The other grids are supported by load().
 *
 * Numerics (fdss.o, -fdefault-real-8): `implicit real*8 (A-H,O-Z)`, but the D-exponent literals (1.D0, 0.05D0, 1.D5,
 * the Q grid 1.25D0 ...) are REAL(16): (1.D0-X)**4 and **7 are evaluated in binary128 (__subtf3, __multf3, the powi
 * chain x2 = x*x, x3 = x*x2, x4 = x2*x2, x7 = x3*x4), and so are the products fFINT(...) * (1.D0-X)**n * X**p, which
 * are rounded to REAL(8) once, on assignment; in fFINT, FAC*(1.D0-D(I)) is a binary128 product rounded to REAL(8).
 * fFINT itself returns REAL(8) (-ff2c passes it as binary128: exact). X**0.5 and X**0.3 are libm pow (no sqrt).
 */
#ifndef SIMC_PHYSICS_SIDIS_FDSS_H
#define SIMC_PHYSICS_SIDIS_FDSS_H

#include <array>
#include <ostream>
#include <string>
#include <vector>

namespace simc {
namespace sidis {

/// @addtogroup physics
/// @{

/// Fortran: fdss.f:fFINT(NARG, ARG, NENT, ENT, TABLE) (l.206-247). ARG(NARG), NENT(NARG), ENT(sum NENT), TABLE
/// column-major; 0-based pointers to the Fortran arrays (ARG(1) = arg[0]). NARG must be 1 or 2 (the Fortran
/// dimensions D(2), NCOMB(2), IENT(2)).
/// @param narg number of dimensions (1 or 2).
/// @param arg the point.
/// @param nent grid nodes per dimension.
/// @param ent the grid nodes.
/// @param table the table (column-major).
/// @return the interpolated value.
double ffint(int narg, const double* arg, const int* nent, const double* ent, const double* table);

/// The nine fDSS outputs (always z times the distribution).
struct FdssValues {
    double U{0}, UB{0}, D{0}, DB{0}, S{0}, SB{0}, C{0}, B{0}, GL{0};  ///< z times the fragmentation function of u, ubar, d, dbar, s, sbar, c, b, gluon.
};

/// @brief One fDSS grid (hadron IH, order IO) with the SAVEd arrays of fdss.f.
///
/// Read-only after load(); eval() is const and thread-safe.
class FdssGrid {
public:
    static constexpr int NPART = 9, NX = 35, NQ = 24, NARG = 2;  ///< dimensions: partons, x nodes, Q2 nodes, interpolation dimensions.
    /// fdss.f:57-63 DATA QS (REAL(16) literals rounded to REAL(8)) and XB (REAL(8) literals)
    static const std::array<double, NQ> kQS;
    /// The x nodes (fdss.f DATA XB).
    static const std::array<double, NX> kXB;

    /// Grid file of (IH, IO), fdss.f:72-95: 1 pion, 2 kaon, 3 proton, 4 charged hadrons; IO 0 = LO, 1 = NLO.
    /// Returns nullptr for a wrong set (Fortran: ' WRONG SET' STOP).
    /// @param ih hadron (1-4).
    /// @param io order (0 LO, 1 NLO).
    /// @return the file name, or nullptr.
    static const char* gridFile(int ih, int io);

    /// The first-call branch of fDSS (FINI = 0, l.70-139): read <dir>/<gridFile(ih, io)> and build the interpolation
    /// arrays. Returns false and sets `err` for a wrong set, a missing file or a short/malformed record (in the Fortran
    /// a missing file is created empty by OPEN and the READ fails with an end-of-file run-time error).
    /// @param dir the data directory.
    /// @param ih hadron (1-4).
    /// @param io order (0 LO, 1 NLO).
    /// @param[out] err error text (optional, may be nullptr).
    /// @return false on error.
    bool load(const std::string& dir, int ih, int io, std::string* err = nullptr);
    /// @return true after a successful load().
    bool loaded() const { return loaded_; }
    /// @return the hadron of the grid.
    int ih() const { return ih_; }
    /// @return the order of the grid.
    int io() const { return io_; }

    /// fDSS(IH, IC, IO, X, Q2, ...) after the initialisation (l.58-69, 140-200): IC = 1 (+), -1 (-), 0 (average);
    /// another IC is a Fortran STOP (' WRONG CHARGE', FortranStop). The out-of-range warnings (X < 0.05 or > 1,
    /// Q2 < 1 or > 1e5; the Fortran extrapolates) are written to `log`.
    /// @param ic charge (1, -1, 0).
    /// @param x z of the hadron.
    /// @param q2 Q2 [GeV^2].
    /// @param log stream for the warnings (not owned), or nullptr.
    /// @return the nine values.
    /// @throws FortranStop for a wrong charge.
    FdssValues eval(int ic, double x, double q2, std::ostream* log = nullptr) const;

    // SAVEd arrays of fDSS (column-major (NX, NQ): element (ix, iq) at (iq-1)*NX + (ix-1))
    std::array<std::vector<double>, 9> F;   ///< XUTOTF XDTOTF XSTOTF XCTOTF XBTOTF XGF XUVALF XDVALF XSVALF
    /// Interpolation nodes of fFINT (x nodes, then log Q2 nodes).
    std::array<double, NX + NQ> ARRF{};
    /// Number of nodes per dimension of fFINT.
    std::array<int, NARG> NA{};
    /// PARTON(9, NQ, NX-1) as read (element (i, n, m) at ((m-1)*NQ + (n-1))*9 + (i-1))
    std::vector<double> parton;

private:
    bool loaded_{false};
    int ih_{0}, io_{0};
};

/// @}

}  // namespace sidis
}  // namespace simc

#endif  // SIMC_PHYSICS_SIDIS_FDSS_H
