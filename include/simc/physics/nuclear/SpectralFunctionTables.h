// include/simc/physics/nuclear/SpectralFunctionTables.h
//
// Kernels of the A(e,e'p) spectral functions (task P8.1). Fortran-structured (1-based arrays, SAVEd locals as
// members); the plug-ins of include/simc/physics/nuclear/SpectralFunctions.h wrap them behind ISpectralFunction.
//
//   SfTable           COMMON /sftable/ (sf_lookup.inc) and the routines of sf_lookup.f:
//                       load()        sf_lookup_init  (l.1-81)    read + normalise a Benhar table
//                       lookup()      sf_lookup       (l.99-190)  bin probability at (Em, Pm), bilinear (as coded)
//                       lookupDiff()  sf_lookup_diff  (l.85-95)   differential spectral function
//                       generateEm()  generate_em     (l.195-245) Em at fixed Pm (2 grnd); used by A(e,e'pi/K) only
//   TheorySfWeight    event.f:complete_main (l.1402-1430): SF_weight of D(e,e'p) / A(e,e'p) from the theory-file
//                     momentum distributions (init.f:theory_init, RunSetup::theory).
//
// fortran_compat (D-031, docs/FORTRAN_DIFFERENCES.md):
//   * FORTRAN-BUG(D-033 F6): sf_lookup leaves Em1, Em2, sf1, sf2 unassigned for Em == Emval(numEm) (and NaN Em) and
//     then uses the SAVEd values of the PREVIOUS call (NaN on the first call of a process). Corrected: Em ==
//     Emval(numEm) uses the last interval; no state is carried between calls (NaN Em -> 0).
//   * FORTRAN-BUG(D-033 F3): sf_lookup_diff divides the interpolated BIN probability by the event's own Pm**2 and by
//     hard-coded widths 5 MeV x 20 MeV/c (and pi = 3.1415926535): +Inf/NaN at Pm = 0, a 1/Pm**2 divergence below the
//     first node, and wrong widths for the 3mod/4 files (dPm 19.733, dEm 1/5/10). Corrected: the differential
//     spectral function of every node, sfval/(4 pi Pmval**2 dPm dEm) with the node's widths from the file, is
//     interpolated with the same (corrected) bilinear scheme.
//   * "logsf" is a linear interpolation (no log is taken): documented, kept (FORTRAN-QUIRK D-033 F6).
#ifndef SIMC_PHYSICS_NUCLEAR_SPECTRALFUNCTIONTABLES_H
#define SIMC_PHYSICS_NUCLEAR_SPECTRALFUNCTIONTABLES_H

#include <string>
#include <vector>

#include "simc/core/RunSetup.h"

namespace simc {

class IRandom;

namespace nuclear {

/// Fortran: sf_lookup.inc COMMON /sftable/ + the SAVEd locals of sf_lookup.f (-fno-automatic).
class SfTable {
public:
    static constexpr int nPmmax = 100;   ///< sf_lookup.inc
    static constexpr int nEmmax = 200;

    /// Fortran: sf_lookup_init(infile, protonflag) (sf_lookup.f:1-81). Reads the table, keeps the proton (SFp) or
    /// neutron (SFn) column, normalises it to 1 (sftotnorm). The 'Bad Pm/Em value' and 'Uncorrected S.F. has
    /// normalization=' lines go to `msgs` (gfortran list-directed text). Returns false (and `err`) if the file cannot
    /// be read (the Fortran aborts).
    bool load(const std::string& path, bool protonflag, std::vector<std::string>* msgs, std::string* err);

    /// Fortran: sf_lookup(Em, Pm, SF) (sf_lookup.f:99-190). Throws FortranStop for the two 'stop' checks.
    double lookup(double Em, double Pm);
    /// Fortran: sf_lookup_diff(Em, Pm, SFd) (sf_lookup.f:85-95). Corrected unless fortran_compat (F3, see above).
    double lookupDiff(double Em, double Pm);
    /// Fortran: generate_em(Pm, Em) (sf_lookup.f:195-245): exactly 2 rng.uniform().
    double generateEm(double Pm, IRandom& rng);

    void setFortranCompat(bool c) { compat_ = c; }
    bool fortranCompat() const { return compat_; }
    bool loaded() const { return numPm_ > 0; }

    // ---- /sftable/ (1-based like the Fortran) ----
    int numPm() const { return numPm_; }
    int numEm() const { return numEm_; }
    double Pmval(int i) const { return Pmval_[i]; }
    double dPm(int i) const { return dPm_[i]; }
    double Emval(int i) const { return Emval_[i]; }
    double dEm(int i) const { return dEm_[i]; }
    double sfval(int iEm, int iPm) const { return sfval_[idx(iEm, iPm)]; }
    double sfnorm(int iPm) const { return sfnorm_[iPm]; }
    double sftotnorm() const { return sftotnorm_; }
    const std::string& file() const { return file_; }

private:
    int idx(int iEm, int iPm) const { return iPm * (nEmmax + 1) + iEm; }
    /// The interpolation of sf_lookup on table `t` (sfval_ or dens_). `stale` = use/keep the SAVEd Em1..sf2.
    double interpolate(const std::vector<double>& t, double Em, double Pm, bool stale);

    bool compat_{false};
    std::string file_;
    int numPm_{0}, numEm_{0};
    std::vector<double> Pmval_, dPm_, Emval_, dEm_, sfnorm_;
    std::vector<double> sfval_;   ///< sfval(iEm, iPm), Fortran column-major, 1-based
    std::vector<double> dens_;    ///< corrected mode: sfval / (4 pi Pmval^2 dPm dEm)
    double sftotnorm_{0};
    // SAVEd locals of sf_lookup (zero-initialised static storage)
    double Em1_{0}, Em2_{0}, sf1_{0}, sf2_{0};
    // SAVEd locals of generate_em: x(nEmmax), y(nEmmax)
    std::vector<double> gx_, gy_;
};

/// Fortran: event.f:complete_main (l.1402-1430), the theory-file branch (doing_deuterium, or doing_heavy without the
/// Benhar SF): SF_weight = sum over the shells i of rho_i(Pm) * nprot_i, with rho_i the linear interpolation of the
/// theory(i, n) table (nint() index, forced into [1, n-1]) and, for doing_heavy, times the Lorentzian Em distribution
/// width/((Em-Em_i)**2+width**2)/pi/Em_int_i (0 below E_Fermi).
double TheorySfWeight(const TheoryData& t, bool doing_heavy, double Em, double Pm);

}  // namespace nuclear
}  // namespace simc

#endif  // SIMC_PHYSICS_NUCLEAR_SPECTRALFUNCTIONTABLES_H
