// include/simc/physics/delta/DeltaModel.h
//
// The doing_delta "cross section" of the Fortran SIMC (task P11.1): physics_delta.f:peedelta, called by
// event.f:complete_main (l.1511-1513) for H/D/A(e,e'p)pi (the proton is detected, the pion is missing):
//     main%sigcc       = peedelta(vertex, main)   ("Need new xsec model.")
//     main%sigcc_recon = 1.0
// and the IPhysicsModel plug-in "delta".
//
// IMPORTANT (D-048 F4, FORTRAN-DIFF-CANDIDATE, kept in both modes): peedelta returns 1.0*jacobian*(gtpr*fac)
// (physics_delta.f:132); the line with the cross section (l.130, sigma_eepi*jacobian*(gtpr*fac)) is commented out. The
// Fortran Delta weight is therefore phase space x virtual-photon flux x the transform_to_cm jacobian dt dphi_cm ->
// dOmega_lab of the proton: it contains NO cross section. SIMC++ keeps this placeholder (no model can be derived from the
// Fortran); the question whether a real Delta cross section should be offered is open with the user/an expert.
//
// What peedelta computes on the way (physics_delta.f:50-128): transform_to_cm with the PROTON as the 'hadron'
// (main%thetacm, phicm, pcm, davejac, johnjac), tfcos, s, main%wcm = sqrt(s), k_eq, and
//     ntup%sigcm1 = ntup%sigcm = sig_blok(thetacm, phicm, main%t/1.d6, vertex%q2/1.d6, s/1.d6, main%epsilon,
//                                         targ%Mtar_struck/1000., which_pion)
// (the Brauel p(e,e'pi+)n fit, physics_pion.f:400-465; pi- branch for which_pion 1/3/11). It is written to the ntuple
// column `sigcm` only; it does not enter the weight.
//   FORTRAN-BUG(D-048 F3): under -fdefault-real-8 the literal 1.d6 is REAL(16): the three quotients are binary128
//   temporaries whose addresses are passed to the REAL*8 dummies of sig_blok, which reads their low 8 bytes (the low 64
//   mantissa bits) as doubles. ntup%sigcm is meaningless (0 in ~65 % of events, negative in ~19 %, NaN/Inf in ~0.4 %).
//   fortran_compat reproduces it (DeltaConfig::sigblok_binary128_quirk); the corrected default passes the doubles.
//   FORTRAN-QUIRK: physics_delta.f:74 `tfcos-1. .lt. 1.d-8` is a binary128 comparison (kept in both modes).
//   FORTRAN-QUIRK(D-048 F9): ntup%sigcm2 is not touched by peedelta (stale; not a column of the delta ntuple).
//
// JSON config: {"verbose": bool} (no data files).
#ifndef SIMC_PHYSICS_DELTA_DELTAMODEL_H
#define SIMC_PHYSICS_DELTA_DELTAMODEL_H

#include <ostream>
#include <string>

#include "simc/core/Event.h"
#include "simc/stages/IPhysicsModel.h"

namespace simc {
namespace delta {

/// The run constants peedelta reads from COMMON blocks.
struct DeltaConfig {
    int which_pion{0};                   ///< dbase.f which_pion (selects the pi- branch of sig_blok)
    double Mtar_struck{0};               ///< targ%Mtar_struck [MeV] (= Mp)
    bool sigblok_binary128_quirk{true};  ///< D-048 F3: true = the Fortran (fortran_compat), false = corrected
};

/// The intermediate quantities of one peedelta call (for the L1 tests; peedelta itself returns only the value).
struct PeedeltaTrace {
    double tfcos{0}, s{0}, k_eq{0}, fac{0}, gtpr{0};
    double t_arg{0}, q2_arg{0}, s_arg{0};   ///< the t, q2_gev, s_gev values sig_blok actually used
};

/// Fortran: physics_delta.f:peedelta (l.1-135). Reads vertex (Ein, nu, q, Q2, e%E, p%P, p%E, up, uq), main%t,
/// main%epsilon, main%phi_pq and main.pfermi (COMMON /pfermi/). Writes main%thetacm, phicm, pcm, davejac, johnjac, wcm,
/// main.ntup.sigcm1 and main.ntup.sigcm. Returns 1.0*jacobian*(gtpr*fac) [D-048 F4: no cross section].
double peedelta(const Event& vertex, MainEvent& main, const DeltaConfig& cfg, std::ostream* log = nullptr,
                PeedeltaTrace* trace = nullptr);

/// What sig_blok reads when the Fortran passes the binary128 quotient x/1.d6 to a REAL*8 dummy: the low 8 bytes of the
/// binary128 value (x86-64, little endian) reinterpreted as a double (D-048 F3).
double binary128LowBytes(double x_over_1d6_numerator);

class DeltaModel : public IPhysicsModel {
public:
    std::string name() const override { return "delta"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    /// peedelta(vertex, main) on a copy of main.
    double sigcc(const Event& vertex, const MainEvent& main) override;
    /// main%sigcc = peedelta(vertex, main): also writes main%thetacm ... and ntup%sigcm1, sigcm.
    double sigccMain(const Event& vertex, MainEvent& main) override;
    /// event.f:1513 main%sigcc_recon = 1.0
    double sigccRecon(const Event&, const MainEvent&) override { return 1.0; }

    DeltaConfig& config() { return cfg_; }
    const DeltaConfig& config() const { return cfg_; }

private:
    bool verbose_{true};
    DeltaConfig cfg_;
};

}  // namespace delta
}  // namespace simc

#endif  // SIMC_PHYSICS_DELTA_DELTAMODEL_H
