// include/simc/physics/rho/RhoPhysics.h
//
// Exclusive rho0 electroproduction of the Fortran SIMC (task P13.1), bit-for-bit ports:
//   peerho      rho_physics.f:1-402    d5sigma/dEe'dOmega_e'dOmega_rho [ub/MeV/sr^2] (PYTHIA/HERMES-type model)
//   generateRho generate_rho.f:1-128   the rho in 4pi in the photon-nucleon CM (Breit-Wigner mass), boosted to the lab
//   rhoDecay    rho_decay.f:1-160       rho0 -> pi+ pi-: the detected pion, thrown in the rho rest frame
// The reaction plug-in that calls generateRho / rhoDecay (event.f:701-705, 420-421) is task P13.2.
//
// Numerics (D-010, D-024): rho_physics.f, generate_rho.f and rho_decay.f are compiled with -fdefault-real-8, so every
// real literal is a double, EXCEPT the DOUBLE PRECISION literals 1.d+06 (rho_physics.f:280, 300) and 1.d3 (l.384),
// which become REAL(16): `sig219/1.d+06` and `davesig/1.d3` are binary128 divisions rounded back to REAL*8
// (__float128 here). x**2 -> x*x, real exponents -> pow, the Fortran operation order is kept.
//
// The COMMON state these routines read and write:
//   /pfermi/ efer, pfer, pferx..z      -> MainEvent::pfermi (FermiMomentum)
//   /gnrl/   Mh, Mh2, Mh2_final        -> RhoMassState (per-event mutable for rho: the Breit-Wigner rho mass from
//                                         generate_rho, the pion mass after rho_decay; owned by the P13.2 reaction)
//   ntup%rhomass, ntup%rhotheta, ntup%sigcm -> NtupVars
//
// Fortran findings (D-057; P13.F1 ids F1..F22, reference/P13.F1_report.md; docs/FORTRAN_DIFFERENCES.md). Corrected by
// default, kept by fortran_compat:
//   F3  generate_rho.f:80        `Mh2_final = Mh` (a mass, not a mass squared) -> Mh2 (GenerateRhoConfig).
//   F6  rho_physics.f:88-93      peerho recomputes the COMMON efer on shell (the A > 1 branch tests the PION flags) and
//                                writes it back -> the generated efer is used and kept (PeerhoConfig::efer_from_common).
//   F9  rho_decay.f:123-124      beam-left arm, pion on the far side of the beam: yptar with the wrong sign.
//   F20 rho_physics.f:388-389    NaN passes the 'crazy sigma' guard -> 0.
//   F21 rho_decay.f:111-112      acos of a rounding-overflowed cosine gives a NaN yptar that is accepted -> clamped.
// Kept in both modes (expert list): F5 generate_rho.f:72 s without the Fermi terms; F7 rho_physics.f:115, 233 t and
//   tmin with the nominal Mrho2 and main%w/t overwritten in GeV; F8 rho_decay.f:57-64 decay angles w.r.t. the lab axes
//   (not the helicity frame); F10 rho_decay.f:110, 142 xptar = tan(asin(px/p)).
// Faithful: F2 (rho_decay never sets success = .true.; the call after a failed try is the reaction's business, P13.2),
//   F22 (dlen unused, the 'going backwards' test unreachable), the write(222,...) diagnostics of rho_physics.f:252, 262,
//   293, 390 go to the optional `diag` stream (no fort.222 file).
#ifndef SIMC_PHYSICS_RHO_RHOPHYSICS_H
#define SIMC_PHYSICS_RHO_RHOPHYSICS_H

#include <ostream>

#include "simc/core/Event.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace rho {

/// Fortran: the mutable part of COMMON /gnrl/ for doing_rho (core/Event.h HadronMass): generate_rho sets the
/// Breit-Wigner rho mass of the try, rho_decay the pion mass; persistent from one try to the next.
using RhoMassState = HadronMass;

// --------------------------------------------------------------------------------------------- peerho
/// The run constants peerho reads from COMMON blocks.
struct PeerhoConfig {
    double Mtar_struck{0};        ///< targ%Mtar_struck [MeV] (= Mp for rho, dbase.f:316)
    double targ_M{0};             ///< targ%M [MeV] (only for the doing_deutpi/doing_hepi branch)
    bool doing_deutpi{false};     ///< rho_physics.f:89 tests the PION flags (always false in a rho run)
    bool doing_hepi{false};
    /// D-057 F6 (FORTRAN-BUG): false = the Fortran (fortran_compat): efer recomputed on shell and written back to
    /// main.pfermi.efer; true = corrected: the COMMON efer of the generation is used and left unchanged.
    bool efer_from_common{true};
    /// D-057 F20 (FORTRAN-BUG): true = a NaN cross section is set to 0 like the other 'crazy' values.
    bool nan_guard{true};
};

/// The intermediate quantities of one peerho call (for the L1 tests; the Fortran returns only peerho).
struct PeerhoTrace {
    double efer{0}, tcos{0}, tsin{0}, tfcos{0}, tfsin{0}, epsi{0};
    double s{0}, nu_eq{0}, t{0}, tmin{0}, tprime{0};
    double qstar{0}, ppicm{0}, thetacm{0}, phicm{0}, phiqn{0};
    double sig0{0}, R{0}, sigt{0}, cdeltatau{0}, brho{0}, sig219{0}, sig{0};
    double gtpr{0}, factor{0}, jacobian{0}, johnjac{0}, davesig{0};
    bool crazy{false};   ///< sigma_eerho > 1e10 or < 0 -> 0 (l.388-393)
};

/// Fortran: rho_physics.f:peerho (l.1-402). d5sigma/dEe'dOmega_e'dOmega_rho [ub/MeV/sr^2].
/// Reads vertex (Ein, nu, q, Q2, e%E, e%theta, p%P, p%E, up, uq), main%phi_pq and main.pfermi (COMMON /pfermi/).
/// Writes main%w (= sqrt(s) in GeV), main%t (GeV^2), main%pcm, thetacm, phicm, davejac, johnjac, main.ntup.sigcm and,
/// in the Fortran mode (cfg.efer_from_common = false), main.pfermi.efer. The fort.222 diagnostics go to `diag`.
double peerho(const Event& vertex, MainEvent& main, const PeerhoConfig& cfg, std::ostream* diag = nullptr,
              PeerhoTrace* trace = nullptr);

// --------------------------------------------------------------------------------------------- generate_rho
struct GenerateRhoConfig {
    double Mrec_struck{0};          ///< targ%Mrec_struck [MeV] (= Mp, dbase.f:317)
    double spec_p_P{0};             ///< spec%p%P [MeV/c] (central hadron momentum)
    /// D-057 F3 (FORTRAN-BUG): false = the Fortran `Mh2_final = Mh` (fortran_compat); true = corrected Mh2_final = Mh2.
    bool mh2_final_fixed{true};
};

/// What generate_rho does to the event.
struct GenerateRhoResult {
    bool success{false};
    int n_uniform{0};   ///< grnd calls: 0 (beta_cm > 1), 1 (below threshold) or 3
};

/// Fortran: generate_rho.f:generate_rho (l.1-128). Reads vertex (nu, Q2, q, uq) and pf (COMMON /pfermi/ efer, pfer,
/// pferx..z). On success writes vertex%p%delta, P, E, up, xptar (= lab theta, acos(pz/p)) and yptar (= lab phi,
/// atan2(py, px)); gnrl.Mh, Mh2, Mh2_final and ntup_rhomass are written as soon as the mass is drawn (also when the
/// threshold test fails afterwards).
/// Random numbers, in this order: Breit-Wigner u1 (Mh = Mrho + 75.1 tan((2 u1 - 1) atan(1000/150.2))), then
/// phi = 2 pi u2, cos(theta) = 2 u3 - 1 (u2, u3 only above threshold; none when |beta_cm| > 1).
GenerateRhoResult generateRho(Event& vertex, const FermiMomentum& pf, RhoMassState& gnrl, double& ntup_rhomass,
                              const GenerateRhoConfig& cfg, IRandom& rng);

// --------------------------------------------------------------------------------------------- rho_decay
struct RhoDecayConfig {
    int hadron_arm{0};          ///< 1 HMS, 3 HRS-R: beam right; 2 SOS, 4 HRS-L, 5 SHMS: beam left; else STOP
    double spec_p_theta{0};     ///< spec%p%theta [rad]
    double spec_p_phi{0};       ///< spec%p%phi [rad] (physics_angles)
    bool yptar_sign_fixed{true};   ///< D-057 F9 (FORTRAN-BUG): false = the Fortran sign for far-side pions, left arms
    bool acos_clamp{true};         ///< D-057 F21 (FORTRAN-BUG): false = the Fortran acos of an argument > 1 (NaN)
};

struct RhoDecayResult {
    bool success{false};
    int n_uniform{0};     ///< grnd calls: 1 + 2 per accept-reject trial
    int n_trials{0};      ///< accept-reject trials of the decay angle
};

/// Fortran: rho_decay.f:rho_decay (l.1-160). The rho in `orig` (p%P, up, Q2) decays into the detected pion:
/// cos(theta*) from sin^2 + 2 eps R cos^2 (R = 0.33 (Q2/Mrho2)^0.61, accept-reject), phi* flat, E* = rhomass/2.
/// gnrl.Mh2 (the rho mass squared of generate_rho) gives the boost; on success orig%up, xptar, yptar, delta, P, E, theta,
/// phi become the pion's and gnrl = (Mpi, Mpi2, Mpi2). ntup.rhotheta is written after the accept-reject loop; ntup.rhomass
/// is read. success_in is returned unchanged when no failure occurs (D-057 F2: the Fortran never sets success = .true.).
/// Random numbers: phi* = 2 pi u1; then per trial cos(theta*) = 2 u - 1 and norm = (1 + 2 eps R) u'.
/// hadron_arm outside 1..5 -> simc::FortranStop ('Unknown spectrometer setup dude!', rho_decay.f:96-98).
/// The 'Decay pion going backwards - fahgettaboutit!' message (l.264) goes to `log` (if non-null).
RhoDecayResult rhoDecay(Event& orig, double p_spec, double epsilon, bool success_in, RhoMassState& gnrl,
                        NtupVars& ntup, const RhoDecayConfig& cfg, IRandom& rng, std::ostream* log = nullptr);

/// R = sigma_L/sigma_T = 0.33 (Q2/Mrho2)^0.61 (rho_physics.f:250, rho_decay.f:48), Q2 in MeV^2.
double rhoR(double Q2);

}  // namespace rho
}  // namespace simc

#endif  // SIMC_PHYSICS_RHO_RHOPHYSICS_H
