/**
 * @file
 * @ingroup core
 * @brief The event data model: kinematics of one try (simc::Event) and everything else the try produces
 *        (simc::MainEvent).
 *
 * A try of the Monte Carlo fills three simc::Event values (the run keeps them in one simc::driver::EventState):
 *  - `vertex`: the kinematics generated at the interaction vertex (before radiation and energy loss);
 *  - `orig`:   the same event after radiation, energy loss and Coulomb corrections, i.e. what enters the spectrometers;
 *  - `recon`:  the kinematics reconstructed from the spectrometer tracks (with the energy-loss correction applied).
 * and one simc::MainEvent with the weights, the target position, the spectrometer and focal-plane coordinates and the
 * ntuple-only values.
 *
 * Units are the internal units of the framework (simc::units): MeV, MeV/c, MeV^2, cm, rad; delta in %; xptar/yptar and
 * the focal-plane dx/dy are slopes (dx/dz, dy/dz). The structs are plain aggregates with zero defaults; a try
 * overwrites only what it computes, so values that a try does not rewrite keep those of the previous try (defined
 * behaviour, as in the original program). Not thread-safe: each worker of the parallel loop has its own events.
 *
 * @note Origin: Fortran SIMC structures.inc (derived types event, event_main, arm, ...), field by field. The COMMON
 *       blocks /ntuple_variables/, /pfermi/ and the per-try part of /gnrl/ (simulate.inc) are carried per event in
 *       MainEvent::ntup, MainEvent::pfermi and MainEvent::gnrl instead of globals.
 */
#ifndef SIMC_CORE_EVENT_H
#define SIMC_CORE_EVENT_H

#include <array>

namespace simc {

/** @addtogroup core
 *  @{ */

/**
 * @brief A 3-vector in the SIMC lab frame (used for unit direction vectors).
 * @note Origin: Fortran SIMC structures.inc type Cartesian.
 */
struct Cartesian {
    /// Lab-frame components (z along the beam; the arm at phi = pi/2 lies at +y, the one at 3 pi/2 at -y); unit
    /// vectors are normalised.
    double x{0}, y{0} /**< @copydoc x */, z{0} /**< @copydoc x */;
};

/**
 * @brief Full description of one particle (electron or hadron arm) of an simc::Event.
 * @note Origin: Fortran SIMC structures.inc type arm_full.
 */
struct ArmFull {
    /// delta = (P - P_central)/P_central [%]; xptar, yptar: target slopes in the spectrometer frame [dx/dz, dy/dz];
    /// z: vertex position along the spectrometer axis ("ytar" coordinate) [cm].
    double delta{0}, xptar{0} /**< @copydoc delta */, yptar{0} /**< @copydoc delta */, z{0} /**< @copydoc delta */;
    /// theta, phi: lab polar and azimuthal angles [rad]; E: total energy [MeV]; P: momentum [MeV/c].
    double theta{0}, phi{0} /**< @copydoc theta */, E{0} /**< @copydoc theta */, P{0} /**< @copydoc theta */;
};

/**
 * @brief The kinematics of one event: beam, both arms, momentum transfer and the derived missing quantities.
 *
 * Used three times per try (vertex, orig, recon; see the file description). The reaction stage fills the vertex
 * event, the radiation and the event completion the orig event, the reconstruction the recon event.
 * @note Origin: Fortran SIMC structures.inc type event.
 */
struct Event {
    /// Incident beam energy at the vertex [MeV].
    double Ein{0};
    /// Missing energy Em [MeV] and missing momentum Pm [MeV/c] of the reaction (definition per reaction family).
    double Em{0}, Pm{0} /**< @copydoc Em */;
    /// Missing energy nu + M_target - E_hadron [MeV] and magnitude of the missing momentum vector [MeV/c].
    double Emiss{0}, Pmiss{0} /**< @copydoc Emiss */;
    /// Missing-momentum vector p_hadron - q in the lab frame [MeV/c].
    double Pmx{0}, Pmy{0} /**< @copydoc Pmx */, Pmz{0} /**< @copydoc Pmx */;
    /// Missing-momentum components parallel to q, perpendicular to q (in plane) and out of the scattering plane [MeV/c].
    double PmPar{0}, PmPer{0} /**< @copydoc PmPar */, PmOop{0} /**< @copydoc PmPar */;
    /// Energy transfer nu [MeV], |q| [MeV/c], Q^2 [MeV^2], recoil kinetic energy Trec [MeV], invariant mass W [MeV],
    /// recoil mass Mrec [MeV].
    double nu{0}, q{0} /**< @copydoc nu */, Q2{0} /**< @copydoc nu */, Trec{0} /**< @copydoc nu */, W{0} /**< @copydoc nu */, Mrec{0} /**< @copydoc nu */;
    /// Virtual-photon polarisation epsilon; angle between hadron and q (theta_pq) [rad]; angle between the target
    /// polarisation and q (theta_tarq) [rad]; azimuth of the hadron around q (phi_pq) [rad]; azimuth of the target
    /// polarisation (phi_targ) [rad].
    double epsilon{0}, theta_pq{0} /**< @copydoc epsilon */, theta_tarq{0} /**< @copydoc epsilon */, phi_pq{0} /**< @copydoc epsilon */, phi_targ{0} /**< @copydoc epsilon */;
    /// Polarised-target angles [rad]: beta (target polarisation in the (q, hadron) frame), phi_s = phi_pq - phi_targ
    /// (Sivers) and phi_c = phi_pq + phi_targ (Collins), each folded into [0, 2 pi).
    double beta{0}, phi_s{0} /**< @copydoc beta */, phi_c{0} /**< @copydoc beta */;
    /// SIDIS variables: hadron energy fraction z, transverse momentum squared pt2 [MeV^2], Bjorken x.
    double zhad{0}, pt2{0} /**< @copydoc zhad */, xbj{0} /**< @copydoc zhad */;
    /// The electron arm (e) and the hadron arm (p).
    ArmFull e, p /**< @copydoc e */;
    /// Unit vectors of the scattered electron (ue), the hadron (up) and the momentum transfer (uq) in the lab.
    Cartesian ue, up /**< @copydoc ue */, uq /**< @copydoc ue */;
};

/**
 * @brief The interaction point of a try and the energy losses on the way in and out of the target.
 * @note Origin: Fortran SIMC structures.inc type event_target.
 */
struct EventTarget {
    /// Vertex position x, y, z [cm] (lab) and the raster offsets rasterx, rastery [cm].
    double x{0}, y{0} /**< @copydoc x */, z{0} /**< @copydoc x */, rasterx{0} /**< @copydoc x */, rastery{0} /**< @copydoc x */;
    std::array<double, 3> teff{};    ///< Material traversed [radiation lengths]: [0] beam, [1] electron, [2] hadron (Fortran teff(1:3))
    std::array<double, 3> Eloss{};   ///< Energy loss [MeV]: [0] beam, [1] electron, [2] hadron (Fortran Eloss(1:3))
    /// Coulomb energy shift of the electron in the field of the nucleus [MeV].
    double Coulomb{0};
};

/**
 * @brief Spectrometer target coordinates of one arm (delta, slopes, ytar).
 * @note Origin: Fortran SIMC structures.inc types arm / arm2.
 */
struct ArmSP {
    /// delta [%], target slopes yptar, xptar [dy/dz, dx/dz], target coordinate z ("ytar") [cm].
    double delta{0}, yptar{0} /**< @copydoc delta */, xptar{0} /**< @copydoc delta */, z{0} /**< @copydoc delta */;
};

/**
 * @brief Spectrometer target coordinates of both arms.
 * @note Origin: Fortran SIMC structures.inc type double_arm.
 */
struct DoubleArm {
    /// Electron arm (e) and hadron arm (p).
    ArmSP e, p /**< @copydoc e */;
};

/**
 * @brief Focal-plane track of one arm.
 * @note Origin: Fortran SIMC structures.inc type arm_FP.
 */
struct ArmFP {
    /// Focal-plane position x, y [cm], slopes dx, dy [dx/dz, dy/dz], path length through the spectrometer [cm].
    double x{0}, dx{0} /**< @copydoc x */, y{0} /**< @copydoc x */, dy{0} /**< @copydoc x */, path{0} /**< @copydoc x */;
};

/**
 * @brief Focal-plane tracks of both arms.
 * @note Origin: Fortran SIMC structures.inc type double_arm_FP.
 */
struct DoubleArmFP {
    /// Electron arm (e) and hadron arm (p).
    ArmFP e, p /**< @copydoc e */;
};

/**
 * @brief Values of a try that only the outputs (ntuple, histograms) read.
 * @note Origin: Fortran SIMC simulate.inc type ntupvars + decdist (l.161-183); only the fields used by the outputs.
 */
struct NtupVars {
    double radphot{0};   ///< total radiated photon energy (MeV)
    /// Which tail radiated: 0 none, 1 incoming electron, 2 scattered electron, 3 hadron, 12 forced two tails.
    double radarm{0};
    /// Resolution factor of the reaction (reaction-specific, written as the ntuple column "resfac").
    double resfac{0};
    /// Model cross section in the centre-of-mass frame of the reaction (units of the model, e.g. [ub/MeV^2/rad]).
    double sigcm{0};
    /// Relative momentum of the recoil system (krel) [MeV/c] (reaction-specific; may keep a previous value).
    double krel{0};
    /// Missing mass on the nucleon (mm) and on the nucleus (mmA) [MeV], Mandelstam t [MeV^2].
    double mm{0}, mmA{0} /**< @copydoc mm */, t{0} /**< @copydoc mm */;
    double decdist{0};   ///< decay distance (cm), Fortran decdist
    double Mh2_final{0}; ///< hadron mass^2 after decay (MeV^2)
    double sigcm1{0};    ///< pion/kaon: first model cross section (physics_pion.f:126 sig_param_2021) [ub/MeV^2/rad]
    double sigcm2{0};    ///< pion/kaon: second model (physics_pion.f:148 MAID, W < 2 GeV only) [ub/MeV^2/rad]
    double dilu{0};      ///< SIDIS dilution (event.f:1525 ntup%dilu = 1.0; semi_dilution is dead code). Not an ntuple column
    double xfermi{0};    ///< SIDIS with Fermi motion: x_Bjorken = Q2/(2 P.q) (semi_physics.f:250, ntuple column 55)
    /// pi0 photons (pizero_decay.f, event.f:899-901): (E, px, py, pz) in the lab [MeV] (task P9.2)
    std::array<double, 4> gamma1{}, gamma2{} /**< @copydoc gamma1 */;
    /// pi0 photon positions at the calorimeter [cm] (simc.f:1500-1560, mc_calo only; never written with HMS/SHMS)
    double xcal_gamma1{0}, ycal_gamma1{0} /**< @copydoc xcal_gamma1 */, xcal_gamma2{0} /**< @copydoc xcal_gamma1 */, ycal_gamma2{0} /**< @copydoc xcal_gamma1 */;
    /// rho: the Breit-Wigner rho mass of the try [MeV] (generate_rho.f:79 ntup%rhomass = Mh; ntuple column Mrho) and
    /// the decay polar angle in the rho rest frame [rad] (rho_decay.f:66 ntup%rhotheta; column Thrho). Task P13.1.
    double rhomass{0}, rhotheta{0} /**< @copydoc rhomass */;
};

/**
 * @brief Fermi motion of the struck nucleon of this try.
 *
 * pfer [MeV/c] is the momentum, (pferx, pfery, pferz) its unit vector in the lab, efer [MeV] its energy (the struck
 * nucleon mass of the target without Fermi motion). Read by the pion/kaon cross sections and the pion ntuple.
 * @note Origin: Fortran SIMC COMMON /pfermi/ (simulate.inc l.212-217), filled by event.f:generate (l.330-371); read by
 *       physics_pion.f:peepi, jacobians.f:transform_to_cm and results_write.f:141-149.
 */
struct FermiMomentum {
    /// Momentum [MeV/c], its lab unit vector (pferx, pfery, pferz), and the nucleon energy efer [MeV].
    double pfer{0}, pferx{0} /**< @copydoc pfer */, pfery{0} /**< @copydoc pfer */, pferz{0} /**< @copydoc pfer */, efer{0} /**< @copydoc pfer */;
};

/**
 * @brief Hadron mass of the try when it changes per try (rho production).
 *
 * Mh, Mh2 = the hadron mass [MeV] and its square [MeV^2]; Mh2_final the mass squared after the decay [MeV^2]. The set-up
 * fixes them once for every reaction except rho production, where the Breit-Wigner rho mass of the try and the decay
 * pion change them per try. They persist across tries.
 * @note Origin: Fortran SIMC the mutable part of COMMON /gnrl/ (simulate.inc:42, 91), set by dbase.f, generate_rho
 *       and rho_decay (task P13.2).
 */
struct HadronMass {
    /// Hadron mass [MeV], its square [MeV^2], the final-state mass squared after the decay [MeV^2].
    double Mh{0}, Mh2{0} /**< @copydoc Mh */, Mh2_final{0} /**< @copydoc Mh */;
};

/**
 * @brief Everything a try produces besides the three event kinematics: weights, target, spectrometer and focal-plane
 *        coordinates, ntuple values.
 * @note Origin: Fortran SIMC structures.inc type event_main.
 */
struct MainEvent {
    /// Event weight (cross section x jacobians x spectral-function and generation weights; normalised at the end of
    /// the run), spectral-function weight, generation weight (phase-space volume factor), jacobian of the generation
    /// variables.
    double weight{0}, SF_weight{0} /**< @copydoc weight */, gen_weight{0} /**< @copydoc weight */, jacobian{0} /**< @copydoc weight */;
    /// Shift of the vertex beam energy from its run average (Coulomb + energy loss) and of the scattered-electron
    /// energy (Coulomb) [MeV]; used to compare generated energies with the generation limits.
    double Ein_shift{0}, Ee_shift{0} /**< @copydoc Ein_shift */;
    /// Cross section of the try (vertex kinematics), of the reconstructed kinematics, and at the central kinematics
    /// (sigcent is never assigned: 0) [units of the physics model].
    double sigcc{0}, sigcc_recon{0} /**< @copydoc sigcc */, sigcent{0} /**< @copydoc sigcc */;
    /// Copies of the corresponding simc::Event quantities for the outputs (epsilon, angles [rad], beta).
    double epsilon{0}, theta_pq{0} /**< @copydoc epsilon */, theta_tarq{0} /**< @copydoc epsilon */, phi_pq{0} /**< @copydoc epsilon */, phi_targ{0} /**< @copydoc epsilon */, beta{0} /**< @copydoc epsilon */;
    /// Meson production: invariant mass w [MeV], Mandelstam t and its minimum tmin [MeV^2], Q^2 [MeV^2].
    double w{0}, t{0} /**< @copydoc w */, tmin{0} /**< @copydoc w */, q2{0} /**< @copydoc w */;
    /// Meson production in the centre-of-mass frame: momentum pcm [MeV/c], angles thetacm, phicm [rad], energy wcm
    /// [MeV].
    double pcm{0}, thetacm{0} /**< @copydoc pcm */, phicm{0} /**< @copydoc pcm */, wcm{0} /**< @copydoc pcm */;
    /// Reaction jacobians written to the ntuple (davejac: the SIDIS/meson jacobian; johnjac: never assigned, 0).
    double davejac{0}, johnjac{0} /**< @copydoc davejac */;
    /// Interaction point, energy losses and Coulomb shift of the try.
    EventTarget target;
    /// Spectrometer target coordinates of both arms: true values (SP) and reconstructed values (RECON).
    DoubleArm SP, RECON /**< @copydoc SP */;
    /// Focal-plane tracks of both arms.
    DoubleArmFP FP;
    /// Recoil kinetic energy of the try [MeV].
    double Trec{0};
    NtupVars ntup;   ///< COMMON /ntuple_variables/ values for this event
    FermiMomentum pfermi;   ///< COMMON /pfermi/ of this event (task P9.1)
    /// COMMON /gnrl/ Mh, Mh2, Mh2_final when they change per try (doing_rho, task P13.2): used by the radiation
    /// (generate_rad), the montecarlo (spectrometer mass, recon energy, Eloss) instead of the set-up values when
    /// gnrl_per_event is true (set by the reaction).
    HadronMass gnrl;
    /// True when #gnrl holds the hadron mass of this try (set by the reaction; rho production only).
    bool gnrl_per_event{false};
};

/**
 * @brief The hadron arm that the reports and limits treat as "generated".
 *
 * That is vertex.p, except for rho production in the default (corrected) mode: there vertex.p holds the rho (with its
 * lab theta/phi in the slope fields) and the detected particle is the decay pion in orig.p.
 * @param rho_corrected True for rho production with `fortran_compat` off.
 * @param vertex The vertex event of the try.
 * @param orig The orig event of the try.
 * @return `orig.p` if @p rho_corrected, else `vertex.p` (a reference into the argument).
 * @note Origin: Fortran SIMC simc.f:259-261, 283-285, 316-317, event.f:limits_update l.22-24; FORTRAN-BUG D-057 F18
 *       (fortran_compat: vertex.p as the Fortran).
 */
inline const ArmFull& GeneratedHadronArm(bool rho_corrected, const Event& vertex, const Event& orig) {
    return rho_corrected ? orig.p : vertex.p;
}

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_EVENT_H
