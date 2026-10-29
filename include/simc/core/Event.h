// include/simc/core/Event.h
//
// Event data model. Mirrors structures.inc (Fortran derived types) field by field.
// Units follow the Fortran: MeV, MeV/c, cm, rad; delta in %; xptar/yptar are slopes
// (dx/dz, dy/dz); focal-plane dx/dy are slopes.
//
// Fortran COMMON /ntuple_variables/ (simulate.inc l.161-183: ntup, decdist) is carried
// per event in MainEvent::ntup instead of a global.
#ifndef SIMC_CORE_EVENT_H
#define SIMC_CORE_EVENT_H

#include <array>

namespace simc {

/// Fortran: type Cartesian
struct Cartesian {
    double x{0}, y{0}, z{0};
};

/// Fortran: type arm_full (delta %, xptar/yptar slopes, z cm, theta/phi rad, E/P MeV)
struct ArmFull {
    double delta{0}, xptar{0}, yptar{0}, z{0};
    double theta{0}, phi{0}, E{0}, P{0};
};

/// Fortran: type event (vertex, orig, recon)
struct Event {
    double Ein{0};
    double Em{0}, Pm{0};
    double Emiss{0}, Pmiss{0};
    double Pmx{0}, Pmy{0}, Pmz{0};
    double PmPar{0}, PmPer{0}, PmOop{0};
    double nu{0}, q{0}, Q2{0}, Trec{0}, W{0}, Mrec{0};
    double epsilon{0}, theta_pq{0}, theta_tarq{0}, phi_pq{0}, phi_targ{0};
    double beta{0}, phi_s{0}, phi_c{0};
    double zhad{0}, pt2{0}, xbj{0};
    ArmFull e, p;
    Cartesian ue, up, uq;
};

/// Fortran: type event_target
struct EventTarget {
    double x{0}, y{0}, z{0}, rasterx{0}, rastery{0};
    std::array<double, 3> teff{};    ///< Fortran teff(1:3) -> [0..2]
    std::array<double, 3> Eloss{};   ///< Fortran Eloss(1:3) -> [0..2]
    double Coulomb{0};
};

/// Fortran: type arm / arm2 (spectrometer quantities: delta %, yptar, xptar slopes, z cm)
struct ArmSP {
    double delta{0}, yptar{0}, xptar{0}, z{0};
};

/// Fortran: type double_arm
struct DoubleArm {
    ArmSP e, p;
};

/// Fortran: type arm_FP (x, y cm; dx, dy slopes; path cm)
struct ArmFP {
    double x{0}, dx{0}, y{0}, dy{0}, path{0};
};

struct DoubleArmFP {
    ArmFP e, p;
};

/// Fortran: type ntupvars + decdist (simulate.inc l.161-183), only the fields used by outputs.
struct NtupVars {
    double radphot{0};   ///< total radiated photon energy (MeV)
    double radarm{0};
    double resfac{0};
    double sigcm{0};
    double krel{0};
    double mm{0}, mmA{0}, t{0};
    double decdist{0};   ///< decay distance (cm), Fortran decdist
    double Mh2_final{0}; ///< hadron mass^2 after decay (MeV^2)
    double sigcm1{0};    ///< pion/kaon: first model cross section (physics_pion.f:126 sig_param_2021) [ub/MeV^2/rad]
    double sigcm2{0};    ///< pion/kaon: second model (physics_pion.f:148 MAID, W < 2 GeV only) [ub/MeV^2/rad]
    double dilu{0};      ///< SIDIS dilution (event.f:1525 ntup%dilu = 1.0; semi_dilution is dead code). Not an ntuple column
    double xfermi{0};    ///< SIDIS with Fermi motion: x_Bjorken = Q2/(2 P.q) (semi_physics.f:250, ntuple column 55)
    /// pi0 photons (pizero_decay.f, event.f:899-901): (E, px, py, pz) in the lab [MeV] (task P9.2)
    std::array<double, 4> gamma1{}, gamma2{};
    /// pi0 photon positions at the calorimeter [cm] (simc.f:1500-1560, mc_calo only; never written with HMS/SHMS)
    double xcal_gamma1{0}, ycal_gamma1{0}, xcal_gamma2{0}, ycal_gamma2{0};
    /// rho: the Breit-Wigner rho mass of the try [MeV] (generate_rho.f:79 ntup%rhomass = Mh; ntuple column Mrho) and
    /// the decay polar angle in the rho rest frame [rad] (rho_decay.f:66 ntup%rhotheta; column Thrho). Task P13.1.
    double rhomass{0}, rhotheta{0};
};

/// Fortran: COMMON /pfermi/ efer, pfer, pferx, pfery, pferz (simulate.inc l.212-217): the Fermi motion of the struck
/// nucleon of this event (event.f:generate l.330-371). pfer [MeV/c], (pferx, pfery, pferz) its unit vector (lab),
/// efer [MeV] its energy (targ%Mtar_struck without Fermi motion). Read by the pion/kaon cross sections
/// (physics_pion.f:peepi, jacobians.f:transform_to_cm) and the pion ntuple (results_write.f:141-149).
struct FermiMomentum {
    double pfer{0}, pferx{0}, pfery{0}, pferz{0}, efer{0};
};

/// Fortran: type event_main
struct MainEvent {
    double weight{0}, SF_weight{0}, gen_weight{0}, jacobian{0};
    double Ein_shift{0}, Ee_shift{0};
    double sigcc{0}, sigcc_recon{0}, sigcent{0};
    double epsilon{0}, theta_pq{0}, theta_tarq{0}, phi_pq{0}, phi_targ{0}, beta{0};
    double w{0}, t{0}, tmin{0}, q2{0};
    double pcm{0}, thetacm{0}, phicm{0}, wcm{0};
    double davejac{0}, johnjac{0};
    EventTarget target;
    DoubleArm SP, RECON;
    DoubleArmFP FP;
    double Trec{0};
    NtupVars ntup;   ///< COMMON /ntuple_variables/ values for this event
    FermiMomentum pfermi;   ///< COMMON /pfermi/ of this event (task P9.1)
};

}  // namespace simc

#endif  // SIMC_CORE_EVENT_H
