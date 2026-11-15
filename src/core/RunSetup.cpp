// src/core/RunSetup.cpp — see include/simc/core/RunSetup.h.
// The derived reaction flags of dbase.f (doing_hydpi, doing_deutkaon, doing_eep, ...: which reaction family on which
// target class), the reaction title of the .hist report and the minimal RunSetup made from the raw parameters. The
// full set-up (limits, theory tables, radiative corrections) is built by setup/SetupBuilder.cpp on top of this.
#include "simc/core/RunSetup.h"

#include <cmath>

namespace simc {

namespace {
constexpr double kPi = 3.141592653589793;
constexpr double kDegrad = 180.0 / kPi;   // constants.inc degrad
int nint(double x) { return static_cast<int>(std::lround(x)); }
// Collects the dbase.f STOP conditions (the caller decides whether a problem is fatal: the set-up refuses the run).
void problem(std::vector<std::string>* v, const std::string& s) { if (v) v->push_back(s); }
}  // namespace

// Fortran: dbase.f l.124-234
ReactionFlags DeriveReactionFlags(const SimcParameters& p, std::vector<std::string>* problems) {
    ReactionFlags f;
    f.doing_pion = p.doing_pion;
    f.doing_kaon = p.doing_kaon;
    f.doing_delta = p.doing_delta;
    f.doing_semi = p.doing_semi;
    f.doing_rho = p.doing_rho;
    f.doing_phsp = p.doing_phsp;
    f.doing_pizero = p.doing_pizero;
    f.doing_hplus = p.doing_hplus;
    f.do_fermi = p.do_fermi;
    f.using_tgt_field = p.using_tgt_field;
    const int A = nint(p.targ.A);   // the target class is decided on nint(targ%A): 1 = H, 2 = D, >= 3 nucleus

    // --- 1. semi-inclusive: doing_pion/doing_kaon with doing_semi select the SIDIS hadron (dbase.f:127-134)
    if (f.doing_pion && f.doing_semi) { f.doing_semipi = true; f.doing_semika = false; f.doing_pion = false; }
    if (f.doing_kaon && f.doing_semi) { f.doing_semika = true; f.doing_semipi = false; f.doing_kaon = false; }

    // --- 2. one reaction family per run (the first true flag of the dbase.f if-chain wins), split by target class;
    //        the coherent channels (which_* >= 10) are run as on a single nucleon (the doing_hyd* flag)
    if (f.doing_pion) {
        if (A == 1 && p.which_pion == 1 && !f.doing_pizero) problem(problems, "STOP Pi- production from Hydrogen not allowed!");
        if (A <= 2 && p.which_pion >= 10) problem(problems, "STOP Coherent production from Hydrogen/Deuterium not allowed!");
        if (A == 3 && p.which_pion == 11) problem(problems, "STOP Coherent Pi- production from 3He not allowed!");  // dbase.f:142-143
        if (A == 4 && p.which_pion >= 10) problem(problems, "STOP Coherent production from 4He not allowed!");     // dbase.f:144-145
        f.doing_hydpi = A == 1;
        f.doing_deutpi = A == 2;
        f.doing_hepi = A >= 3;
        if (p.which_pion >= 10) { f.doing_hydpi = true; f.doing_deutpi = false; f.doing_hepi = false; }
    } else if (f.doing_kaon) {
        if (A == 1 && p.which_kaon == 2) problem(problems, "STOP Sigma- production from Hydrogen not allowed!");
        if (A <= 2 && p.which_kaon >= 10) problem(problems, "STOP Coherent production from Hydrogen/Deuterium not allowed!");
        f.doing_hydkaon = A == 1;
        f.doing_deutkaon = A == 2;
        f.doing_hekaon = A >= 3;
        if (p.which_kaon >= 10) { f.doing_hydkaon = true; f.doing_deutkaon = false; f.doing_hekaon = false; }
    } else if (f.doing_delta) {
        f.doing_hyddelta = A == 1;
        f.doing_deutdelta = A == 2;
        f.doing_hedelta = A >= 3;
    } else if (f.doing_semi) {
        f.doing_hydsemi = A == 1;
        f.doing_deutsemi = A == 2;
        f.doing_hesemi = A >= 3;
        if (f.doing_hydsemi && f.do_fermi) f.do_fermi = false;
    } else if (f.doing_rho) {
        f.doing_hydrho = A == 1;
        f.doing_deutrho = A == 2;
        f.doing_herho = A == 3;   // FORTRAN-DIFF-CANDIDATE(D-018): dbase.f:213 uses .eq.3 (other he* flags .ge.3). Intended? (rho A >= 4 is refused by the set-up)
    } else {   // no production flag: (e,e'p)
        f.doing_eep = true;
        f.doing_hyd_elast = A == 1;
        f.doing_deuterium = A == 2;
        f.doing_heavy = A >= 3;
    }
    // --- 3. dbase.f:227-234: phsp switches the reaction flags off but keeps doing_hyd_elast & co. (the elastic
    //        kinematics of complete_ev still run for H + phsp, see FORTRAN_DIFFERENCES.md D-031)
    if (f.doing_phsp) {
        f.doing_eep = f.doing_pion = f.doing_kaon = f.doing_delta = f.doing_rho = false;
    }
    return f;
}

// Fortran: simc.f:report l.664-790 (the title line of the .hist; also used in the 'reaction not available' message)
std::string ReactionTitle(const ReactionFlags& f, const SimcParameters& p) {
    const bool plus = f.doing_hplus;
    if (f.doing_eep) {
        if (f.doing_hyd_elast) return "H(e,e'p)";
        if (f.doing_deuterium) return "D(e,e'p)";
        if (f.doing_heavy) return "A(e,e'p)";
    } else if (f.doing_semi) {
        const char* h = f.doing_semipi ? (plus ? "pi+" : "pi-") : (plus ? "k+" : "k-");
        const std::string t = p.targ.A == 1 ? "H" : (p.targ.A == 2 ? "D" : "A");
        return t + "(e,e'" + h + ")X";
    } else if (f.doing_rho) {
        return "H(e,e'rho)";
    } else if (f.doing_delta) {
        return f.doing_hyddelta ? "H(e,e'p)pi" : (f.doing_deutdelta ? "D(e,e'p)pi" : "A(e,e'p)pi");
    } else if (f.doing_pion) {
        return f.doing_deutpi ? "D(e,e'pi)" : ((f.doing_hydpi && p.targ.A == 1) ? "H(e,e'pi)" : "A(e,e'pi)");
    } else if (f.doing_kaon) {
        return f.doing_deutkaon ? "D(e,e'K)" : ((f.doing_hydkaon && p.targ.A == 1) ? "H(e,e'K)" : "A(e,e'K)");
    } else if (f.doing_phsp) {
        return "PHASE SPACE";
    }
    return "unknown";
}

RunSetup RunSetup::FromParameters(const SimcParameters& p, const std::string& run_name) {
    RunSetup s;
    s.params = p;
    s.flags = DeriveReactionFlags(p);
    s.run_name = run_name;
    // dbase.f:239-244: central angles in rad (the input is in degrees; the sign of theta only selects the side of the
    // beam, which is carried by the arm phi)
    s.e_theta = std::abs(p.spec.e.theta) / kDegrad;
    s.e_cos_th = std::cos(s.e_theta);
    s.e_sin_th = std::sin(s.e_theta);
    s.p_theta = std::abs(p.spec.p.theta) / kDegrad;
    s.p_cos_th = std::cos(s.p_theta);
    s.p_sin_th = std::sin(s.p_theta);
    return s;
}

}  // namespace simc
