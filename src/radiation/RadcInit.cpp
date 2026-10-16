// src/radiation/RadcInit.cpp — init.f:radc_init, radc_init_ev, basicrad_init_ev and RadcRun::FromSetup.
#include <cmath>
#include <ostream>

#include "simc/core/SimcConstants.h"
#include "simc/radiation/kernel/Brem.h"
#include "simc/radiation/kernel/Radc.h"

namespace simc {
namespace radiation {

RadcRun RadcRun::FromSetup(const RunSetup& s) {
    RadcRun r;
    const SimcParameters& p = s.params;
    r.rad_flag = s.sw.rad_flag;
    r.extrad_flag = s.sw.extrad_flag;   // after the radc_init default (SetupBuilder, init.f:621-626)
    r.intcor_mode = p.intcor_mode;
    r.use_expon = p.use_expon;
    r.one_tail = p.one_tail;
    r.use_offshell_rad = p.use_offshell_rad;
    r.hardwired_rad = s.sw.hardwired_rad;
    for (int i = 1; i <= 3; ++i) r.doing_tail[i] = s.sw.doing_tail[i - 1];
    r.etatzai = s.rad.etatzai;
    r.Egamma_tot_max = s.rad.Egamma_tot_max;
    r.Egamma1_max = s.rad.Egamma1_max;
    r.Egamma2_max = s.rad.Egamma2_max;
    r.Egamma3_max = s.rad.Egamma3_max;
    r.brem.produce_output = s.sw.produce_output;
    r.brem.exponentiate = s.sw.exponentiate;
    r.brem.include_hard = s.sw.include_hard;
    r.brem.calculate_spence = s.sw.calculate_spence;
    r.Egamma_gen_max = p.Egamma_gen_max;
    r.dE_edge_test = p.dE_edge_test;
    r.doing_eep = s.flags.doing_eep;
    r.doing_hyd_elast = s.flags.doing_hyd_elast;
    r.doing_deuterium = s.flags.doing_deuterium;
    r.doing_heavy = s.flags.doing_heavy;
    r.doing_pion = s.flags.doing_pion;
    r.doing_kaon = s.flags.doing_kaon;
    r.doing_rho = s.flags.doing_rho;
    r.doing_delta = s.flags.doing_delta;
    r.delta_tail1_case4 = s.flags.doing_delta && !s.fortran_compat;   // D-048 F2
    r.doing_semi = s.flags.doing_semi;
    r.targ_A = s.targ.A;
    r.Mh = s.Mh;
    r.Mh2 = s.Mh2;
    r.spec_e_P = s.spec.e.P;
    r.spec_p_P = s.spec.p.P;
    r.edge_e_E = s.edge.e.E;
    r.edge_p_E = s.edge.p.E;
    r.edge_Em = s.edge.Em;
    r.VERTEXedge_Em = s.VERTEXedge.Em;
    r.VERTEXedge_Pm = s.VERTEXedge.Pm;
    r.VERTEXedge_Trec = s.VERTEXedge.Trec;
    r.gen_sumEgen = s.gen.sumEgen;
    return r;
}

// Fortran: init.f:radc_init (l.576-651)
RadcInitResult radcInit(double targ_Z, double targ_L1, double targ_L2, int rad_flag, int extrad_flag, int use_expon,
                        bool debug1) {
    RadcInitResult r;
    r.extrad_flag = extrad_flag;
    if (extrad_flag == 0) {   // l.621-629
        if (rad_flag == 0) {
            r.extrad_flag = 3;
        } else if (rad_flag == 1 || rad_flag == 2 || rad_flag == 3) {
            r.extrad_flag = 1;
        }
    } else if (extrad_flag < 0) {
        throw FortranStop("Imbecile! check your stupid setting of EXTRAD_FLAG");
    }
    r.etatzai = (12.0 + (targ_Z + 1.) / (targ_Z * targ_L1 + targ_L2)) / 9.0;   // l.634
    r.brem.produce_output = debug1;                                            // l.638-647
    r.brem.exponentiate = use_expon == 1;
    r.brem.include_hard = true;
    r.brem.calculate_spence = true;
    return r;
}

// Fortran: init.f:radc_init_ev (l.655-728)
void radcInitEv(const RadcRun& run, RadcState& st, const MainEvent& main, const Event& vertex, std::ostream* log) {
    constexpr double zero = 0.0;
    double r;

    st.etta = 1.0;   // l.682
    for (int i = 1; i <= 2; ++i) st.bt[i] = run.etatzai * main.target.teff[i - 1];   // l.686-688

    // l.692-696. FORTRAN-DIFF-CANDIDATE(D-019, Q7): doing_tail(3) is passed as doing_proton (init.f:693), so switching
    // off tail 3 also removes the ep interference term from lambda(1:2). Question: intended?
    // FORTRAN-DIFF-CANDIDATE(D-005): this overwrites the .inp/JSON input lambda(1:3) on every event (init.f:692-695), so
    // that input has no effect. Question: remove the input, or use it when set?
    for (int i = 1; i <= 3; ++i)
        st.lambda[i] = lambdaDave(i, 1, run.doing_tail[3], vertex.Ein, vertex.e.E, vertex.p.E, vertex.p.P,
                                  vertex.e.theta, st.lambda_warned, log);
    st.rad_proton_this_ev = st.lambda[3] > 0;

    // l.700-715
    const double Ecutoff = 450.;
    if (run.intcor_mode == 0) {
        // FORTRAN-DIFF-CANDIDATE(D-019, Q3): schwinger does not set dsoft_prime; the SAVEd value of an earlier
        // brem/bremos call (0 in a SIMC run with intcor_mode = 0) is used for g(4) below (init.f:701-715), so g(4)
        // loses its internal part (0.019443 vs 0.098637 at the P1.2 kinematics). Question: is intcor_mode = 0 meant to
        // have a derivative term?
        r = schwinger(Ecutoff, vertex, true, run.targ_A, st.etta, run.use_expon, st.initev_dsoft, st.initev_dhard);
    } else {
        if (!run.use_offshell_rad) {
            r = brem(run.brem, vertex.Ein, vertex.e.E, Ecutoff, st.rad_proton_this_ev, st.initev_dsoft,
                     st.initev_dhard, st.initev_dsoft_prime, log);
        } else {
            r = bremos(run.brem, Ecutoff, zero, zero, vertex.Ein, vertex.e.P * vertex.ue.x, vertex.e.P * vertex.ue.y,
                       vertex.e.P * vertex.ue.z, zero, zero, zero, vertex.p.P * vertex.up.x, vertex.p.P * vertex.up.y,
                       vertex.p.P * vertex.up.z, vertex.p.E, st.rad_proton_this_ev, st.initev_dsoft, st.initev_dhard,
                       st.initev_dsoft_prime, log);
        }
    }
    (void)r;
    st.hardcorfac = 1. / (1. - st.initev_dhard);
    st.g[4] = -st.initev_dsoft_prime * Ecutoff + st.bt[1] + st.bt[2];

    basicradInitEv(st, vertex.Ein, vertex.e.E, vertex.p.E);   // l.719

    for (int i = 1; i <= 3; ++i) st.frac[i] = st.g[i] / st.g[0];   // l.723-725
}

// Fortran: init.f:basicrad_init_ev (l.732-813). gamma = intrinsic = libm tgamma (D-019).
void basicradInitEv(RadcState& st, double e1, double e2, double e3) {
    using constants::euler;
    using constants::Mp;
    constexpr double one = 1.;
    const double e[4] = {0, e1, e2, e3};
    auto& g = st.g;
    auto& c = st.c;
    auto& c_int = st.c_int;
    auto& c_ext = st.c_ext;
    const auto& lambda = st.lambda;
    const auto& bt = st.bt;

    // l.758-761
    g[1] = lambda[1] + bt[1];
    g[2] = lambda[2] + bt[2];
    g[3] = lambda[3];
    g[0] = g[1] + g[2] + g[3];

    // Internal constants (l.765-781)
    c_int[1] = lambda[1] / std::pow(e[1] * e[2], lambda[1] / 2.);
    c_int[2] = lambda[2] / std::pow(e[1] * e[2], lambda[2] / 2.);
    c_int[3] = lambda[3] / std::pow(Mp * e[3], lambda[3] / 2.);
    for (int i = 1; i <= 3; ++i) c_int[i] = c_int[i] * std::exp(-euler * lambda[i]) / std::tgamma(one + lambda[i]);
    st.g_int = lambda[1] + lambda[2] + lambda[3];
    c_int[0] = c_int[1] * c_int[2] * st.g_int / lambda[1] / lambda[2];
    if (lambda[3] > 0) c_int[0] = c_int[0] * c_int[3] / lambda[3];
    c_int[0] = c_int[0] * std::tgamma(one + lambda[1]) * std::tgamma(one + lambda[2]) * std::tgamma(one + lambda[3]) /
               std::tgamma(one + st.g_int);

    // External constants (l.785-791)
    for (int i = 1; i <= 2; ++i) c_ext[i] = bt[i] / std::pow(e[i], bt[i]) / std::tgamma(one + bt[i]);
    c_ext[3] = 0.0;
    st.g_ext = bt[1] + bt[2];
    c_ext[0] = c_ext[1] * c_ext[2] * st.g_ext / bt[1] / bt[2];
    c_ext[0] = c_ext[0] * std::tgamma(one + bt[1]) * std::tgamma(one + bt[2]) / std::tgamma(one + st.g_ext);

    // Internal + external constants (l.795-799)
    for (int i = 1; i <= 2; ++i)
        c[i] = c_int[i] * c_ext[i] * g[i] / lambda[i] / bt[i] * std::tgamma(one + lambda[i]) *
               std::tgamma(one + bt[i]) / std::tgamma(one + g[i]);
    c[3] = c_int[3];

    // Combined tails (l.803-810)
    c[0] = c[1] * c[2] * g[0] / g[1] / g[2];
    if (g[3] > 0) c[0] = c[0] * c[3] / g[3];
    c[0] = c[0] * std::tgamma(one + g[1]) * std::tgamma(one + g[2]) * std::tgamma(one + g[3]) / std::tgamma(one + g[0]);
    c[4] = g[4] / std::pow(e1 * e2, g[4]) / std::tgamma(one + g[4]);
    if (g[3] > 0) c[4] = c[4] / std::pow(e3, g[4]);
}

}  // namespace radiation
}  // namespace simc
