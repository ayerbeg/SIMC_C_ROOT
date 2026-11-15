// src/radiation/Radc.cpp — see include/simc/radiation/kernel/Radc.h.
// The event side of the SIMC radiative corrections (radc.f, peaking approximation): generate_rad picks one of the
// three radiating tails (incoming electron, scattered electron, outgoing hadron), draws the photon energies from the
// basicrad power law, rebuilds the event kinematics with complete_ev and computes the radiative weight
// (peaked_rad_weight). Called through the IRadiator plug-in SimcRadiator by the reactions' generate step.
// Fortran: simc_gfortran/radc.f. Operation order follows the Fortran source; comments give the line numbers.
#include "simc/radiation/kernel/Radc.h"

#include <algorithm>
#include <cmath>
#include <ostream>

#include "simc/core/SimcConstants.h"
#include "simc/radiation/kernel/Brem.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace radiation {

namespace {
inline double sq(double x) { return x * x; }
// Fortran MIN/MAX of two reals as gfortran evaluates them for ordinary (non-NaN) values.
inline double fmin2(double a, double b) { return b < a ? b : a; }
inline double fmax2(double a, double b) { return b > a ? b : a; }
}  // namespace

// Fortran: radc.f:basicrad (l.3-88). One photon energy from the distribution g*E^(g-1) between emin and emax
// (inverse-transform sampling), with its weight relative to the full distribution.
void basicrad(const RadcState& st, int& itail, double Egamma_lo, double Egamma_hi, double& Egamma, double& weight,
              double& val_reciprocal, IRandom& rng) {
    // Initialize (l.47-49)
    Egamma = 0.0;
    weight = 0.0;
    val_reciprocal = 0.0;

    // ... is radiation in this direction turned off? (l.53-57)
    // FORTRAN-QUIRK(D-019, Q6): itail 0 (peaked basis) -> 4; it is reset to 0 only on the normal exit.
    if (itail == 0) itail = 4;
    if (st.g[itail] <= 0) {
        weight = 1.0;
        return;
    }

    // ... do we have a chance? (l.61)
    if (Egamma_hi <= Egamma_lo || Egamma_hi <= 0) return;

    // l.65-74
    const double power_hi = std::pow(Egamma_hi, st.g[itail]);
    double power_lo = 0.0;
    if (Egamma_lo > 0) power_lo = std::pow(Egamma_lo, st.g[itail]);
    const double ymin = power_lo / power_hi;
    const double y = ymin + rng.uniform() * (1. - ymin);
    const double x = std::pow(y, 1. / st.g[itail]);
    Egamma = x * Egamma_hi;

    // l.78-84
    if (Egamma > 0) val_reciprocal = std::pow(Egamma, 1. - st.g[itail]) * (power_hi - power_lo) / st.g[itail];
    weight = st.c[itail] / st.g[itail] * (power_hi - power_lo);

    if (itail == 4) itail = 0;
}

// Fortran: radc.f:gamma (l.92-116). DEAD CODE in SIMC (D-019): kept only for the L1 table.
// FORTRAN-QUIRK: the recursion factor is (y+1+i) instead of (y+i), wrong for x >= 2.
double gammaRadc(double x) {
    double gamma = 1.0;
    const int n = static_cast<int>(std::lround((x - 1) - 0.5));   // nint
    const double y = x - 1 - n;
    if (n != 0) {
        const int s = n >= 0 ? 1 : -1;   // sign(1,n)
        for (int i = s; s > 0 ? i <= n : i >= n; i += s) {
            const double f = y + 1 + i;
            gamma = gamma * (s == 1 ? f : 1. / f);   // (y+1+i)**s, s = +-1 (__powidf2)
        }
    }
    // y**k with a constant exponent: gfortran's powi_table chains (y3 = y2*y, y4 = y2*y2, y5 = y3*y2)
    const double y2 = y * y, y3 = y2 * y, y4 = y2 * y2, y5 = y3 * y2;
    gamma = gamma * (1. - 0.5748646 * y + 0.9512363 * y2 - 0.6998588 * y3 + 0.4245549 * y4 - 0.1010678 * y5);
    return gamma;
}

// Fortran: radc.f:peaked_rad_weight (l.523-646). The weight of a radiated photon of energy Egamma in the peaking
// approximation: the ratio of the full (internal + external) spectrum to the basicrad spectrum it was drawn from.
double peakedRadWeight(const RadcRun& run, RadcState& st, const Event& vertex, double Egamma, double emin,
                       double emax, double basicrad_val_reciprocal, double basicrad_weight, std::ostream* log) {
    (void)basicrad_val_reciprocal;   // l.545: only touched to avoid an unused-variable warning
    constexpr double zero = 0.0;
    const double ein = vertex.Ein;
    const double eout = vertex.e.E;
    const double eul = 0.577215665;
    double r, dhard = 0;
    double dsoft_int_primemin = 0, dsoft_int_primemax = 0;

    // ... External (l.566-585)
    // FORTRAN-DIFF-CANDIDATE(D-019, Q4): phi_ext (rad_flag 0), dsoft_ext*, and the Friedrich terms are computed and
    // never used in the returned weight, so for rad_flag = 0 the weight does not depend on Egamma (radc.f:566-642).
    // Question: is the rad_flag = 0 peaked weight meant to include the external (phi_ext) factor?
    double phi_ext = 1.0;
    if (run.extrad_flag <= 2) {
        phi_ext = extradPhi(run, st, 0, ein, eout, Egamma);
        if (run.rad_flag == 1) return basicrad_weight * phi_ext;
        double dsoft_extmin = 0;
        if (emin > 0) dsoft_extmin = std::log(st.g_ext / st.c_ext[0] / std::pow(emin, st.g_ext));
        const double dsoft_extmax = std::log(st.g_ext / st.c_ext[0] / std::pow(emax, st.g_ext));
        const double dsoft_ext_prime = -st.g_ext / Egamma;
        (void)dsoft_extmin;
        (void)dsoft_extmax;
        (void)dsoft_ext_prime;
    } else {
        const double t1 = st.bt[1] / run.etatzai;
        const double t2 = st.bt[2] / run.etatzai;
        double dsoft_ext1, dsoft_ext1_prime, dsoft_ext2, dsoft_ext2_prime;
        extradFriedrich(run.etatzai, ein, Egamma, t1, dsoft_ext1, dsoft_ext1_prime);
        extradFriedrich(run.etatzai, eout, Egamma, t2, dsoft_ext2, dsoft_ext2_prime);
    }

    // ... Internal (l.590-630)
    if (run.rad_flag == 0) {
        if (!run.use_offshell_rad) {
            if (emin > 0) {
                r = brem(run.brem, ein, eout, emin, st.rad_proton_this_ev, st.prw_dsoft_intmin, dhard,
                         dsoft_int_primemin, log);
            } else {
                st.prw_dsoft_intmin = 1.0;
            }
            r = brem(run.brem, ein, eout, emax, st.rad_proton_this_ev, st.prw_dsoft_intmax, dhard,
                     dsoft_int_primemax, log);
        } else {
            if (emin > 0) {
                r = bremos(run.brem, emin, zero, zero, ein, vertex.e.E * vertex.ue.x, vertex.e.E * vertex.ue.y,
                           vertex.e.E * vertex.ue.z, zero, zero, zero, vertex.p.P * vertex.up.x,
                           vertex.p.P * vertex.up.y, vertex.p.P * vertex.up.z, vertex.p.E, st.rad_proton_this_ev,
                           st.prw_dsoft_intmin, dhard, dsoft_int_primemin, log);
            } else {
                st.prw_dsoft_intmin = 1.0;
            }
            r = bremos(run.brem, emax, zero, zero, ein, vertex.e.E * vertex.ue.x, vertex.e.E * vertex.ue.y,
                       vertex.e.E * vertex.ue.z, zero, zero, zero, vertex.p.P * vertex.up.x, vertex.p.P * vertex.up.y,
                       vertex.p.P * vertex.up.z, vertex.p.E, st.rad_proton_this_ev, st.prw_dsoft_intmax, dhard,
                       dsoft_int_primemax, log);
        }
        (void)r;
    } else {
        // FORTRAN-DIFF-CANDIDATE(D-019, Q4): dsoft_int / dsoft_int_prime are not used below; for rad_flag = 1 with
        // extrad_flag = 3 the weight uses the SAVEd dsoft_intmin/max of an earlier call (0 in a fresh process).
        // Question: should this branch compute dsoft_intmin/max (radc.f:590-620)?
        const double dsoft_int = std::log(st.g_int / st.c_int[0] / std::pow(Egamma, st.g_int));
        const double dsoft_int_prime = -st.g_int / Egamma;
        (void)dsoft_int;
        (void)dsoft_int_prime;
    }

    // ... All together now (l.634-643). gamma = intrinsic = libm tgamma (D-019), through the memo st.gamma (same bits).
    double w;
    if (emin > 0) {
        w = st.c_ext[0] / st.g_ext *
            (std::exp(-st.prw_dsoft_intmax) * std::pow(emax, st.g_ext) -
             std::exp(-st.prw_dsoft_intmin) * std::pow(emin, st.g_ext));
    } else {
        w = st.c_ext[0] / st.g_ext * (std::exp(-st.prw_dsoft_intmax) * std::pow(emax, st.g_ext));
    }
    const double g4 = st.g[4];
    w = w * std::exp(-eul * g4) / st.gamma(1. + g4) * st.gamma(1. + g4 - st.bt[1] - st.bt[2]) *
        st.gamma(1. + st.bt[1]) * st.gamma(1. + st.bt[2]) / st.gamma(1. + g4);
    if (w < 0) w = 0;
    return w;
}

// Fortran: radc.f:extrad_friedrich (l.650-664). External bremsstrahlung correction (Friedrich) for radiator thickness trad.
void extradFriedrich(double etatzai, double Ei, double Ecutoff, double trad, double& dbrem, double& dbrem_prime) {
    const double x = Ecutoff / Ei;
    dbrem = trad * (-(etatzai - 0.5) - etatzai * std::log(x) + etatzai * x - 0.5 * sq(x));
    dbrem_prime = -trad / Ei * (etatzai / x - etatzai + x);
}

// Fortran: radc.f:extrad_phi (l.668-707). The external-bremsstrahlung shape factor phi(Egamma/E) of tail itail.
double extradPhi(const RadcRun& run, const RadcState& st, int itail, double E1, double E2, double Egamma) {
    const double E[3] = {0, E1, E2};
    double phi = 1.0;   // CASE 1
    if (run.extrad_flag == 2) {
        if (itail == 0) {
            phi = 1. - (st.bt[1] / E[1] + st.bt[2] / E[2]) / (st.g[1] + st.g[2]) * Egamma;
        } else if (itail == 1 || itail == 2) {
            phi = 1. - st.bt[itail] / E[itail] / st.g[itail] * Egamma;
        }
    } else if (run.extrad_flag == 3) {
        if (itail == 0)
            throw FortranStop(
                "Idiot! a multiplicative factor EXTRAD_PHI is not defined for peaking approx and EXTRAD_FLAG>2!");
        if (itail == 1 || itail == 2) {
            const double x = Egamma / E[itail];
            const double t = st.bt[itail] / run.etatzai;
            phi = phi * (1. - x + sq(x) / run.etatzai) *
                  std::exp(t * ((run.etatzai - 0.5) - run.etatzai * x + sq(x) / 2.)) * st.gamma(1. + st.bt[itail]);
        }
    }
    return phi;
}

// Fortran: radc.f:schwinger (l.711-739). Schwinger's soft (dsoft) and hard (dhard) internal correction.
// FORTRAN-QUIRK: exp(dsoft) when use_expon == 0 (opposite of brem's convention); only dhard is used in SIMC.
double schwinger(double Ecutoff, const Event& vertex, bool include_hard, double targ_A, double etta, int use_expon,
                 double& dsoft, double& dhard) {
    using constants::alpi;
    using constants::amu;
    const double lq = std::log(vertex.Q2 / constants::Me2) - 1.0;
    const double s2 = sq(std::sin(vertex.e.theta / 2.));
    const double b = 1. + 2. * vertex.nu * s2 / (targ_A * amu);
    const double spence_ = spen(1. - s2) - 2.5893784;
    dsoft = alpi * lq * std::log(vertex.Ein / (sq(etta)) * vertex.e.E * b / (sq(Ecutoff)));
    dhard = -alpi * (2.166666 * lq + spence_ - sq(std::log(vertex.Ein / vertex.e.E)) / 2.0);
    double s;
    if (use_expon == 0)
        s = std::exp(dsoft);
    else
        s = 1. + dsoft;
    if (include_hard) s = s / (1. - dhard);
    return s;
}

// Fortran: radc.f:spen (l.743-764)
double spen(double x) {
    double y = 1.0, s = 0.0;
    int i = 0;
    while (i <= 100 && std::fabs(y) > std::fabs(s) * 1.e-4) {
        i = i + 1;
        y = x * y;
        s = s + y / static_cast<double>(i * i);
    }
    return s;
}

// Fortran: radc.f:lambda_dave (l.768-821). The effective radiator thickness lambda of each tail (internal + external).
double lambdaDave(int itail, int plus_flag, bool doing_proton, double e1, double e2, double e3, double p3, double th,
                  bool& warned, std::ostream* log) {
    using constants::alpi;
    using constants::Me;
    double lambda_dave = 0.0;
    double plus_term = 0.0;
    if (plus_flag == 1 && itail < 3) {
        plus_term = std::log((1. - std::cos(th)) / 2.);
        if (doing_proton) plus_term = plus_term + 2. * std::log(e1 / e2);
    }
    if (itail == 1) {
        lambda_dave = alpi * (2. * std::log(2. * e1 / Me) - 1. + plus_term);
    } else if (itail == 2) {
        lambda_dave = alpi * (2. * std::log(2. * e2 / Me) - 1. + plus_term);
    } else if (itail == 3) {
        if (doing_proton) {
            lambda_dave = alpi * ((e3 / p3) * std::log((e3 + p3) / (e3 - p3)) - 2.);
            if (lambda_dave < 0) {
                lambda_dave = 0.0;
                if (!warned) {
                    if (log) {
                        *log << " +-----------------------------------+\n"
                             << " | WARNING: lambda_dave(3)<0, think  |\n"
                             << " | about running with proton rad OFF |\n"
                             << " +-----------------------------------+\n";
                    }
                    warned = true;
                }
            }
        } else {
            lambda_dave = 0.0;
        }
    }
    return lambda_dave;
}

// Fortran: radc.f:generate_rad (l.120-519)
bool generateRad(const RadcRun& run, RadcState& st, MainEvent& main, Event& vertex, Event& orig, IRandom& rng,
                 const CompleteEvFn& complete_ev, std::ostream* log) {
    using constants::Mp;
    int peaked_basis_flag;
    double basicrad_weight = 0, basicrad_val_reciprocal = 0;

    // Initialize (l.179-194)
    for (int i = 1; i <= 3; ++i) {
        st.lim[i] = 0;
        st.phot[i] = 0;
    }
    double rad_weight = 1;
    for (int i = 1; i <= 3; ++i) st.Egamma_used[i] = 0.0;
    main.ntup.radphot = 0.;
    main.ntup.radarm = 0.;
    const bool force2 = false;   // the force2 block is commented out in the Fortran (l.223-226)

    // Which tail has to be radiated? (l.198-216)
    // FORTRAN-DIFF-CANDIDATE(D-019): with one_tail != 0 ntail is still drawn from frac(1:3), which includes the
    // disabled tails; a disabled tail is then skipped (doing_tail tests at l.234/358), the event passes with no photon
    // and gen_weight = 1/hardcorfac. Question: should frac be restricted to the enabled tails?
    peaked_basis_flag = 1;
    if (run.rad_flag <= 1) {
        peaked_basis_flag = 0;
        const double x = rng.uniform();
        if (x >= st.frac[1] + st.frac[2]) {
            st.ntail = 3;
        } else if (x >= st.frac[1]) {
            st.ntail = 2;
        } else {
            st.ntail = 1;
        }
    } else if (run.rad_flag == 2) {
        st.ntail = static_cast<int>(rng.uniform() * 3.) + 1;
        if (st.ntail == 4) st.ntail = 3;
    } else if (run.rad_flag == 3) {
        st.ntail = 0;
    } else {
        throw FortranStop("Idiot! rad_flag is set stupidly");
    }
    const int ntail = st.ntail;

    // l.229-230
    const double max_delta_Trec =
        fmax2((vertex.Trec - run.VERTEXedge_Trec.min), (run.VERTEXedge_Trec.max - vertex.Trec));

    // RADIATE TAIL #1: the incident electron tail (l.234-337)
    if (run.doing_tail[1] && (ntail == 0 || ntail == 1)) {
        if (run.doing_heavy) {   // CASE 1 (l.249-257)
            st.Egamma_min[1] = vertex.Em - run.VERTEXedge_Em.max - max_delta_Trec;
            st.Egamma_max[1] = vertex.Em - run.VERTEXedge_Em.min + max_delta_Trec;
            if (ntail != 0) st.Egamma_max[1] = fmin2(st.Egamma_max[1], vertex.Em - run.edge_Em.min + max_delta_Trec);
        } else if (run.doing_hyd_elast) {   // CASE 2 (l.269-276)
            double ebeam_max = Mp * run.edge_e_E.max / (Mp - run.edge_e_E.max * (1. - vertex.ue.z));
            if (ebeam_max < 0) ebeam_max = 1.e10;
            const double ebeam_min = Mp * run.edge_e_E.min / (Mp - run.edge_e_E.min * (1. - vertex.ue.z));
            st.Egamma_min[1] = vertex.Ein - ebeam_max;
            st.Egamma_max[1] = vertex.Ein - ebeam_min;
            st.Egamma_max[1] = fmin2(st.Egamma_max[1], run.edge_Em.max);
        } else if (run.doing_deuterium) {   // CASE 3 (l.280-282)
            st.Egamma_max[1] = fmin2(run.Egamma1_max, run.gen_sumEgen.max - vertex.e.E);
            // FORTRAN-DIFF-CANDIDATE(D-019): radc.f:283-286, for ntail = 0 Egamma_min(1) keeps its stale COMMON value.
            // Question: should it be set for ntail = 0 too (D(e,e'p) only)?
            if (ntail != 0) st.Egamma_min[1] = run.gen_sumEgen.min - vertex.e.E;
        } else if (run.doing_pion || run.doing_kaon || run.doing_rho || run.doing_semi ||   // CASE 4 (l.288-293)
                   run.delta_tail1_case4) {
            // FORTRAN-BUG(D-048 F2): original Fortran radc.f:288-289:
            // `else if (doing_pion .or. doing_kaon .or. doing_rho .or. doing_semi) then` has no doing_delta, so for
            // doing_delta Egamma_min(1)/Egamma_max(1) are never assigned (0), basicrad returns weight 0 (radc.f:61) and
            // EVERY try with ntail = 1 fails (about 44 % of the tries at the Delta validation kinematics); corrected because
            // doing_delta generates the electron energy in the pion sumEgen window (init.f:455-463, the "generated
            // ELECTRON energy limits" branch, the same as pion production), so the pion limit applies. fortran_compat
            // keeps the Fortran.
            st.Egamma_min[1] = 0.;
            st.Egamma_max[1] = run.gen_sumEgen.max - vertex.e.E;
        }

        // l.297-307
        st.Egamma_max[1] = fmin2(st.Egamma_max[1], run.Egamma1_max);
        st.Egamma_min[1] = st.Egamma_min[1] - run.dE_edge_test;
        st.Egamma_max[1] = st.Egamma_max[1] + run.dE_edge_test;
        st.lim[1] = st.Egamma_max[1];
        if (run.hardwired_rad) {
            st.Egamma_max[1] = run.Egamma_gen_max;
            if (st.lim[1] > run.Egamma_gen_max && log)
                *log << " SIMC found an Egamma(1) limit above your hardwired value\n";
        }

        // ... radiate (l.312-317)
        int it = 1 * peaked_basis_flag;
        basicrad(st, it, st.Egamma_min[1], st.Egamma_max[1], st.Egamma_used[1], basicrad_weight,
                 basicrad_val_reciprocal, rng);
        if (basicrad_weight <= 0) return false;

        // ... adjust kinematics and call complete_ev (l.322-328)
        vertex.Ein = vertex.Ein - st.Egamma_used[1];
        const bool evsuccess = complete_ev(main, vertex, rng);
        if (!evsuccess) return false;
        // FORTRAN-DIFF-CANDIDATE(D-019, Q5): complete_ev -> radc_init_ev has just recomputed hardcorfac, frac, g, c for
        // the reduced Ein; they are used from here on (peaked_rad_weight, tails 2/3, the final division radc.f:512).
        // Question: is the post-radiation hardcorfac the intended normalisation?
        rad_weight = rad_weight * basicrad_weight;
        if (run.rad_flag <= 1) {
            // FORTRAN-DIFF-CANDIDATE(D-019, Q5): radc.f:330-333 overwrites rad_weight (rad_weight*basicrad_weight is
            // discarded) instead of multiplying. Question: intended (peaked_rad_weight already contains basicrad_weight)?
            rad_weight = peakedRadWeight(run, st, vertex, st.Egamma_used[1], st.Egamma_min[1], st.Egamma_max[1],
                                         basicrad_val_reciprocal, basicrad_weight, log);
        } else {
            rad_weight = rad_weight * extradPhi(run, st, 1, vertex.Ein, vertex.e.E, st.Egamma_used[1]);
        }
    }

    // VERTEXedges (SF limits) (l.342-354)
    if (run.doing_heavy) {
        if (vertex.Em < run.VERTEXedge_Em.min || vertex.Em > run.VERTEXedge_Em.max || vertex.Pm < run.VERTEXedge_Pm.min ||
            vertex.Pm > run.VERTEXedge_Pm.max)
            return false;
    }

    // RADIATE TAIL #2: the scattered electron tail (l.358-410)
    if (run.doing_tail[2] && (ntail == 0 || ntail == 2)) {
        st.Egamma_min[2] = vertex.e.E - run.edge_e_E.max;
        st.Egamma_max[2] = vertex.e.E - run.edge_e_E.min;
        if (run.doing_eep) {
            st.Egamma_max[2] =
                fmin2(st.Egamma_max[2], (run.edge_Em.max - vertex.Em) - st.Egamma_used[1] + max_delta_Trec);
            if (ntail != 0 || !st.rad_proton_this_ev)
                st.Egamma_min[2] =
                    fmax2(st.Egamma_min[2], (run.edge_Em.min - vertex.Em) - st.Egamma_used[1] - max_delta_Trec);
        }
        st.Egamma_max[2] = fmin2(st.Egamma_max[2], run.Egamma_tot_max - st.Egamma_used[1]);
        st.Egamma_min[2] = st.Egamma_min[2] - run.dE_edge_test;
        st.Egamma_max[2] = st.Egamma_max[2] + run.dE_edge_test;
        st.lim[2] = st.Egamma_max[2];
        if (run.hardwired_rad) {
            st.Egamma_max[2] = run.Egamma_gen_max;
            if (st.lim[2] > run.Egamma_gen_max && log)
                *log << " SIMC found an Egamma(2) limit above your hardwired value\n";
        }
        int it = 2 * peaked_basis_flag;
        basicrad(st, it, st.Egamma_min[2], st.Egamma_max[2], st.Egamma_used[2], basicrad_weight,
                 basicrad_val_reciprocal, rng);
        if (basicrad_weight <= 0) return false;
        rad_weight = rad_weight * basicrad_weight;
        if (run.rad_flag <= 1) {
            rad_weight = peakedRadWeight(run, st, vertex, st.Egamma_used[2], st.Egamma_min[2], st.Egamma_max[2],
                                         basicrad_val_reciprocal, basicrad_weight, log);
        } else {
            rad_weight = rad_weight * extradPhi(run, st, 2, vertex.Ein, vertex.e.E, st.Egamma_used[2]);
        }
    }

    // RADIATE TAIL #3: the scattered proton tail (l.414-461)
    if (st.rad_proton_this_ev && (ntail == 0 || ntail == 3)) {
        st.Egamma_min[3] = vertex.p.E - run.edge_p_E.max;
        st.Egamma_max[3] = vertex.p.E - run.edge_p_E.min;
        if (run.doing_eep) {
            st.Egamma_max[3] = fmin2(st.Egamma_max[3], (run.edge_Em.max - vertex.Em) - st.Egamma_used[1] -
                                                           st.Egamma_used[2] + max_delta_Trec);
            st.Egamma_min[3] = fmax2(st.Egamma_min[3], (run.edge_Em.min - vertex.Em) - st.Egamma_used[1] -
                                                           st.Egamma_used[2] - max_delta_Trec);
        }
        st.Egamma_max[3] = fmin2(st.Egamma_max[3], run.Egamma_tot_max - st.Egamma_used[1] - st.Egamma_used[2]);
        st.Egamma_min[3] = st.Egamma_min[3] - run.dE_edge_test;
        st.Egamma_max[3] = st.Egamma_max[3] + run.dE_edge_test;
        st.lim[3] = st.Egamma_max[3];
        if (run.hardwired_rad) {
            st.Egamma_max[3] = run.Egamma_gen_max;
            if (st.lim[3] > run.Egamma_gen_max && log)
                *log << " SIMC found an Egamma(3\t) limit above your hardwired value\n";   // TAB as in radc.f:444
        }
        int it = 3 * peaked_basis_flag;
        basicrad(st, it, st.Egamma_min[3], st.Egamma_max[3], st.Egamma_used[3], basicrad_weight,
                 basicrad_val_reciprocal, rng);
        if (basicrad_weight <= 0) return false;
        rad_weight = rad_weight * basicrad_weight;
        if (run.rad_flag <= 1) {
            rad_weight = peakedRadWeight(run, st, vertex, st.Egamma_used[3], st.Egamma_min[3], st.Egamma_max[3],
                                         basicrad_val_reciprocal, basicrad_weight, log);
        } else {
            rad_weight = rad_weight * extradPhi(run, st, 3, vertex.Ein, vertex.e.E, st.Egamma_used[3]);
        }
    }

    // ORIG = vertex + radiation (l.476-500)
    orig = vertex;
    orig.Ein = vertex.Ein + st.Egamma_used[1];
    orig.e.E = vertex.e.E - st.Egamma_used[2];
    if (orig.e.E <= 0e0) return false;
    orig.e.P = orig.e.E;
    orig.e.delta = (orig.e.P - run.spec_e_P) / run.spec_e_P * 100.;
    orig.p.E = vertex.p.E - st.Egamma_used[3];
    // COMMON /gnrl/ Mh, Mh2: per try for doing_rho (the Breit-Wigner rho mass of generate_rho)
    const double Mh = main.gnrl_per_event ? main.gnrl.Mh : run.Mh;
    const double Mh2 = main.gnrl_per_event ? main.gnrl.Mh2 : run.Mh2;
    if (orig.p.E <= Mh) return false;
    orig.p.P = std::sqrt(sq(orig.p.E) - Mh2);
    orig.p.delta = (orig.p.P - run.spec_p_P) / run.spec_p_P * 100.;

    // l.503-508
    st.phot[1] = st.Egamma_used[1];
    st.phot[2] = st.Egamma_used[2];
    st.phot[3] = st.Egamma_used[3];
    main.ntup.radphot = st.phot[1] + st.phot[2] + st.phot[3];
    main.ntup.radarm = ntail;
    if (force2) main.ntup.radarm = 12;

    // l.512-515
    main.gen_weight = main.gen_weight * rad_weight / st.hardcorfac;
    if (force2) main.gen_weight = main.gen_weight * st.frac[2];
    return true;
}

}  // namespace radiation
}  // namespace simc
