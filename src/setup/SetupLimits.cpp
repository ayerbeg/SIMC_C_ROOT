// src/setup/SetupLimits.cpp — BuildRunSetup phase 3 (task P15.3 split of SetupBuilder.cpp, no logic change):
// init.f:limits_init (l.91-575, called at dbase.f:685): slops, edges, extreme_trip_thru_target through the target
// model, VERTEXedge, radiation limits, generation limits, and the diagnostic histogram axes (init.f:520-535).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <fstream>
#include <sstream>

#include "simc/core/SimcConstants.h"
#include "simc/core/SpectrometerId.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/core/Units.h"
#include "SetupInternal.h"

namespace simc {
namespace setup_detail {

namespace K = constants;

bool BuildLimits(RunSetup& s, const SimcParameters& p, ITargetLossProvider* losses, double benharEmMax) {
    auto msg = [&](const std::string& m) { s.messages.push_back(m); };
    const ReactionFlags& f = s.flags;
    RunSwitches& sw = s.sw;
    TargetInfo& t = s.targ;
    GenLimits& gen = s.gen;
    // =========================================================================================
    // init.f:limits_init (l.91-575), called at dbase.f:685
    // =========================================================================================
    Slops& slop = s.slop;
    struct SlopParam { double d, t, p; };
    // simulate.inc:17-36
    // init.f:151-196 (the Fortran also lists the retired code 6 with the SHMS, init.f:176/199; rejected at input).
    auto slopFor = [](int arm_code, SlopParam& sp) {
        const std::optional<SpectrometerId> id = SpectrometerFromCode(arm_code);
        if (!id) return false;
        switch (*id) {
            case SpectrometerId::HMS: sp = {0.5, 0.005, 0.005}; return true;
            case SpectrometerId::SOS: sp = {1.0, 0.008, 0.008}; return true;
            case SpectrometerId::HRSR: sp = {0.5, 0.005, 0.005}; return true;
            case SpectrometerId::HRSL: sp = {0.5, 0.005, 0.005}; return true;
            case SpectrometerId::SHMS: sp = {0.5, 0.005, 0.005}; return true;
            default: return false;   // calorimeters: slop stays 0
        }
    };
    SlopParam sp;
    if (sw.using_E_arm_montecarlo && slopFor(p.electron_arm, sp)) {   // init.f:151-173
        slop.MC.e.delta.used = sp.d;
        slop.MC.e.yptar.used = sp.t;
        slop.MC.e.xptar.used = sp.p;
    }
    if (sw.using_P_arm_montecarlo && slopFor(p.hadron_arm, sp)) {     // init.f:174-196
        slop.MC.p.delta.used = sp.d;
        slop.MC.p.yptar.used = sp.t;
        slop.MC.p.xptar.used = sp.p;
    }

    // init.f:201-212 add slop to SPedges
    auto widen = [](ArmCuts& c, const Slops::Arm& a) {
        c.delta.min = c.delta.min - a.delta.used;
        c.delta.max = c.delta.max + a.delta.used;
        c.yptar.min = c.yptar.min - a.yptar.used;
        c.yptar.max = c.yptar.max + a.yptar.used;
        c.xptar.min = c.xptar.min - a.xptar.used;
        c.xptar.max = c.xptar.max + a.xptar.used;
    };
    widen(s.SPedge.e, slop.MC.e);
    widen(s.SPedge.p, slop.MC.p);

    // init.f:216-225 TRUE edges
    EdgeTrue& edge = s.edge;
    const double dE = p.dE_edge_test;
    const SpecArm& se = s.spec.e;
    const SpecArm& spp = s.spec.p;
    edge.e.E.min = (1. + s.SPedge.e.delta.min / 100.) * se.P + t.Coulomb.min - dE;
    edge.e.E.max = (1. + s.SPedge.e.delta.max / 100.) * se.P + t.Coulomb.max + dE;
    Cut pp;
    pp.min = (1. + s.SPedge.p.delta.min / 100.) * spp.P - dE;
    pp.max = (1. + s.SPedge.p.delta.max / 100.) * spp.P + dE;
    pp.min = std::max(0.001e0, pp.min);
    edge.p.E.min = std::sqrt(pp.min * pp.min + s.Mh2);
    edge.p.E.max = std::sqrt(pp.max * pp.max + s.Mh2);

    // init.f:228-246 extreme theta/z for extreme_trip_thru_target
    ExtremeTripArgs xa;
    xa.electron_arm = p.electron_arm;
    xa.hadron_arm = p.hadron_arm;
    xa.fortran_compat = s.fortran_compat;
    xa.Ebeam = s.Ebeam;
    {
        const ArmCuts& a = s.SPedge.e;
        xa.the.max = std::acos((se.cos_th - se.sin_th * a.yptar.max) /
                               std::sqrt(1. + a.yptar.max * a.yptar.max + a.xptar.max * a.xptar.max));
        xa.the.min = std::acos((se.cos_th - se.sin_th * a.yptar.min) / std::sqrt(1. + a.yptar.min * a.yptar.min));
    }
    {
        const ArmCuts& a = s.SPedge.p;
        xa.thp.max = std::acos((spp.cos_th - spp.sin_th * a.yptar.max) /
                               std::sqrt(1. + a.yptar.max * a.yptar.max + a.xptar.max * a.xptar.max));
        xa.thp.min = std::acos((spp.cos_th - spp.sin_th * a.yptar.min) / std::sqrt(1. + a.yptar.min * a.yptar.min));
    }
    xa.z.min = -0.5 * t.length;
    xa.z.max = 0.5 * t.length;
    xa.pe = edge.e.E;
    xa.pp = pp;
    xa.m = s.Mh;
    losses->extremeLosses(xa, t);   // init.f:245
    s.target_losses_from_provider = true;
    if (!sw.using_Eloss) {   // init.f:247-252
        for (auto& e : t.Eloss) {
            e.min = 0.0;
            e.max = 0.0;
        }
    }
    if (!sw.mc_smear) t.musc_max = {0., 0., 0.};   // init.f:254-258

    // init.f:260-271
    edge.e.E.min = edge.e.E.min + t.Eloss[1].min;
    edge.e.E.max = edge.e.E.max + t.Eloss[1].max;
    edge.p.E.min = edge.p.E.min + t.Eloss[2].min;
    edge.p.E.max = edge.p.E.max + t.Eloss[2].max;
    edge.e.yptar.min = s.SPedge.e.yptar.min - t.musc_max[1];
    edge.e.yptar.max = s.SPedge.e.yptar.max + t.musc_max[1];
    edge.e.xptar.min = s.SPedge.e.xptar.min - t.musc_max[1];
    edge.e.xptar.max = s.SPedge.e.xptar.max + t.musc_max[1];
    edge.p.yptar.min = s.SPedge.p.yptar.min - t.musc_max[2];
    edge.p.yptar.max = s.SPedge.p.yptar.max + t.musc_max[2];
    edge.p.xptar.min = s.SPedge.p.xptar.min - t.musc_max[2];
    edge.p.xptar.max = s.SPedge.p.xptar.max + t.musc_max[2];

    // init.f:283 measured beam energy used in recon
    s.Ebeam_vertex_ave = s.Ebeam + t.Coulomb.ave - t.Eloss[0].ave;

    // init.f:287-298 slops in energies and momenta
    const double Ebeam_max = s.Ebeam + s.dEbeam / 2. - t.Eloss[0].min + t.Coulomb.max;
    const double Ebeam_min = s.Ebeam - s.dEbeam / 2. - t.Eloss[0].max + t.Coulomb.min;
    const double slop_Coulomb = t.Coulomb.max - t.Coulomb.ave;
    double slop_Ebeam = s.dEbeam / 2. + slop_Coulomb;
    double slop_Ee = slop.MC.e.delta.used / 100. * se.P + slop_Coulomb;
    const double r = std::sqrt(edge.p.E.max * edge.p.E.max - s.Mh2);
    double slop_Ep = std::sqrt((r + slop.MC.p.delta.used / 100. * spp.P) * (r + slop.MC.p.delta.used / 100. * spp.P) +
                               s.Mh2) -
                     edge.p.E.max;
    slop_Ebeam = slop_Ebeam + (t.Eloss[0].max - t.Eloss[0].min);
    slop_Ee = slop_Ee + (t.Eloss[1].max - t.Eloss[1].min);
    slop_Ep = slop_Ep + (t.Eloss[2].max - t.Eloss[2].min);

    if (f.doing_heavy) {   // init.f:300-305 'reconstructed' Em cuts
        slop.total.Em.used = slop_Ebeam + slop_Ee + slop_Ep + dE;
        edge.Em.min = s.cuts.Em.min - slop.total.Em.used;
        edge.Em.max = s.cuts.Em.max + slop.total.Em.used;
        edge.Em.min = std::max(0.e0, edge.Em.min);
    }

    // init.f:325-371 VERTEXedge Em/Pm
    EdgeTrue& V = s.VERTEXedge;
    const TheoryData& th = s.theory;
    auto pdistMax = [&]() { return s.pdist.pval.empty() ? 0.0 : s.pdist.pval.back(); };   // pval(nump)
    if (f.doing_hyd_elast) {
        V.Em.min = 0.0;
        V.Em.max = 0.0;
        V.Pm.min = 0.0;
        V.Pm.max = 0.0;
    } else if (s.flags.doing_deuterium) {
        V.Em.min = K::Mp + K::Mn - t.M;
        V.Em.max = K::Mp + K::Mn - t.M;
        V.Pm.min = 0.0;
        V.Pm.max = th.Pm.empty() ? 0.0 : std::max(std::abs(th.Pm[0].min), std::abs(th.Pm[0].max));
    } else if (f.doing_heavy) {
        V.Pm.min = 0.0;
        V.Pm.max = 0.0;
        if (p.use_benhar_sf) {
            V.Pm.max = 790.0;
        } else {
            for (int i = 0; i < th.nrhoPm; ++i) {
                const double t1 = std::max(std::abs(th.Pm[i].min), std::abs(th.Pm[i].max));
                V.Pm.max = std::max(V.Pm.max, t1);
            }
        }
        V.Em.min = th.E_Fermi;
        V.Em.max = 1000.;   // Need Egamma_tot_max for good limit.
    } else if (f.doing_hydpi || f.doing_hydkaon || f.doing_hyddelta || f.doing_hydrho) {
        V.Em.min = 0.0;
        V.Em.max = 0.0;
        V.Pm.min = 0.0;
        V.Pm.max = 0.0;
    } else if (f.doing_deutpi || f.doing_deutkaon || f.doing_deutdelta || f.doing_deutrho) {
        V.Em.min = K::Mp + K::Mn - t.M;
        V.Em.max = K::Mp + K::Mn - t.M;
        V.Pm.min = 0.0;
        V.Pm.max = pdistMax();
    } else if (f.doing_hepi || f.doing_hekaon || f.doing_hedelta || f.doing_herho) {
        V.Em.min = t.Mtar_struck + t.Mrec - t.M;
        V.Em.max = benharEmMax;   // Emval(numEm) (sf_lookup_init, dbase.f:594-621)
        V.Pm.min = 0.0;
        V.Pm.max = pdistMax();
    } else if (f.doing_semi) {
        V.Em.min = 0.0;
        V.Em.max = 0.0;
        V.Pm.min = 0.0;
        V.Pm.max = 0.0;
    }
    msg("E_bind = " + std::to_string(V.Em.min) + " MeV in limits_init (QF only)");

    // init.f:379-392 recoiling (A-1) system
    if (f.doing_hyd_elast || f.doing_hydpi || f.doing_hydkaon || f.doing_hyddelta || f.doing_hydrho || f.doing_semi) {
        V.Mrec.min = 0.0;
        V.Mrec.max = 0.0;
        V.Trec.min = 0.0;
        V.Trec.max = 0.0;
    } else {
        V.Mrec.min = t.M - t.Mtar_struck + V.Em.min;
        V.Mrec.max = t.M - t.Mtar_struck + V.Em.max;
        V.Trec.min = std::sqrt(V.Mrec.max * V.Mrec.max + V.Pm.min * V.Pm.min) - V.Mrec.max;
        V.Trec.max = std::sqrt(V.Mrec.min * V.Mrec.min + V.Pm.max * V.Pm.max) - V.Mrec.min;
    }

    // init.f:399-406
    if (f.doing_eep || f.doing_semi) {
        V.Trec_struck.min = 0.;
        V.Trec_struck.max = 0.;
    } else {
        V.Trec_struck.min = 0.;
        V.Trec_struck.max = Ebeam_max + t.Mtar_struck - t.Mrec_struck - edge.e.E.min - edge.p.E.min - V.Em.min -
                            V.Trec.min;
    }

    // init.f:415-432 radiation limits
    RadSetup& rad = s.rad;
    rad.Egamma_tot_max = Ebeam_max + t.Mtar_struck - t.Mrec_struck - edge.e.E.min - edge.p.E.min - V.Em.min -
                         V.Trec.min - V.Trec_struck.min;
    if (f.doing_heavy) {
        const double t2 = (edge.Em.max - V.Em.min) + (V.Trec.max - V.Trec.min);
        rad.Egamma_tot_max = std::min(rad.Egamma_tot_max, t2);
    }
    if (sw.hardwired_rad) rad.Egamma_tot_max = p.Egamma_gen_max;
    if (!sw.using_rad) rad.Egamma_tot_max = 0.0;
    if (sw.doing_tail[0]) rad.Egamma1_max = rad.Egamma_tot_max;
    if (sw.doing_tail[1]) rad.Egamma2_max = rad.Egamma_tot_max;
    if (sw.doing_tail[2]) rad.Egamma3_max = rad.Egamma_tot_max;
    if (f.doing_heavy) {
        V.Em.min = std::max(V.Em.min, edge.Em.min - rad.Egamma_tot_max);
        V.Em.max = std::min(V.Em.max, edge.Em.max);
    }

    // init.f:441-475 generated summed energies
    if (f.doing_hyd_elast) {
        gen.sumEgen.min = 0.0;
        gen.sumEgen.max = 0.0;
    } else if (f.doing_heavy) {
        gen.sumEgen.max = Ebeam_max + t.Mtar_struck - V.Trec.min - V.Em.min;
        gen.sumEgen.min = Ebeam_min + t.Mtar_struck - V.Trec.max - V.Em.max - rad.Egamma1_max;
        gen.sumEgen.max = std::min(gen.sumEgen.max, edge.e.E.max + edge.p.E.max + rad.Egamma_tot_max);
        gen.sumEgen.min = std::max(gen.sumEgen.min, edge.e.E.min + edge.p.E.min);
    } else if (f.doing_semi) {
        gen.sumEgen.max = Ebeam_max + t.Mtar_struck - t.Mrec_struck;
        gen.sumEgen.min = edge.e.E.min + edge.p.E.min;
    } else {   // generated ELECTRON energy limits
        gen.sumEgen.max = Ebeam_max + t.Mtar_struck - t.Mrec_struck - edge.p.E.min - V.Em.min - V.Trec.min -
                          V.Trec_struck.min;
        gen.sumEgen.min = Ebeam_min + t.Mtar_struck - t.Mrec_struck - edge.p.E.max - V.Em.max - V.Trec.max -
                          V.Trec_struck.max - rad.Egamma_tot_max;
        gen.sumEgen.max = std::min(gen.sumEgen.max, edge.e.E.max + rad.Egamma2_max);
        gen.sumEgen.min = std::max(gen.sumEgen.min, edge.e.E.min);
    }
    gen.sumEgen.min = gen.sumEgen.min - dE;
    gen.sumEgen.max = gen.sumEgen.max + dE;
    gen.sumEgen.min = std::max(0.e0, gen.sumEgen.min);

    // init.f:480-498 E arm generation limits
    if (f.doing_hyd_elast) {
        gen.e.E.min = edge.e.E.min;
        gen.e.E.max = edge.e.E.max + rad.Egamma2_max;
    } else if (s.flags.doing_deuterium || f.doing_pion || f.doing_kaon || f.doing_rho || f.doing_delta) {
        gen.e.E.min = gen.sumEgen.min;
        gen.e.E.max = gen.sumEgen.max;
    } else if (f.doing_heavy || f.doing_semi) {
        gen.e.E.min = gen.sumEgen.min - edge.p.E.max - rad.Egamma3_max;
        gen.e.E.max = gen.sumEgen.max - edge.p.E.min;
    }
    gen.e.E.min = std::max(gen.e.E.min, edge.e.E.min);
    gen.e.E.max = std::min(gen.e.E.max, edge.e.E.max + rad.Egamma2_max);
    gen.e.delta.min = (gen.e.E.min / se.P - 1.) * 100.;
    gen.e.delta.max = (gen.e.E.max / se.P - 1.) * 100.;
    gen.e.yptar.min = edge.e.yptar.min;
    gen.e.yptar.max = edge.e.yptar.max;
    gen.e.xptar.min = edge.e.xptar.min;
    gen.e.xptar.max = edge.e.xptar.max;

    // init.f:503-517 P arm generation limits
    if (f.doing_hyd_elast || s.flags.doing_deuterium || f.doing_pion || f.doing_kaon || f.doing_rho || f.doing_delta) {
        gen.p.E.min = edge.p.E.min;
        gen.p.E.max = edge.p.E.max + rad.Egamma3_max;
    } else if (f.doing_heavy || f.doing_semi) {
        gen.p.E.min = gen.sumEgen.min - edge.e.E.max - rad.Egamma2_max;
        gen.p.E.max = gen.sumEgen.max - edge.e.E.min;
    }
    gen.p.E.min = std::max(gen.p.E.min, edge.p.E.min);
    gen.p.E.max = std::min(gen.p.E.max, edge.p.E.max + rad.Egamma3_max);
    gen.p.delta.min = (std::sqrt(gen.p.E.min * gen.p.E.min - s.Mh2) / spp.P - 1.) * 100.;
    gen.p.delta.max = (std::sqrt(gen.p.E.max * gen.p.E.max - s.Mh2) / spp.P - 1.) * 100.;
    gen.p.yptar.min = edge.p.yptar.min;
    gen.p.yptar.max = edge.p.yptar.max;
    gen.p.xptar.min = edge.p.xptar.min;
    gen.p.xptar.max = edge.p.xptar.max;

    // init.f:520-535 axis specs for the diagnostic histograms (H%gen = H%geni = H%RECON)
    constexpr double nHbins = 50;   // histograms.inc
    DiagHistAxes& H = s.hist;
    H.e.delta.min = gen.e.delta.min;
    H.e.yptar.min = gen.e.yptar.min;
    H.e.xptar.min = -gen.e.xptar.max;
    H.e.delta.bin = (gen.e.delta.max - gen.e.delta.min) / nHbins;
    H.e.yptar.bin = (gen.e.yptar.max - gen.e.yptar.min) / nHbins;
    H.e.xptar.bin = (gen.e.xptar.max - gen.e.xptar.min) / nHbins;
    H.p.delta.min = gen.p.delta.min;
    H.p.yptar.min = gen.p.yptar.min;
    H.p.xptar.min = -gen.p.xptar.max;
    H.p.delta.bin = (gen.p.delta.max - gen.p.delta.min) / nHbins;
    H.p.yptar.bin = (gen.p.yptar.max - gen.p.yptar.min) / nHbins;
    H.p.xptar.bin = (gen.p.xptar.max - gen.p.xptar.min) / nHbins;
    H.Em.min = V.Em.min;
    H.Em.bin = (std::max(100.e0, V.Em.max) - V.Em.min) / nHbins;
    H.Pm.min = V.Pm.min;
    H.Pm.bin = (std::max(100.e0, V.Pm.max) - V.Pm.min) / nHbins;
    return true;
}

}  // namespace setup_detail
}  // namespace simc
