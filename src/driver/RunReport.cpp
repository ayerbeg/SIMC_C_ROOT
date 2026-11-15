// src/driver/RunReport.cpp — see include/simc/driver/RunReport.h.
// Fortran: simc.f l.360-432 (normalisation), l.462-532 (.geni), report (l.644-1142).
// The .hist writer is in RunReportHist.cpp (records: RunReportRecords.cpp), the .geni writer in
// RunReportGeni.cpp; this file keeps the normalisation, FortranCtime and the summary sections.
#include "simc/driver/RunReport.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "simc/core/Errors.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SpectrometerId.h"
#include "simc/driver/FortranWriter.h"
#include "simc/core/Units.h"
#include "RunReportInternal.h"

namespace simc {
namespace driver {

using fortran::List;
using fortran::Record;
using constants::degrad;
using namespace report_detail;


// Parallel event loop: the counters and sums of a worker added to the run's (in worker order; RunReport.h).
void MergeTotals(RunTotals& a, const RunTotals& b) {
    a.ntried = a.ntried + b.ntried;
    a.nevent = a.nevent + b.nevent;
    a.ncontribute = a.ncontribute + b.ncontribute;
    a.ncontribute_no_rad_proton = a.ncontribute_no_rad_proton + b.ncontribute_no_rad_proton;
    a.npasscuts = a.npasscuts + b.npasscuts;
    a.wtcontribute = a.wtcontribute + b.wtcontribute;
    a.sum_sigcc = a.sum_sigcc + b.sum_sigcc;
    for (auto [x, y] : {std::pair<TwoArmSums*, const TwoArmSums*>{&a.sumerr, &b.sumerr}, {&a.sumerr2, &b.sumerr2}})
        for (auto [u, v] : {std::pair<TwoArmSums::Arm*, const TwoArmSums::Arm*>{&x->e, &y->e}, {&x->p, &y->p}}) {
            u->delta = u->delta + v->delta;
            u->xptar = u->xptar + v->xptar;
            u->yptar = u->yptar + v->yptar;
            u->ytar = u->ytar + v->ytar;
        }
    reactions::mergeContrib(a.contrib, b.contrib);
    reactions::mergeSlops(a.slop, b.slop);
}

// Fortran: simc.f l.360-432. The weight normalisation: normfac = luminosity * genvol * nevent / ntried, so that
// sum(weight) * normfac is a counting rate; plus the mean and r.m.s. of the reconstruction errors (resolutions).
Normalisation Normalise(const RunSetup& s, const RunTotals& t, double genvol) {
    Normalisation n;
    n.targetfac = s.targetfac;
    n.luminosity = s.luminosity;
    n.charge = s.params.EXPER.charge;
    // l.366 normfac = luminosity/ntried*nevent
    n.normfac = n.luminosity / static_cast<double>(t.ntried) * static_cast<double>(t.nevent);
    // l.374-394 (IReaction::genVolume)
    n.genvol = genvol;
    // l.396-398
    n.normfac = n.normfac * n.genvol;
    if (s.flags.doing_phsp) n.normfac = 1.0;
    n.wtcontribute = t.wtcontribute * n.normfac;

    // l.405-432 resolutions
    if (t.npasscuts > 1) {
        const double tmpnum = static_cast<double>(t.npasscuts);
        auto one = [&](const TwoArmSums::Arm& s1, const TwoArmSums::Arm& s2, TwoArmSums::Arm& ave, TwoArmSums::Arm& res) {
            ave.delta = s1.delta / tmpnum;
            ave.xptar = s1.xptar / tmpnum;
            ave.yptar = s1.yptar / tmpnum;
            ave.ytar = s1.ytar / tmpnum;
            auto r = [&](double a, double b) {
                const double m = a / tmpnum;
                return std::sqrt(fmaxF(0., (b / tmpnum) - m * m));
            };
            res.delta = r(s1.delta, s2.delta);
            res.xptar = r(s1.xptar, s2.xptar);
            res.yptar = r(s1.yptar, s2.yptar);
            res.ytar = r(s1.ytar, s2.ytar);
        };
        one(t.sumerr.e, t.sumerr2.e, n.aveerr.e, n.resol.e);
        one(t.sumerr.p, t.sumerr2.p, n.aveerr.p, n.resol.p);
    }
    return n;
}

std::string FortranCtime(long epoch_seconds) {
    std::time_t tt = static_cast<std::time_t>(epoch_seconds);
    char buf[64];
    std::tm tmv;
    localtime_r(&tt, &tmv);
    std::strftime(buf, sizeof buf, "%a %b %e %H:%M:%S %Y", &tmv);   // = ctime() without the newline
    return buf;
}

// ------------------------------------------------------------------------------------------ summary sections
// The RunSummary handed to the sinks at the end of the run: the counters and the normalisation (the numbers of the
// .hist), plus `extra` sections (central event, target, radiation, recon summary, limits, slops, kinematics) that the
// summary JSON writes as they are, so that simc_compare finds every .hist block in the C++ summary.
RunSummary MakeRunSummary(const RunSetup& s, const RunReport& r) {
    const RunTotals& t = r.totals;
    const Normalisation& n = r.norm;
    RunSummary sum;
    sum.run_name = r.run_name;
    sum.ngen = r.ngen;
    sum.ntried = t.ntried;
    sum.ncontribute = t.ncontribute;
    sum.ncontribute_no_rad_proton = t.ncontribute_no_rad_proton;
    sum.nevent = t.nevent;
    sum.wtcontribute = n.wtcontribute;
    sum.sum_sigcc = t.sum_sigcc;
    sum.central_sigcc = r.central.sigcc;
    sum.charge = n.charge;
    sum.targetfac = n.targetfac;
    sum.luminosity = n.luminosity;
    sum.genvol = n.genvol;
    sum.normfac = n.normfac;
    sum.random_seed = r.random_seed;
    sum.stops = r.stops;

    Json& x = sum.extra;
    const CentralEvent& c = r.central;
    x["central"] = {{"Q2_GeV2", c.Q2 / 1.e6}, {"q", c.q}, {"nu", c.nu}, {"Em", c.Em}, {"Pm", c.Pm}, {"W", c.W}, {"MM", c.MM},
                    {"e", {{"delta", c.e.delta}, {"xptar_mr", c.e.xptar}, {"yptar_mr", c.e.yptar}}},
                    {"p", {{"delta", c.p.delta}, {"xptar_mr", c.p.xptar}, {"yptar_mr", c.p.yptar}}}};
    const TargetInfo& tg = s.targ;
    auto tf = [](const TargetFunction& v) { return Json::array({v.ave, v.min, v.max}); };
    x["target"] = {{"A", tg.A}, {"Z", tg.Z}, {"mass_amu", tg.mass_amu}, {"mass_MeV", tg.M}, {"Mrec_amu", tg.mrec_amu},
                   {"Mrec_MeV", tg.Mrec}, {"Mtar_struck", tg.Mtar_struck}, {"Mrec_struck", tg.Mrec_struck},
                   {"rho", tg.rho}, {"thick_g_cm2", tg.thick}, {"angle_deg", tg.angle * degrad},
                   {"abundancy", tg.abundancy}, {"X0", tg.X0}, {"X0_cm", tg.X0_cm}, {"length", tg.length},
                   {"zoffset", tg.zoffset}, {"xoffset", tg.xoffset}, {"yoffset", tg.yoffset},
                   {"musc_nsig_max", tg.musc_nsig_max}, {"musc_max_beam_mr", tg.musc_max[0] * units::kMradPerRad},
                   {"musc_max_e_mr", tg.musc_max[1] * units::kMradPerRad}, {"musc_max_p_mr", tg.musc_max[2] * units::kMradPerRad},
                   {"losses", {{"Coulomb", tf(tg.Coulomb)}, {"Eloss_beam", tf(tg.Eloss[0])}, {"Eloss_e", tf(tg.Eloss[1])},
                               {"Eloss_p", tf(tg.Eloss[2])}, {"teff_beam", tf(tg.teff[0])}, {"teff_e", tf(tg.teff[1])},
                               {"teff_p", tf(tg.teff[2])}}}};
    if (s.sw.using_rad) {
        const CentralRad& cr = c.rad;
        x["radiative"] = {{"rad_flag", s.sw.rad_flag}, {"extrad_flag", s.sw.extrad_flag}, {"one_tail", s.params.one_tail},
                          {"intcor_mode", s.params.intcor_mode},
                          {"dE_edge_test", s.params.dE_edge_test},
                          {"Egamma_max", s.rad.Egamma_tot_max},
                          {"central", {{"hardcorfac", cr.hardcorfac}, {"etatzai", cr.etatzai},
                                       {"frac", {cr.frac[0], cr.frac[1], cr.frac[2]}},
                                       {"lambda", {cr.lambda[0], cr.lambda[1], cr.lambda[2]}},
                                       {"bt", {cr.bt[0], cr.bt[1]}},
                                       {"c_int", {cr.c_int[0], cr.c_int[1], cr.c_int[2], cr.c_int[3]}},
                                       {"c_ext", {cr.c_ext[0], cr.c_ext[1], cr.c_ext[2], cr.c_ext[3]}},
                                       {"c", {cr.c[0], cr.c[1], cr.c[2], cr.c[3]}},
                                       {"g_int", cr.g_int}, {"g_ext", cr.g_ext},
                                       {"g", {cr.g[0], cr.g[1], cr.g[2], cr.g[3]}}}}};
    }
    auto rec = [&](double a, double b, const char* u) { return Json{{"ave_error", a}, {"resolution", b}, {"units", u}}; };
    const TwoArmSums& ae = n.aveerr;
    const TwoArmSums& rs = n.resol;
    x["recon_summary"] = {
        {"e", {{"delta", rec(units::kPermillePerPercent * ae.e.delta, units::kPermillePerPercent * rs.e.delta, "x10^-3")}, {"xptar", rec(units::kMradPerRad * ae.e.xptar, units::kMradPerRad * rs.e.xptar, "mr")},
               {"yptar", rec(units::kMradPerRad * ae.e.yptar, units::kMradPerRad * rs.e.yptar, "mr")}, {"ytar", rec(units::kMmPerCm * ae.e.ytar, units::kMmPerCm * rs.e.ytar, "mm")}}},
        {"p", {{"delta", rec(units::kPermillePerPercent * ae.p.delta, units::kPermillePerPercent * rs.p.delta, "x10-3")}, {"xptar", rec(units::kMradPerRad * ae.p.xptar, units::kMradPerRad * rs.p.xptar, "mr")},
               {"yptar", rec(units::kMradPerRad * ae.p.yptar, units::kMradPerRad * rs.p.yptar, "mr")}, {"ytar", rec(units::kMmPerCm * ae.p.ytar, units::kMmPerCm * rs.p.ytar, "mm")}}}};

    // limits tables (labels as hist2json writes them: white space collapsed)
    const reactions::Contrib& cb = t.contrib;
    auto row = [](const char* label, double umin, double umax, double lo, double hi, const char* u) {
        return Json{{"label", label}, {"used", {umin, umax}}, {"found", {lo, hi}}, {"units", u}};
    };
    const GenLimits& g = s.gen;
    const EdgeTrue& V = s.VERTEXedge;
    const EdgeTrue& E = s.edge;
    const auto& SE = s.SPedge;
    x["limits"]["vertex"] = Json::array({
        row("E arm delta", g.e.delta.min, g.e.delta.max, cb.gen.e.delta.lo, cb.gen.e.delta.hi, "%"),
        row("E arm yptar", g.e.yptar.min * 1000., g.e.yptar.max * 1000., cb.gen.e.yptar.lo * 1000., cb.gen.e.yptar.hi * 1000., "mr"),
        row("E arm xptar", g.e.xptar.min * 1000., g.e.xptar.max * 1000., cb.gen.e.xptar.lo * 1000., cb.gen.e.xptar.hi * 1000., "mr"),
        row("P arm delta", g.p.delta.min, g.p.delta.max, cb.gen.p.delta.lo, cb.gen.p.delta.hi, "%"),
        row("P arm yptar", g.p.yptar.min * 1000., g.p.yptar.max * 1000., cb.gen.p.yptar.lo * 1000., cb.gen.p.yptar.hi * 1000., "mr"),
        row("P arm xptar", g.p.xptar.min * 1000., g.p.xptar.max * 1000., cb.gen.p.xptar.lo * 1000., cb.gen.p.xptar.hi * 1000., "mr"),
        row("sumEgen", g.sumEgen.min, g.sumEgen.max, cb.gen.sumEgen.lo, cb.gen.sumEgen.hi, "MeV"),
        row("Trec", V.Trec.min, V.Trec.max, cb.vertex.Trec.lo, cb.vertex.Trec.hi, "MeV"),
        row("Em", V.Em.min, V.Em.max, cb.vertex.Em.lo, cb.vertex.Em.hi, "MeV"),
        row("Pm", V.Pm.min, V.Pm.max, cb.vertex.Pm.lo, cb.vertex.Pm.hi, "MeV/c")});
    x["limits"]["original"] = Json::array({
        row("E arm E", E.e.E.min, E.e.E.max, cb.tru.e.E.lo, cb.tru.e.E.hi, "MeV"),
        row("E arm yptar", E.e.yptar.min * 1000., E.e.yptar.max * 1000., cb.tru.e.yptar.lo * 1000., cb.tru.e.yptar.hi * 1000., "mr"),
        row("E arm xptar", E.e.xptar.min * 1000., E.e.xptar.max * 1000., cb.tru.e.xptar.lo * 1000., cb.tru.e.xptar.hi * 1000., "mr"),
        row("P arm E", E.p.E.min, E.p.E.max, cb.tru.p.E.lo, cb.tru.p.E.hi, "MeV"),
        row("P arm yptar", E.p.yptar.min * 1000., E.p.yptar.max * 1000., cb.tru.p.yptar.lo * 1000., cb.tru.p.yptar.hi * 1000., "mr"),
        row("P arm xptar", E.p.xptar.min * 1000., E.p.xptar.max * 1000., cb.tru.p.xptar.lo * 1000., cb.tru.p.xptar.hi * 1000., "mr"),
        row("Em", fmaxF(-999999.999e0, E.Em.min), fminF(999999.999e0, E.Em.max), cb.tru.Em.lo, cb.tru.Em.hi, "MeV"),
        row("Pm", 0., 0., cb.tru.Pm.lo, cb.tru.Pm.hi, "MeV"),
        row("Trec", 0., 0., cb.tru.Trec.lo, cb.tru.Trec.hi, "MeV")});
    x["limits"]["spectrometer"] = Json::array({
        row("E arm delta", SE.e.delta.min, SE.e.delta.max, cb.SP.e.delta.lo, cb.SP.e.delta.hi, "%"),
        row("E arm yptar", SE.e.yptar.min * 1000., SE.e.yptar.max * 1000., cb.SP.e.yptar.lo * 1000., cb.SP.e.yptar.hi * 1000., "mr"),
        row("E arm xptar", SE.e.xptar.min * 1000., SE.e.xptar.max * 1000., cb.SP.e.xptar.lo * 1000., cb.SP.e.xptar.hi * 1000., "mr"),
        row("P arm delta", SE.p.delta.min, SE.p.delta.max, cb.SP.p.delta.lo, cb.SP.p.delta.hi, "%"),
        row("P arm yptar", SE.p.yptar.min * 1000., SE.p.yptar.max * 1000., cb.SP.p.yptar.lo * 1000., cb.SP.p.yptar.hi * 1000., "mr"),
        row("P arm xptar", SE.p.xptar.min * 1000., SE.p.xptar.max * 1000., cb.SP.p.xptar.lo * 1000., cb.SP.p.xptar.hi * 1000., "mr")});
    if (s.sw.using_rad) {
        x["limits"]["radiation"] = Json::array({
            row("Egamma(1)", 0., s.rad.Egamma1_max, cb.rad.Egamma[0].lo, cb.rad.Egamma[0].hi, "MeV"),
            row("Egamma(2)", 0., s.rad.Egamma2_max, cb.rad.Egamma[1].lo, cb.rad.Egamma[1].hi, "MeV"),
            row("Egamma(3)", 0., s.rad.Egamma3_max, cb.rad.Egamma[2].lo, cb.rad.Egamma[2].hi, "MeV"),
            row("Egamma_total", 0., s.rad.Egamma_tot_max, cb.rad.Egamma_total.lo, cb.rad.Egamma_total.hi, "MeV")});
    }
    const Slops& sl = t.slop;
    auto srow = [](const char* label, double used, double lo, double hi, const char* u) {
        return Json{{"label", label}, {"used", used}, {"min", lo}, {"max", hi}, {"units", u}};
    };
    x["slop"]["MC"] = Json::array({
        srow("E arm delta", sl.MC.e.delta.used, sl.MC.e.delta.lo, sl.MC.e.delta.hi, "%"),
        srow("E arm yptar", sl.MC.e.yptar.used * 1000., sl.MC.e.yptar.lo * 1000., sl.MC.e.yptar.hi * 1000., "mr"),
        srow("E arm xptar", sl.MC.e.xptar.used * 1000., sl.MC.e.xptar.lo * 1000., sl.MC.e.xptar.hi * 1000., "mr"),
        srow("P arm delta", sl.MC.p.delta.used, sl.MC.p.delta.lo, sl.MC.p.delta.hi, "%"),
        srow("P arm yptar", sl.MC.p.yptar.used * 1000., sl.MC.p.yptar.lo * 1000., sl.MC.p.yptar.hi * 1000., "mr"),
        srow("P arm xptar", sl.MC.p.xptar.used * 1000., sl.MC.p.xptar.lo * 1000., sl.MC.p.xptar.hi * 1000., "mr")});
    x["slop"]["total"] = Json::array({srow("Em", sl.total.Em.used, sl.total.Em.lo, sl.total.Em.hi, "MeV"),
                                      srow("Pm", 0., sl.total.Pm.lo, sl.total.Pm.hi, "MeV/c")});
    // kinematics in the internal (converted) units, like the .hist prints them
    auto arm = [](const SpecArm& a) {
        return Json{{"angle_deg", a.theta * degrad}, {"momentum", a.P}, {"offset_x", a.offset.x}, {"offset_y", a.offset.y},
                    {"offset_z", a.offset.z}, {"offset_xptar_mr", a.offset.xptar}, {"offset_yptar_mr", a.offset.yptar}};
    };
    x["kinematics"] = {{"Ebeam", s.Ebeam}, {"dE_over_E", s.dEbeam / s.Ebeam}, {"x_width", s.gen.xwid},
                       {"y_width", s.gen.ywid}, {"fr_pattern", s.targ.fr_pattern}, {"fr1", s.targ.fr1},
                       {"fr2", s.targ.fr2}, {"e", arm(s.spec.e)}, {"p", arm(s.spec.p)}};
    // (the extra sections replace BuildSummaryJson's sections of the same name: give them complete)
    x["counters"] = {{"ngen", r.ngen},
                     {"ntried", t.ntried},
                     {"ncontribute", t.ncontribute},
                     {"nco_no_rad_prot", t.ncontribute_no_rad_proton},
                     {"pct_no_rad_proton", (100. * static_cast<double>(t.ncontribute_no_rad_proton)) /
                                               fmaxF(static_cast<double>(t.ncontribute), 0.1e0)}};
    x["run"] = {{"tag", r.run_name},
                {"producer", "simc++"},
                {"random_seed", r.random_seed},
                {"rng", r.rng_type},                  // the random source (mt64 / ranlux / ...)
                {"fortran_compat", s.fortran_compat},  // simc_compare reads the pass_cuts mode from it
                {"reaction", ReactionTitle(s.flags, s.params)},
                {"begin_time", r.begin_time},
                {"end_time", r.end_time}};
    if (s.flags.doing_heavy) x["run"]["theory_file"] = s.theory.file;
    if (r.threads > 1) x["run"]["threads"] = r.threads;   // (absent for the sequential loop)
    return sum;
}

}  // namespace driver
}  // namespace simc
