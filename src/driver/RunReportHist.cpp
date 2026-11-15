// src/driver/RunReportHist.cpp — part of the report writer (RunReport.h): the .hist text, record by record as the
// Fortran WRITE statements produce it (input echo, reaction, normalisation, central values, resolutions)
// (Fortran: simc.f:report, l.644-1142).
#include "simc/driver/RunReport.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <sstream>
#include <stdexcept>

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


// Fortran: simc.f:report (l.644-1142)
std::string FormatHist(const RunSetup& s, const RunReport& r) {
    const SimcParameters& p = s.params;
    const ReactionFlags& f = s.flags;
    const RunTotals& t = r.totals;
    const Normalisation& n = r.norm;
    const reactions::Contrib& c = t.contrib;
    const CentralEvent& central = r.central;
    Out o;

    // Run time (l.663-664)
    o.line("");
    // character*30 timestring (left-justified, blank-padded) written with a30
    auto t30 = [](std::string v) { v.resize(30, ' '); return v; };
    { Record x; x.x(1).lit("BEGIN Time: ").a(t30(r.begin_time), 30); o.line(x); }
    { Record x; x.x(1).lit("END Time:   ").a(t30(r.end_time), 30); o.line(x); }

    // Kinematics (l.667-815)
    o.list("KINEMATICS:");
    reactionLines(o, s);
    o.line(kin1("Ebeam", s.Ebeam, "MeV"));
    o.line(kin1("(dE/E)beam", s.dEbeam / s.Ebeam, "(full wid)"));
    o.line(kin1("x-width", s.gen.xwid, "cm"));
    o.line(kin1("y-width", s.gen.ywid, "cm"));
    { Record x; x.x(9).a("fr_pattern", 12).lit(" = ").i(s.targ.fr_pattern, 15).x(2).a("1=square,2=circ", 16); o.line(x); }
    o.line(kin1("fr1", s.targ.fr1, "cm"));
    o.line(kin1("fr2", s.targ.fr2, "cm"));
    o.list(" ");
    { Record x; x.x(9).x(18).a("____E arm____", 15).x(2).a("____P arm____", 15).x(2); o.line(x); }
    o.line(kin2("angle", s.spec.e.theta * degrad, s.spec.p.theta * degrad, "deg"));
    o.line(kin2("momentum", s.spec.e.P, s.spec.p.P, "MeV/c"));
    o.line(kin2("x offset", s.spec.e.offset.x, s.spec.p.offset.x, "cm"));
    o.line(kin2("y offset", s.spec.e.offset.y, s.spec.p.offset.y, "cm"));
    o.line(kin2("z offset", s.spec.e.offset.z, s.spec.p.offset.z, "cm"));
    // FORTRAN-BUG(D-018): original Fortran simc.f:815-818: `write(iun,...) 'xptar offset', spec%e%offset%xptar,
    // spec%p%offset%xptar, 'mr'` (and yptar) -- the values are in rad (dbase.f:472-475 divides the mrad input by 1000);
    // corrected because the label says mr: the corrected .hist prints mrad (= the input values).
    const double offFac = s.fortran_compat ? 1.0 : units::kMradPerRad;
    o.line(kin2("xptar offset", s.spec.e.offset.xptar * offFac, s.spec.p.offset.xptar * offFac, "mr"));
    o.line(kin2("yptar offset", s.spec.e.offset.yptar * offFac, s.spec.p.offset.yptar * offFac, "mr"));

    o.list("                      VALUES FOR \"CENTRAL\" EVENT:");
    o.line(kin2("delta", central.e.delta, central.p.delta, "%"));
    o.line(kin2("xptar", central.e.xptar, central.p.xptar, "mr"));
    o.line(kin2("yptar", central.e.yptar, central.p.yptar, "mr"));
    o.line(cen1("Q2", central.Q2 / 1.e6, "(GeV/c)^2"));
    o.line(cen1("q", central.q, "MeV/c"));
    o.line(cen1("nu", central.nu, "MeV"));
    o.line(cen1("recon Em", central.Em, "MeV"));
    o.line(cen1("recon Pm", central.Pm, "MeV/c"));
    o.line(cen1("recon  W", central.W, "MeV/c"));
    o.line(cen1("recon MM", central.MM, "MeV/c"));

    // Target (l.837-875)
    const TargetInfo& tg = s.targ;
    o.list("TARGET specs:");
    o.line(targ2("A", tg.A, " ", "Z", tg.Z, " "));
    o.line(targ2("mass", tg.mass_amu, "amu", "mass", tg.M, "MeV"));
    o.line(targ2("Mrec", tg.mrec_amu, "amu", "Mrec", tg.Mrec, "MeV"));
    o.line(targ2("Mtar_struc", tg.Mtar_struck, "MeV", "Mrec_struc", tg.Mrec_struck, "MeV"));
    o.line(targ2("rho", tg.rho, "g/cm3", "thick", tg.thick, "g/cm2"));
    o.line(targ2("angle", tg.angle * degrad, "deg", "abundancy", tg.abundancy, "%"));
    o.line(targ2("X0", tg.X0, "g/cm2", "X0_cm", tg.X0_cm, "cm"));
    o.line(targ2("length", tg.length, "cm", "zoffset", tg.zoffset, "cm"));
    o.line(targ2("xoffset", tg.xoffset, "cm", "yoffset", tg.yoffset, "cm"));
    { Record x; x.t(12).a("__ave__", 15).a("__lo__", 15).a("__hi__", 15); o.line(x); }
    o.line(targ3("Coulomb", tg.Coulomb, "MeV"));
    o.line(targ3("Eloss_beam", tg.Eloss[0], "MeV"));
    o.line(targ3("Eloss_e", tg.Eloss[1], "MeV"));
    o.line(targ3("Eloss_p", tg.Eloss[2], "MeV"));
    o.line(targ3("teff_beam", tg.teff[0], "radlen"));
    o.line(targ3("teff_e", tg.teff[1], "radlen"));
    o.line(targ3("teff_p", tg.teff[2], "radlen"));
    o.line(targ1("musc_nsig_max", tg.musc_nsig_max, " "));
    o.line(targ1("musc_max_beam", tg.musc_max[0] * units::kMradPerRad, "mr"));
    o.line(targ1("musc_max_e", tg.musc_max[1] * units::kMradPerRad, "mr"));
    o.line(targ1("musc_max_p", tg.musc_max[2] * units::kMradPerRad, "mr"));

    // Flags (l.878-913). Values as the Fortran COMMON holds them after dbase.f.
    o.list("FLAGS:");
    o.line(flags({{"doing_eep", f.doing_eep}, {"doing_kaon", f.doing_kaon}, {"doing_pion", f.doing_pion}}));
    o.line(flags({{"doing_semi", f.doing_semi}, {"doing_rho", f.doing_rho}, {"doing_hplus", f.doing_hplus}}));
    o.line(flags({{"doing_semipi", f.doing_semipi}, {"doing_semika", f.doing_semika}, {"doing_pizero", f.doing_pizero}}));
    o.line(flags({{"doing_delta", f.doing_delta}, {"doing_phsp", f.doing_phsp}}));
    {
        Record x;
        x.x(5).x(2).a("which_pion", 19).lit("=").i(p.which_pion, 2);
        x.x(2).a("which_kaon", 19).lit("=").i(p.which_kaon, 2);
        x.x(2).a("pizero_ngamma", 19).lit("=").i(p.pizero_ngamma, 2);
        o.line(x);
    }
    o.line(flags({{"doing_hyd_elast", f.doing_hyd_elast}, {"doing_deuterium", f.doing_deuterium}, {"doing_heavy", f.doing_heavy}}));
    o.line(flags({{"doing_hydpi", f.doing_hydpi}, {"doing_deutpi", f.doing_deutpi}, {"doing_hepi", f.doing_hepi}}));
    o.line(flags({{"doing_hydkaon", f.doing_hydkaon}, {"doing_deutkaon", f.doing_deutkaon}, {"doing_hekaon", f.doing_hekaon}}));
    o.line(flags({{"doing_hydsemi", f.doing_hydsemi}, {"doing_deutsemi", f.doing_deutsemi}, {"do_fermi", f.do_fermi}}));
    o.line(flags({{"doing_hydrho", f.doing_hydrho}, {"doing_deutrho", f.doing_deutrho}, {"doing_herho", f.doing_herho}}));
    {
        Record x;
        x.x(5).x(2).a("mc_smear", 19).lit("=").l(s.sw.mc_smear, 2);
        x.x(2).a("electron_arm", 19).lit("=").i(p.electron_arm, 2);
        x.x(2).a("hadron_arm", 19).lit("=").i(p.hadron_arm, 2);
        o.line(x);
    }
    // FORTRAN-BUG(D-018): original Fortran simc.f:897-898: `write(iun,'(5x,3(2x,a19,''='',l2)))' ... 'deForest_flag',
    // deForest_flag` -- an integer (-1, 0, 1) written with the L edit descriptor (T if non-zero, so -1 and +1 look the
    // same); corrected because the value is not a logical: printed with i2.
    if (s.fortran_compat) {
        o.line(flags({{"using_Eloss", s.sw.using_Eloss}, {"using_Coulomb", s.sw.using_Coulomb}, {"deForest_flag", p.deForest_flag != 0}}));
    } else {
        Record x = flags({{"using_Eloss", s.sw.using_Eloss}, {"using_Coulomb", s.sw.using_Coulomb}});
        x.x(2).a("deForest_flag", 19).lit("=").i(p.deForest_flag, 2);
        o.line(x);
    }
    o.line(flags({{"correct_Eloss", s.sw.correct_Eloss}, {"correct_raster", p.correct_raster}, {"doing_decay", p.doing_decay}}));
    o.line(flags({{"using_E_arm_montecarlo", s.sw.using_E_arm_montecarlo},
                  {"using_P_arm_montecarlo", s.sw.using_P_arm_montecarlo},
                  {"use_benhar_sf", p.use_benhar_sf}}));
    // simc.f:907-909 (the Fortran also lists the retired code 6; rejected at input, D-032)
    if (Is(p.electron_arm, SpectrometerId::SHMS) || Is(p.hadron_arm, SpectrometerId::SHMS)) {
        Record x;
        // FORTRAN-BUG(D-018): original Fortran simc.f:907: `write(iun,'(7x,a19,''='',l2)') 'use_first_cer',use_first_cer`
        // -- an integer written with L; corrected: i2 (fortran_compat keeps L).
        if (s.fortran_compat) x.x(7).a("use_first_cer", 19).lit("=").l(p.use_first_cer != 0, 2);
        else x.x(7).a("use_first_cer", 19).lit("=").i(p.use_first_cer, 2);
        o.line(x);
    }
    { Record x; x.x(7).a("ctau", 11).lit("=").f(p.ctau, 10, 3).a("cm", 4); o.line(x); }
    if (p.use_benhar_sf) { Record x; x.x(7).a("transparency", 12).lit("=").f(p.transparency, 8, 4); o.line(x); }

    // Counters (l.915-921)
    o.list("COUNTERS:");
    { Record x; x.x(12).lit("Ngen (request) = ").i(r.ngen, 10); o.line(x); }
    { Record x; x.x(12).lit("Ntried         = ").i(t.ntried, 10); o.line(x); }
    { Record x; x.x(12).lit("Ncontribute    = ").i(t.ncontribute, 10); o.line(x); }
    { Record x; x.x(12).lit("Nco_no_rad_prot= ").i(t.ncontribute_no_rad_proton, 10); o.line(x); }
    {
        const double pct = (100. * static_cast<double>(t.ncontribute_no_rad_proton)) /
                           fmaxF(static_cast<double>(t.ncontribute), 0.1e0);
        Record x;
        x.x(12).lit("-> %no_rad_prot= ").f(pct, 10, 3);
        o.line(x);
    }
    o.line("");
    { Record x; x.x(1).lit("INTEGRATED WEIGHTS (number of counts in delta/Em cuts!):"); o.line(x); }
    {
        Record x;
        x.x(1).lit("              MeV: wtcontr= ").e(n.wtcontribute / static_cast<double>(t.nevent), 16, 8);
        o.line(x);
    }

    // Radiative corrections (l.924-951)
    o.list("RADIATIVE CORRECTIONS:");
    if (!s.sw.using_rad) {
        Record x;
        x.x(1).a("using_rad", 14).lit("=").l(s.sw.using_rad, 3);
        o.line(x);
    } else {
        {
            Record x;
            // FORTRAN-BUG(D-018): original Fortran simc.f:927: `write(iun,'(4(x,a14,''='',l3))') 'use_expon',use_expon,...`
            // -- an integer written with L; corrected: i3 (fortran_compat keeps L).
            if (s.fortran_compat) x.x(1).a("use_expon", 14).lit("=").l(p.use_expon != 0, 3);
            else x.x(1).a("use_expon", 14).lit("=").i(p.use_expon, 3);
            x.x(1).a("include_hard", 14).lit("=").l(s.sw.include_hard, 3);
            x.x(1).a("calc_spence", 14).lit("=").l(s.sw.calculate_spence, 3);
            o.line(x);
        }
        {
            Record x;
            x.x(1).a("using_rad", 14).lit("=").l(s.sw.using_rad, 3);
            x.x(1).a("use_offshell_rad", 14).lit("=").l(p.use_offshell_rad, 3);
            o.line(x);
        }
        {
            Record x;
            x.x(1).a("rad_flag", 14).lit("=").i(s.sw.rad_flag, 3);
            x.x(1).a("extrad_flag", 14).lit("=").i(s.sw.extrad_flag, 3);
            x.x(1).a("one_tail", 14).lit("=").i(p.one_tail, 3);
            x.x(1).a("intcor_mode", 14).lit("=").i(p.intcor_mode, 3);
            o.line(x);
        }
        { Record x; x.x(1).a("dE_edge_test", 14).lit("=").f(s.params.dE_edge_test, 11, 3); o.line(x); }
        { Record x; x.x(1).a("Egamma_max", 14).lit("=").f(s.rad.Egamma_tot_max, 11, 3); o.line(x); }
        const CentralRad& cr = central.rad;
        o.list("Central Values:");
        o.line(rad4("hardcorfac", {cr.hardcorfac}));
        o.line(rad4("etatzai", {cr.etatzai}));
        o.line(rad4("frac(1:3)", {cr.frac[0], cr.frac[1], cr.frac[2]}));
        o.line(rad4("lambda(1:3)", {cr.lambda[0], cr.lambda[1], cr.lambda[2]}));
        o.line(rad4("bt(1:2)", {cr.bt[0], cr.bt[1]}));
        o.line(rad4("c_int(0:3)", {cr.c_int[0], cr.c_int[1], cr.c_int[2], cr.c_int[3]}));
        o.line(rad4("c_ext(0:3)", {cr.c_ext[0], cr.c_ext[1], cr.c_ext[2], cr.c_ext[3]}));
        o.line(rad4("c(0:3)", {cr.c[0], cr.c[1], cr.c[2], cr.c[3]}));
        o.line(rad4("g_int", {cr.g_int}));
        o.line(rad4("g_ext", {cr.g_ext}));
        o.line(rad4("g(0:3)", {cr.g[0], cr.g[1], cr.g[2], cr.g[3]}));
    }

    // Miscellaneous (l.954-973)
    o.list("MISCELLANEOUS:");
    o.list("Note that central.sigcc is for central delta,theta,phi in both spectrometers");
    o.list(" and may give non-physical kinematics, esp. for Hydrogen");
    o.list("Note also that AVE.sigcc is really AVER.weight (the two arenot exactly equal)");
    o.line(misc("CENTRAL.sigcc", central.sigcc, " "));
    o.line(misc("AVERAGE.sigcc", t.sum_sigcc / static_cast<double>(t.nevent), " "));
    o.line(misc("charge", n.charge, "mC"));
    o.line(misc("targetfac", n.targetfac, " "));
    o.line(misc("luminosity", n.luminosity, "ub^-1"));
    {
        const double h = constants::hbarc / 100000.;
        o.line(misc("luminosity", n.luminosity * (h * h), "GeV^2"));
    }
    o.line(misc("genvol", n.genvol, " "));
    o.line(misc("normfac", n.normfac, " "));
    if (f.doing_heavy) {
        Record x;
        x.x(12).lit("Theory file:  ").lit(s.theory.file);
        o.line(x);
    }
    if (r.random_seed != 0) { Record x; x.x(15).lit("Random Seed = ").i(r.random_seed, 10); o.line(x); }
    if (!r.random_state_file.empty() && r.random_state_file.find_first_not_of(' ') != std::string::npos) {
        Record x;
        x.lit(" Random State Save File:  ").lit(r.random_state_file.substr(0, r.random_state_file.find(' ')));
        o.line(x);
    }

    // Resolution summary (l.977-986)
    o.list("RECON SUMMARY:            Ave.Error   Resolution");
    const TwoArmSums& ae = n.aveerr;
    const TwoArmSums& rs = n.resol;
    o.line(recon("Electron arm: delta =", units::kPermillePerPercent * ae.e.delta, units::kPermillePerPercent * rs.e.delta, "x10^-3"));
    o.line(recon("xptar =", units::kMradPerRad * ae.e.xptar, units::kMradPerRad * rs.e.xptar, "mr"));
    o.line(recon("yptar =", units::kMradPerRad * ae.e.yptar, units::kMradPerRad * rs.e.yptar, "mr"));
    o.line(recon("ytar  =", units::kMmPerCm * ae.e.ytar, units::kMmPerCm * rs.e.ytar, "mm"));
    o.line(recon("Hadron arm:   delta =", units::kPermillePerPercent * ae.p.delta, units::kPermillePerPercent * rs.p.delta, "x10-3"));
    o.line(recon("xptar =", units::kMradPerRad * ae.p.xptar, units::kMradPerRad * rs.p.xptar, "mr"));
    o.line(recon("yptar =", units::kMradPerRad * ae.p.yptar, units::kMradPerRad * rs.p.yptar, "mr"));
    o.line(recon("ytar  =", units::kMmPerCm * ae.p.ytar, units::kMmPerCm * rs.p.ytar, "mm"));

    // Input spectrometer limits (l.991-1003)
    o.list("Input Spectrometer Limits:");
    o.line(spedge("SPedge.e.delta.min/max", s.SPedge.e.delta, "%"));
    o.line(spedge("SPedge.e.yptar.min/max", s.SPedge.e.yptar, "rad"));
    o.line(spedge("SPedge.e.xptar.min/max", s.SPedge.e.xptar, "rad"));
    o.line(spedge("SPedge.p.delta.min/max", s.SPedge.p.delta, "%"));
    o.line(spedge("SPedge.p.yptar.min/max", s.SPedge.p.yptar, "rad"));
    o.line(spedge("SPedge.p.xptar.min/max", s.SPedge.p.xptar, "rad"));

    // Generated quantities (l.1010-1050)
    const GenLimits& g = s.gen;
    o.list("Limiting VERTEX values (vertex.e/p.*,Em,Pm,Trec)");
    o.list("   USED limits are gen.e/p.*, and VERTEXedge.Em,Pm,Trec");
    o.line(limHeader1());
    o.line(limHeader2());
    o.line(lim("E arm  delta", g.e.delta.min, g.e.delta.max, c.gen.e.delta.lo, c.gen.e.delta.hi, "%"));
    o.line(lim("E arm  yptar", g.e.yptar.min * 1000., g.e.yptar.max * 1000., c.gen.e.yptar.lo * 1000., c.gen.e.yptar.hi * 1000., "mr"));
    o.line(lim("E arm  xptar", g.e.xptar.min * 1000., g.e.xptar.max * 1000., c.gen.e.xptar.lo * 1000., c.gen.e.xptar.hi * 1000., "mr"));
    o.line(lim("P arm  delta", g.p.delta.min, g.p.delta.max, c.gen.p.delta.lo, c.gen.p.delta.hi, "%"));
    o.line(lim("P arm  yptar", g.p.yptar.min * 1000., g.p.yptar.max * 1000., c.gen.p.yptar.lo * 1000., c.gen.p.yptar.hi * 1000., "mr"));
    o.line(lim("P arm  xptar", g.p.xptar.min * 1000., g.p.xptar.max * 1000., c.gen.p.xptar.lo * 1000., c.gen.p.xptar.hi * 1000., "mr"));
    o.line(lim("sumEgen", g.sumEgen.min, g.sumEgen.max, c.gen.sumEgen.lo, c.gen.sumEgen.hi, "MeV"));
    const EdgeTrue& V = s.VERTEXedge;
    o.line(lim("Trec", V.Trec.min, V.Trec.max, c.vertex.Trec.lo, c.vertex.Trec.hi, "MeV"));
    o.line(lim("Em", V.Em.min, V.Em.max, c.vertex.Em.lo, c.vertex.Em.hi, "MeV"));
    o.line(lim("Pm", V.Pm.min, V.Pm.max, c.vertex.Pm.lo, c.vertex.Pm.hi, "MeV/c"));
    if ((f.doing_deuterium || f.doing_pion || f.doing_kaon || f.doing_delta) && s.sw.using_rad)
        o.list("      *** NOTE: sumEgen.min only used in GENERATE_RAD");

    // True quantities (l.1053-1075)
    const EdgeTrue& E = s.edge;
    o.list("Limiting ORIGINAL values: orig.e/p.*,Em,Pm,Trec (no edge.* limits for Pm,Trec)");
    o.line(limHeader1());
    o.line(limHeader2());
    o.line(lim("E arm   E", E.e.E.min, E.e.E.max, c.tru.e.E.lo, c.tru.e.E.hi, "MeV"));
    o.line(lim("E arm  yptar", E.e.yptar.min * 1000., E.e.yptar.max * 1000., c.tru.e.yptar.lo * 1000., c.tru.e.yptar.hi * 1000., "mr"));
    o.line(lim("E arm  xptar", E.e.xptar.min * 1000., E.e.xptar.max * 1000., c.tru.e.xptar.lo * 1000., c.tru.e.xptar.hi * 1000., "mr"));
    o.line(lim("P arm      E", E.p.E.min, E.p.E.max, c.tru.p.E.lo, c.tru.p.E.hi, "MeV"));
    o.line(lim("P arm  yptar", E.p.yptar.min * 1000., E.p.yptar.max * 1000., c.tru.p.yptar.lo * 1000., c.tru.p.yptar.hi * 1000., "mr"));
    o.line(lim("P arm  xptar", E.p.xptar.min * 1000., E.p.xptar.max * 1000., c.tru.p.xptar.lo * 1000., c.tru.p.xptar.hi * 1000., "mr"));
    o.line(lim("Em", fmaxF(-999999.999e0, E.Em.min), fminF(999999.999e0, E.Em.max), c.tru.Em.lo, c.tru.Em.hi, "MeV"));
    o.line(lim("Pm", 0., 0., c.tru.Pm.lo, c.tru.Pm.hi, "MeV"));
    o.line(lim("Trec", 0., 0., c.tru.Trec.lo, c.tru.Trec.hi, "MeV"));

    // Spectrometer quantities (l.1078-1097)
    o.list("Limiting SPECTROMETER values");
    o.line(limHeader1());
    o.line(limHeader2());
    const auto& SE = s.SPedge;
    o.line(lim("E arm delta", SE.e.delta.min, SE.e.delta.max, c.SP.e.delta.lo, c.SP.e.delta.hi, "%"));
    o.line(lim("E arm  yptar", SE.e.yptar.min * 1000., SE.e.yptar.max * 1000., c.SP.e.yptar.lo * 1000., c.SP.e.yptar.hi * 1000., "mr"));
    o.line(lim("E arm  xptar", SE.e.xptar.min * 1000., SE.e.xptar.max * 1000., c.SP.e.xptar.lo * 1000., c.SP.e.xptar.hi * 1000., "mr"));
    o.line(lim("P arm  delta", SE.p.delta.min, SE.p.delta.max, c.SP.p.delta.lo, c.SP.p.delta.hi, "%"));
    o.line(lim("P arm  yptar", SE.p.yptar.min * 1000., SE.p.yptar.max * 1000., c.SP.p.yptar.lo * 1000., c.SP.p.yptar.hi * 1000., "mr"));
    o.line(lim("P arm  xptar", SE.p.xptar.min * 1000., SE.p.xptar.max * 1000., c.SP.p.xptar.lo * 1000., c.SP.p.xptar.hi * 1000., "mr"));

    // Radiation (l.1100-1113)
    if (s.sw.using_rad) {
        o.list("Limiting RADIATION values CONTRIBUTING to the (Em,Pm) distributions:");
        o.line(limHeader1());
        o.line(limHeader2());
        o.line(lim("Egamma(1)", 0., s.rad.Egamma1_max, c.rad.Egamma[0].lo, c.rad.Egamma[0].hi, "MeV"));
        o.line(lim("Egamma(2)", 0., s.rad.Egamma2_max, c.rad.Egamma[1].lo, c.rad.Egamma[1].hi, "MeV"));
        o.line(lim("Egamma(3)", 0., s.rad.Egamma3_max, c.rad.Egamma[2].lo, c.rad.Egamma[2].hi, "MeV"));
        o.line(lim("Egamma_total", 0., s.rad.Egamma_tot_max, c.rad.Egamma_total.lo, c.rad.Egamma_total.hi, "MeV"));
    }

    // Slops (l.1116-1137)
    const Slops& sl = t.slop;
    o.list("ACTUAL and LIMITING SLOP values used/obtained:");
    { Record x; x.t(25).a("__used__", 10).a("__min__", 10).a("__max__", 10); o.line(x); }
    o.line(slopLine("slop.MC  ", "E arm delta", sl.MC.e.delta.used, sl.MC.e.delta.lo, sl.MC.e.delta.hi, "%"));
    o.line(slopLine(" ", "E arm yptar", sl.MC.e.yptar.used * 1000., sl.MC.e.yptar.lo * 1000., sl.MC.e.yptar.hi * 1000., "mr"));
    o.line(slopLine(" ", "E arm xptar", sl.MC.e.xptar.used * 1000., sl.MC.e.xptar.lo * 1000., sl.MC.e.xptar.hi * 1000., "mr"));
    o.line(slopLine(" ", "P arm delta", sl.MC.p.delta.used, sl.MC.p.delta.lo, sl.MC.p.delta.hi, "%"));
    o.line(slopLine(" ", "P arm yptar", sl.MC.p.yptar.used * 1000., sl.MC.p.yptar.lo * 1000., sl.MC.p.yptar.hi * 1000., "mr"));
    o.line(slopLine(" ", "P arm xptar", sl.MC.p.xptar.used * 1000., sl.MC.p.xptar.lo * 1000., sl.MC.p.xptar.hi * 1000., "mr"));
    o.line(slopLine("slop.total", "Em", sl.total.Em.used, sl.total.Em.lo, sl.total.Em.hi, "MeV"));
    o.line(slopLine(" ", "Pm", 0., sl.total.Pm.lo, sl.total.Pm.hi, "MeV/c"));

    // write(iun,'(/)')
    o.line("");
    o.line("");
    return o.str();
}


}  // namespace driver
}  // namespace simc
