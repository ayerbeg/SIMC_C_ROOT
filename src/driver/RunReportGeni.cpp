// src/driver/RunReportGeni.cpp — part of the report writer (RunReport.h): the .geni text, the per-spectrometer stop
// counters (where the rays were lost: slits, magnets, detectors)
// (Fortran: simc.f l.462-532).
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

namespace {
long cnt(const RunReport& r, const std::string& spec, const std::string& key) {
    for (const auto& sc : r.stops)
        if (sc.spectrometer == spec) {
            auto it = sc.counts.find(key);
            if (it == sc.counts.end()) throw std::runtime_error("FormatGeni: no counter '" + key + "' for " + spec);
            return it->second;
        }
    throw std::runtime_error("FormatGeni: no stop counters for " + spec);
}
bool has(const RunReport& r, const std::string& spec) {
    for (const auto& sc : r.stops)
        if (sc.spectrometer == spec) return true;
    return false;
}
}  // namespace

// Fortran: simc.f l.462-532
std::string FormatGeni(const RunSetup& s, const RunReport& r) {
    Out o;
    const int ea = s.params.electron_arm, ha = s.params.hadron_arm;
    auto row = [&](const std::string& label, const std::string& spec, std::initializer_list<const char*> keys) {
        Record x = List(label);
        for (const char* k : keys) x.listInt(cnt(r, spec, k));
        o.line(x);
    };
    // The Fortran blocks (by arm code) are written for the Fortran spectrometers; an arm served by another ISpectrometer
    // plug-in (e.g. examples/ideal_box) gets a generic block at the end.
    auto uses = [&](SpectrometerId id) { return Is(ea, id) || Is(ha, id); };   // D-032 spectrometer table
    if (uses(SpectrometerId::HMS) && has(r, "hms")) {
        row("HMS Trials:           ", "hms", {"trials"});
        row("Slit hor/vert/corners ", "hms", {"slit_hor", "slit_vert", "slit_oct"});
        // FORTRAN-BUG(D-037 F2): original Fortran simc.f:462-532: hSTOP_coll (mc_hms.f:218) and sSTOP_coll
        // (mc_shms.f:518) are written to no output, so with the collimator simulation the .geni stop counts do not add
        // up to the trials (1253 of 3932 HMS pions in reference/runs/p6f2_pi_hms_cap). Corrected: an extra line with
        // the collimator stops (absorbed or ranged out) whenever the collimator simulation is on. fortran_compat writes
        // the Fortran .geni.
        if (!s.fortran_compat && (s.params.using_HMScoll || cnt(r, "hms", "coll") != 0))
            row("Collimator stops (pi/mu)", "hms", {"coll"});
        row("Q1 entrance/mid/exit  ", "hms", {"Q1_in", "Q1_mid", "Q1_out"});
        row("Q2 entrance/mid/exit  ", "hms", {"Q2_in", "Q2_mid", "Q2_out"});
        row("Q3 entrance/mid/exit  ", "hms", {"Q3_in", "Q3_mid", "Q3_out"});
        row("Dipole entrance/exit  ", "hms", {"D1_in", "D1_out"});
        row("Events reaching hut   ", "hms", {"hut"});
        row("DC1, DC2, Scin, Cal   ", "hms", {"dc1", "dc2", "scin", "cal"});
        row("Successes             ", "hms", {"successes"});
        o.line("");
    }
    if (uses(SpectrometerId::SOS) && has(r, "sos")) {   // simc.f:475-486
        row("SOS Trials:           ", "sos", {"trials"});
        // FORTRAN-BUG(D-011, D-061 F8): original Fortran simc.f:477: `write(7,*) 'Slit hor/vert/corners ',sSTOP_slit_vert,
        // sSTOP_slit_hor,sSTOP_slit_oct`; corrected because the label (and the HMS/SHMS lines) put slit_hor first: the
        // SOS horizontal and vertical slit stops were printed swapped. fortran_compat writes the Fortran order.
        row("Slit hor/vert/corners ", "sos",
            {s.fortran_compat ? "slit_vert" : "slit_hor", s.fortran_compat ? "slit_hor" : "slit_vert", "slit_oct"});
        row("Quad entrance/mid/exit", "sos", {"quad_in", "quad_mid", "quad_out"});
        row("D1 entrance/exit      ", "sos", {"bm01_in", "bm01_out"});
        row("D2 entrance/exit      ", "sos", {"bm02_in", "bm02_out"});
        row("Vacuum exit           ", "sos", {"exit"});
        row("Events reaching hut   ", "sos", {"hut"});
        row("DC1, DC2, Scin, Cal   ", "sos", {"dc1", "dc2", "scin", "cal"});
        row("Successes             ", "sos", {"successes"});
        o.line("");
    }
    // simc.f:487-509: HRS-R then HRS-L. The "S1, S2, Cal" label contains a TAB (simc.f:496, 508; kept).
    // FORTRAN-BUG(D-061 F8): the Fortran HRS blocks have no "Successes" line and the HRS-L block no trailing blank line;
    // corrected: both added (as in the HMS/SOS/SHMS blocks). fortran_compat writes the Fortran blocks.
    auto hrsBlock = [&](const char* key, const char* trials, bool blank) {
        row(trials, key, {"trials"});
        // FORTRAN-BUG(D-011, D-061 F8): original Fortran simc.f:489/501: `write(7,*) 'Slit hor/vert         ',rSTOP_slit_vert,
        // rSTOP_slit_hor`; corrected because the label puts slit_hor first. fortran_compat writes the Fortran order.
        row("Slit hor/vert         ", key,
            {s.fortran_compat ? "slit_vert" : "slit_hor", s.fortran_compat ? "slit_hor" : "slit_vert"});
        row("Q1 entrance/mid/exit  ", key, {"Q1_in", "Q1_mid", "Q1_out"});
        row("Q2 entrance/mid/exit  ", key, {"Q2_in", "Q2_mid", "Q2_out"});
        row("Dipole entrance/exit  ", key, {"D1_in", "D1_out"});
        row("Q3 entrance/mid/exit  ", key, {"Q3_in", "Q3_mid", "Q3_out"});
        row("Events reaching hut   ", key, {"hut"});
        row("VDC1, VDC2            ", key, {"dc1", "dc2"});
        row("S1, S2, Cal\t    ", key, {"s1", "s2", "cal"});
        if (!s.fortran_compat) row("Successes             ", key, {"successes"});
        if (blank || !s.fortran_compat) o.line("");
    };
    if (uses(SpectrometerId::HRSR) && has(r, "hrsr")) hrsBlock("hrsr", "HRSr Trials:          ", true);
    if (uses(SpectrometerId::HRSL) && has(r, "hrsl")) hrsBlock("hrsl", "HRSl Trials:          ", false);
    if (uses(SpectrometerId::SHMS) && has(r, "shms")) {   // simc.f:510-511 also lists the retired code 6
        row("SHMS Trials:          ", "shms", {"trials"});
        row("HB phys entrance/mag entr/mag exit/phys exit  ", "shms", {"hb_in", "hb_men", "hb_mex", "hb_out"});
        row("Slit hor/vert/corners ", "shms", {"slit_hor", "slit_vert", "slit_oct"});
        if (!s.fortran_compat && (s.params.using_SHMScoll || cnt(r, "shms", "coll") != 0))   // FORTRAN-BUG(D-037 F2)
            row("Collimator stops (pi/mu)", "shms", {"coll"});
        row("Q1 phys entrance/mag entr/mid/mag exit/phys exit  ", "shms", {"q1_in", "q1_men", "q1_mid", "q1_mex", "q1_out"});
        row("Q2 phys entrance/mag entr/mid/mag exit/phys exit  ", "shms", {"q2_in", "q2_men", "q2_mid", "q2_mex", "q2_out"});
        // FORTRAN-BUG(D-011): original Fortran simc.f:519-520: `write(7,*) 'Q3 phys entrance/mag entr/mid/mag exit/phys
        // exit ',shmsSTOP_q3_in,shmsSTOP_q3_men,shmsSTOP_q2_mid,shmsSTOP_q3_mex,shmsSTOP_q3_out`; corrected because the
        // Q3 "mid" column must show shmsSTOP_q3_mid (the Q2 value was printed twice and q3_mid never).
        row("Q3 phys entrance/mag entr/mid/mag exit/phys exit       ", "shms",
            {"q3_in", "q3_men", s.fortran_compat ? "q2_mid" : "q3_mid", "q3_mex", "q3_out"});
        row("D1 entrance/flare/mid 1-2   ", "shms", {"D1_in", "D1_flr", "D1_mid1", "D1_mid2"});
        row("D1 mid 3-5            ", "shms", {"D1_mid3", "D1_mid4", "D1_mid5"});
        row("D1 mid 6-7/mag exit/phys exit         ", "shms", {"D1_mid6", "D1_mid7", "D1_mex", "D1_out"});
        row("Events reaching hut   ", "shms", {"hut"});
        // FORTRAN-BUG(D-011): original Fortran simc.f:526: `write(7,*) 'DC1, DC2, Scin, Cal   ',shmsSTOP_dc1,shmsSTOP_dc2`;
        // corrected because the label names four counters but two are printed (scin/cal are on the next line).
        row(s.fortran_compat ? "DC1, DC2, Scin, Cal   " : "DC1, DC2              ", "shms", {"dc1", "dc2"});
        row("S1, S2, S3, Cal       ", "shms", {"s1", "s2", "s3", "cal"});
        row("Successes             ", "shms", {"successes"});
        o.line("");
    }
    // simc.f:531-535: the calorimeter (either side); 'Calo Trials' counts mc_calo calls, i.e. PHOTONS for a
    // pi0 hadron arm. No blank line after it (the Fortran block has none).
    if ((uses(SpectrometerId::CALO_HMS_SIDE) || uses(SpectrometerId::CALO_SHMS_SIDE)) && has(r, "calo")) {
        row("Calo Trials:          ", "calo", {"trials"});
        row("Extent hor/vert       ", "calo", {"slit_hor", "slit_vert"});
        row("Successes             ", "calo", {"successes"});
    }
    // Every non-Fortran spectrometer: one "key value" line per counter.
    for (const auto& sc : r.stops) {
        if (sc.spectrometer == "hms" || sc.spectrometer == "shms" || sc.spectrometer == "sos" || sc.spectrometer == "hrsr" ||
            sc.spectrometer == "hrsl" || sc.spectrometer == "calo")
            continue;
        o.list(sc.spectrometer + " stop counters (" + sc.role + " arm, SIMC++ plug-in):");
        for (const auto& kv : sc.counts) {
            Record x = List(" " + kv.first);
            x.listInt(kv.second);
            o.line(x);
        }
        o.line("");
    }
    return o.str();
}


}  // namespace driver
}  // namespace simc
