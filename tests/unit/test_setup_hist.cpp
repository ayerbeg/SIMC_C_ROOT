// tests/unit/test_setup_hist.cpp
//
// Task P2.1: "setup vs .hist". Builds the RunSetup from a Fortran reference .inp (InpReader ->
// BuildRunSetup) and compares every derived quantity the Fortran prints in the .hist
// (simc.f:report l.646-1140) that the setup can compute. The C++ value is rendered with the SAME
// Fortran edit descriptor (f15.4, e12.6, f12.3, ...) and must give the identical text.
//
// The energy-loss numbers (target_init, extreme_trip_thru_target) come from the default provider of
// BuildRunSetup, material::HallCCryoTarget (P2.2). The real inputs are read by InpReader with the CTP
// float-literal semantics (D-018). Quantities of the event loop / central event (calculate_central,
// genvol, normfac, counters, contrib lo/hi, recon summary, central radiative constants) are DEFERRED.
//
// Usage: test_setup_hist <run.inp> <run.hist> [--data-dir DIR]
// Exit: 0 all checks pass, 1 failure, 77 reference missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "simc/core/SimcConstants.h"
#include "simc/input/InpReader.h"
#include "simc/output/SummaryJson.h"
#include "simc/setup/SetupBuilder.h"

using namespace simc;

namespace {

// ------------------------------------------------------------------ Fortran edit descriptors
std::string fmtF(double v, int w, int d) {
    char b[64];
    std::snprintf(b, sizeof b, "%*.*f", w, d, v);
    std::string s = b;
    if (static_cast<int>(s.size()) > w) return std::string(w, '*');
    return s;
}
// gfortran Ew.d: 0.dddddE+xx (d significant digits after "0.")
std::string fmtE(double v, int w, int d) {
    char b[64];
    std::string out;
    if (v == 0.0) {
        out = std::string(std::signbit(v) ? "-" : "") + "0." + std::string(d, '0') + "E+00";
    } else {
        std::snprintf(b, sizeof b, "%.*E", d - 1, v);   // m.mmmmmE+xx
        std::string s = b;
        const bool neg = s[0] == '-';
        if (neg) s = s.substr(1);
        const auto epos = s.find('E');
        std::string mant = s.substr(0, epos);
        int ex = std::atoi(s.c_str() + epos + 1) + 1;
        mant.erase(1, 1);   // drop '.'
        char eb[16];
        std::snprintf(eb, sizeof eb, "E%c%02d", ex < 0 ? '-' : '+', std::abs(ex));
        out = std::string(neg ? "-" : "") + "0." + mant + eb;
    }
    if (static_cast<int>(out.size()) > w) return std::string(w, '*');
    return out;
}
std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r") - a + 1);
}
std::string collapse(const std::string& s) {
    std::istringstream is(s);
    std::string t, out;
    while (is >> t) out += (out.empty() ? "" : " ") + t;
    return out;
}
std::vector<std::string> tokens(const std::string& s) {
    std::istringstream is(s);
    std::vector<std::string> v;
    std::string t;
    while (is >> t) v.push_back(t);
    return v;
}
bool isNumber(const std::string& t) {
    char* e = nullptr;
    std::strtod(t.c_str(), &e);
    return e != t.c_str() && *e == '\0';
}

// ------------------------------------------------------------------ .hist access
class Hist {
public:
    bool load(const std::string& path) {
        std::ifstream f(path);
        if (!f) return false;
        std::string l, sec = "head";
        while (std::getline(f, l)) {
            const std::string c = collapse(l);
            static const std::vector<std::pair<std::string, std::string>> heads = {
                {"KINEMATICS:", "kinematics"}, {"VALUES FOR \"CENTRAL\" EVENT:", "central"},
                {"TARGET specs:", "target"},   {"FLAGS:", "flags"},
                {"COUNTERS:", "counters"},     {"RADIATIVE CORRECTIONS:", "radiative"},
                {"MISCELLANEOUS:", "misc"},    {"RECON SUMMARY:", "recon"},
                {"Input Spectrometer Limits:", "spedge"}, {"Limiting VERTEX", "vertex"},
                {"Limiting ORIGINAL", "original"}, {"Limiting SPECTROMETER", "spectrometer"},
                {"Limiting RADIATION", "radiation"}, {"ACTUAL and LIMITING SLOP", "slop"}};
            bool header = false;
            for (const auto& h : heads)
                if (c.rfind(h.first, 0) == 0) { sec = h.second; header = true; }
            if (!header) lines_.push_back({sec, c});
        }
        return true;
    }
    /// "key = v0 v1 ..." (key may contain blanks; a18/a12 labels). occurrence is 1-based.
    bool kv(const std::string& sec, const std::string& key, int idx, std::string& out, int occ = 1) const {
        int n = 0;
        for (const auto& [s, c] : lines_) {
            if (s != sec) continue;
            const auto eq = c.find('=');
            if (eq == std::string::npos) continue;
            if (trim(c.substr(0, eq)) != key) continue;
            if (++n != occ) continue;
            auto v = tokens(c.substr(eq + 1));
            if (idx >= static_cast<int>(v.size())) return false;
            out = v[idx];
            return true;
        }
        return false;
    }
    /// "name= value" items (FLAGS / RADIATIVE blocks, several per line) and 9911 "key = value" pairs.
    bool item(const std::string& sec, const std::string& key, std::string& out, int occ = 1) const {
        int n = 0;
        for (const auto& [s, c] : lines_) {
            if (s != sec) continue;
            auto t = tokens(c);
            for (size_t i = 0; i < t.size(); ++i) {
                std::string name, val;
                if (t[i] == "=" && i >= 1 && i + 1 < t.size()) { name = t[i - 1]; val = t[i + 1]; }
                else if (t[i].size() > 1 && t[i].back() == '=' && i + 1 < t.size()) {
                    name = t[i].substr(0, t[i].size() - 1); val = t[i + 1];
                }
                else if (t[i].find('=') != std::string::npos && t[i].back() != '=' && t[i].front() != '=') {
                    const auto e = t[i].find('=');   // "name=value" glued
                    name = t[i].substr(0, e); val = t[i].substr(e + 1);
                }
                if (name == key && ++n == occ) { out = val; return true; }
            }
        }
        return false;
    }
    /// Table rows (formats 9912, 9917, 9918): leading label tokens, then numbers.
    bool row(const std::string& sec, const std::string& label, int idx, std::string& out) const {
        for (const auto& [s, c] : lines_) {
            if (s != sec) continue;
            auto t = tokens(c);
            size_t i = 0;
            std::string lab;
            while (i < t.size() && !isNumber(t[i])) lab += (lab.empty() ? "" : " ") + t[i++];
            if (lab.rfind("slop.MC ", 0) == 0) lab = lab.substr(8);
            if (lab.rfind("slop.total ", 0) == 0) lab = lab.substr(11);
            if (lab != label) continue;
            if (i + idx >= t.size()) return false;
            out = t[i + idx];
            return true;
        }
        return false;
    }

private:
    std::vector<std::pair<std::string, std::string>> lines_;
};

// ------------------------------------------------------------------ checks
enum class Dep { None, Loss, Provider };

struct Result {
    int pass{0}, fail{0}, deferred{0};
    std::vector<std::string> deferred_items;
};

class Checker {
public:
    Checker(const Hist& h, bool loss_ok, bool provider) : h_(h), loss_ok_(loss_ok), provider_(provider) {}
    Result res;

    bool available(Dep d, const std::string& name) {
        if (d == Dep::None) return true;
        if (d == Dep::Loss && loss_ok_) return true;
        if (d == Dep::Provider && provider_) return true;
        ++res.deferred;
        res.deferred_items.push_back(name + " [needs ITargetLossProvider]");
        return false;
    }
    void compare(const std::string& name, bool found, const std::string& fortran, const std::string& cpp) {
        if (!found) {
            ++res.fail;
            std::printf("  FAIL %-40s not found in .hist\n", name.c_str());
            return;
        }
        if (trim(fortran) == trim(cpp)) {
            ++res.pass;
            std::printf("  ok   %-40s %s\n", name.c_str(), trim(cpp).c_str());
        } else {
            ++res.fail;
            std::printf("  FAIL %-40s Fortran '%s'  C++ '%s'\n", name.c_str(), trim(fortran).c_str(), trim(cpp).c_str());
        }
    }
    // kinds of lookups
    void kv(const std::string& sec, const std::string& key, int idx, const std::string& cpp, Dep d = Dep::None,
            int occ = 1) {
        const std::string name = sec + ": " + key + (idx ? "[" + std::to_string(idx) + "]" : "") +
                                 (occ > 1 ? "#" + std::to_string(occ) : "");
        if (!available(d, name)) return;
        std::string f;
        compare(name, h_.kv(sec, key, idx, f, occ), f, cpp);
    }
    void item(const std::string& sec, const std::string& key, const std::string& cpp, Dep d = Dep::None, int occ = 1) {
        const std::string name = sec + ": " + key + (occ > 1 ? "#" + std::to_string(occ) : "");
        if (!available(d, name)) return;
        std::string f;
        compare(name, h_.item(sec, key, f, occ), f, cpp);
    }
    void row(const std::string& sec, const std::string& label, int idx, const std::string& cpp, Dep d = Dep::None) {
        const std::string name = sec + ": " + label + "[" + std::to_string(idx) + "]";
        if (!available(d, name)) return;
        std::string f;
        compare(name, h_.row(sec, label, idx, f), f, cpp);
    }
    void defer(const std::string& what) {
        ++res.deferred;
        res.deferred_items.push_back(what);
    }

private:
    const Hist& h_;
    bool loss_ok_, provider_;
};

std::string L(bool b) { return b ? "T" : "F"; }

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s run.inp run.hist [--data-dir DIR]\n", argv[0]);
        return 1;
    }
    const std::string inp = argv[1], histfile = argv[2];
    std::string data_dir;
    for (int i = 3; i + 1 < argc; i += 2)
        if (std::string(argv[i]) == "--data-dir") data_dir = argv[i + 1];
    Hist H;
    if (!H.load(histfile)) {
        std::printf("SKIP: %s missing\n", histfile.c_str());
        return 77;
    }
    SimcParameters p;
    InpReader rd(p);
    if (!rd.loadRun(inp)) {
        std::printf("SKIP: cannot read %s\n", inp.c_str());
        return 77;
    }
    SetupOptions opt;
    opt.fortran_compat = true;   // D-031: Fortran references -> exact Fortran behaviour
    opt.data_dir = data_dir;
    const RunSetup s = BuildRunSetup(p, opt);   // default provider: material::HallCCryoTarget
    std::printf("setup vs .hist: %s\n  inp %s\n", histfile.c_str(), inp.c_str());
    for (const auto& m : s.stops) std::printf("  setup STOP: %s\n", m.c_str());
    if (!s.ok()) return 1;

    const RunSwitches& sw = s.sw;
    const ReactionFlags& f = s.flags;
    const TargetInfo& t = s.targ;
    const double degrad = constants::degrad;
    Checker C(H, s.target_losses_from_provider, s.target_losses_from_provider);
    const Dep LO = Dep::Loss, PR = Dep::Provider;
    // Eloss ave rows: exact zero if !using_Eloss
    const Dep ELAVE = sw.using_Eloss ? PR : Dep::None;
    const Dep MUSC = sw.mc_smear ? PR : Dep::None;

    // ---- KINEMATICS (simc.f:791-812)
    C.kv("kinematics", "Ebeam", 0, fmtF(s.Ebeam, 15, 4));
    C.kv("kinematics", "(dE/E)beam", 0, fmtF(s.dEbeam / s.Ebeam, 15, 4));
    C.kv("kinematics", "x-width", 0, fmtF(s.gen.xwid, 15, 4));
    C.kv("kinematics", "y-width", 0, fmtF(s.gen.ywid, 15, 4));
    C.kv("kinematics", "fr_pattern", 0, std::to_string(t.fr_pattern));
    C.kv("kinematics", "fr1", 0, fmtF(t.fr1, 15, 4));
    C.kv("kinematics", "fr2", 0, fmtF(t.fr2, 15, 4));
    C.kv("kinematics", "angle", 0, fmtF(s.spec.e.theta * degrad, 15, 4));
    C.kv("kinematics", "angle", 1, fmtF(s.spec.p.theta * degrad, 15, 4));
    C.kv("kinematics", "momentum", 0, fmtF(s.spec.e.P, 15, 4));
    C.kv("kinematics", "momentum", 1, fmtF(s.spec.p.P, 15, 4));
    C.kv("kinematics", "x offset", 0, fmtF(s.spec.e.offset.x, 15, 4));
    C.kv("kinematics", "x offset", 1, fmtF(s.spec.p.offset.x, 15, 4));
    C.kv("kinematics", "y offset", 0, fmtF(s.spec.e.offset.y, 15, 4));
    C.kv("kinematics", "y offset", 1, fmtF(s.spec.p.offset.y, 15, 4));
    C.kv("kinematics", "z offset", 0, fmtF(s.spec.e.offset.z, 15, 4));
    C.kv("kinematics", "z offset", 1, fmtF(s.spec.p.offset.z, 15, 4));
    C.kv("kinematics", "xptar offset", 0, fmtF(s.spec.e.offset.xptar, 15, 4));   // printed in rad (label 'mr')
    C.kv("kinematics", "xptar offset", 1, fmtF(s.spec.p.offset.xptar, 15, 4));
    C.kv("kinematics", "yptar offset", 0, fmtF(s.spec.e.offset.yptar, 15, 4));
    C.kv("kinematics", "yptar offset", 1, fmtF(s.spec.p.offset.yptar, 15, 4));
    C.defer("central event: delta/xptar/yptar, Q2, q, nu, recon Em/Pm/W/MM [calculate_central, P2.3]");

    // ---- TARGET (simc.f:824-869), format 9911 e12.6 / 9912 f15.5 / 9913 f15.5
    C.item("target", "A", fmtE(t.A, 12, 6));
    C.item("target", "Z", fmtE(t.Z, 12, 6));
    C.item("target", "mass", fmtE(t.mass_amu, 12, 6), Dep::None, 1);
    C.item("target", "mass", fmtE(t.M, 12, 6), Dep::None, 2);
    C.item("target", "Mrec", fmtE(t.mrec_amu, 12, 6), Dep::None, 1);
    C.item("target", "Mrec", fmtE(t.Mrec, 12, 6), Dep::None, 2);
    C.item("target", "Mtar_struc", fmtE(t.Mtar_struck, 12, 6));
    C.item("target", "Mrec_struc", fmtE(t.Mrec_struck, 12, 6));
    C.item("target", "rho", fmtE(t.rho, 12, 6));
    C.item("target", "thick", fmtE(t.thick, 12, 6));
    C.item("target", "angle", fmtE(t.angle * degrad, 12, 6));
    C.item("target", "abundancy", fmtE(t.abundancy, 12, 6));
    C.item("target", "X0", fmtE(t.X0, 12, 6));
    C.item("target", "X0_cm", fmtE(t.X0_cm, 12, 6));
    C.item("target", "length", fmtE(t.length, 12, 6));
    C.item("target", "zoffset", fmtE(t.zoffset, 12, 6));
    C.item("target", "xoffset", fmtE(t.xoffset, 12, 6));
    C.item("target", "yoffset", fmtE(t.yoffset, 12, 6));
    C.row("target", "Coulomb", 0, fmtF(t.Coulomb.ave, 15, 5));
    C.row("target", "Coulomb", 1, fmtF(t.Coulomb.min, 15, 5));
    C.row("target", "Coulomb", 2, fmtF(t.Coulomb.max, 15, 5));
    const char* elossNames[3] = {"Eloss_beam", "Eloss_e", "Eloss_p"};
    const char* teffNames[3] = {"teff_beam", "teff_e", "teff_p"};
    const char* muscNames[3] = {"musc_max_beam", "musc_max_e", "musc_max_p"};
    for (int i = 0; i < 3; ++i) {
        C.row("target", elossNames[i], 0, fmtF(t.Eloss[i].ave, 15, 5), ELAVE);
        C.row("target", elossNames[i], 1, fmtF(t.Eloss[i].min, 15, 5), sw.using_Eloss ? PR : Dep::None);
        C.row("target", elossNames[i], 2, fmtF(t.Eloss[i].max, 15, 5), sw.using_Eloss ? PR : Dep::None);
    }
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) {
            const double v = k == 0 ? t.teff[i].ave : (k == 1 ? t.teff[i].min : t.teff[i].max);
            C.row("target", teffNames[i], k, fmtF(v, 15, 5), PR);
        }
    C.row("target", "musc_nsig_max", 0, fmtF(t.musc_nsig_max, 15, 5), PR);
    for (int i = 0; i < 3; ++i) C.row("target", muscNames[i], 0, fmtF(t.musc_max[i] * 1000., 15, 5), MUSC);

    // ---- FLAGS (simc.f:872-912)
    const std::map<std::string, bool> fl = {
        {"doing_eep", f.doing_eep},         {"doing_kaon", f.doing_kaon},
        {"doing_pion", f.doing_pion},       {"doing_semi", f.doing_semi},
        {"doing_rho", f.doing_rho},         {"doing_hplus", f.doing_hplus},
        {"doing_semipi", f.doing_semipi},   {"doing_semika", f.doing_semika},
        {"doing_pizero", f.doing_pizero},   {"doing_delta", f.doing_delta},
        {"doing_phsp", f.doing_phsp},       {"doing_hyd_elast", f.doing_hyd_elast},
        {"doing_deuterium", f.doing_deuterium}, {"doing_heavy", f.doing_heavy},
        {"doing_hydpi", f.doing_hydpi},     {"doing_deutpi", f.doing_deutpi},
        {"doing_hepi", f.doing_hepi},       {"doing_hydkaon", f.doing_hydkaon},
        {"doing_deutkaon", f.doing_deutkaon}, {"doing_hekaon", f.doing_hekaon},
        {"doing_hydsemi", f.doing_hydsemi}, {"doing_deutsemi", f.doing_deutsemi},
        {"do_fermi", f.do_fermi},           {"doing_hydrho", f.doing_hydrho},
        {"doing_deutrho", f.doing_deutrho}, {"doing_herho", f.doing_herho},
        {"mc_smear", sw.mc_smear},          {"using_Eloss", sw.using_Eloss},
        {"using_Coulomb", sw.using_Coulomb}, {"correct_Eloss", sw.correct_Eloss},
        {"correct_raster", s.params.correct_raster}, {"doing_decay", s.params.doing_decay},
        {"using_E_arm_monteca", sw.using_E_arm_montecarlo},   // a19 truncation in the .hist
        {"using_P_arm_monteca", sw.using_P_arm_montecarlo},
        {"use_benhar_sf", s.params.use_benhar_sf}};
    for (const auto& [k, v] : fl) C.item("flags", k, L(v));
    C.item("flags", "which_pion", std::to_string(s.params.which_pion));
    C.item("flags", "which_kaon", std::to_string(s.params.which_kaon));
    C.item("flags", "pizero_ngamma", std::to_string(s.params.pizero_ngamma));
    C.item("flags", "electron_arm", std::to_string(s.params.electron_arm));
    C.item("flags", "hadron_arm", std::to_string(s.params.hadron_arm));
    // deForest_flag and use_first_cer are INTEGERS printed with an L edit descriptor (simc.f:896, 906);
    // only the value 0 -> 'F' is unambiguous.
    if (s.params.deForest_flag == 0) C.item("flags", "deForest_flag", "F");
    else C.defer("flags: deForest_flag (integer printed with L2, value /= 0)");
    const bool shms = s.params.electron_arm == 5 || s.params.hadron_arm == 5 || s.params.electron_arm == 6 ||
                      s.params.hadron_arm == 6;
    if (shms) {
        if (s.params.use_first_cer == 0) C.item("flags", "use_first_cer", "F");
        else C.defer("flags: use_first_cer (integer printed with L2, value /= 0)");
    }
    C.item("flags", "ctau", fmtF(s.params.ctau, 10, 3));

    // ---- COUNTERS
    C.kv("counters", "Ngen (request)", 0, std::to_string(s.params.ngen));
    C.defer("counters: Ntried, Ncontribute, Nco_no_rad_prot, %no_rad_prot, wtcontr [event loop]");

    // ---- RADIATIVE CORRECTIONS (simc.f:918-944)
    C.item("radiative", "using_rad", L(sw.using_rad));
    if (sw.using_rad) {
        C.item("radiative", "use_expon", L(s.params.use_expon != 0));   // integer printed with L3
        C.item("radiative", "include_hard", L(sw.include_hard));
        C.item("radiative", "calc_spence", L(sw.calculate_spence));
        C.item("radiative", "use_offshell_r", L(s.params.use_offshell_rad));
        C.item("radiative", "rad_flag", std::to_string(sw.rad_flag));
        C.item("radiative", "extrad_flag", std::to_string(sw.extrad_flag));
        C.item("radiative", "one_tail", std::to_string(s.params.one_tail));
        C.item("radiative", "intcor_mode", std::to_string(s.params.intcor_mode));
        C.item("radiative", "dE_edge_test", fmtF(s.params.dE_edge_test, 11, 3));
        C.item("radiative", "Egamma_max", fmtF(s.rad.Egamma_tot_max, 11, 3), LO);
        C.kv("radiative", "etatzai", 0, fmtF(s.rad.etatzai, 11, 3));
        C.defer("radiative central values: hardcorfac, frac, lambda, bt, c_int, c_ext, c, g_int, g_ext, g "
                "[radc_init_ev at the central event, P2.3]");
    }

    // ---- MISCELLANEOUS (simc.f:948-973)
    C.kv("misc", "charge", 0, fmtE(s.params.EXPER.charge, 16, 6));
    C.kv("misc", "targetfac", 0, fmtE(s.targetfac, 16, 6));
    C.kv("misc", "luminosity", 0, fmtE(s.luminosity, 16, 6), Dep::None, 1);
    C.kv("misc", "luminosity", 0, fmtE(s.luminosity * ((constants::hbarc / 100000.) * (constants::hbarc / 100000.)), 16, 6), Dep::None, 2);
    if (s.params.random_seed != 0) C.kv("misc", "Random Seed", 0, std::to_string(s.params.random_seed));
    C.defer("misc: CENTRAL.sigcc [calculate_central], AVERAGE.sigcc, genvol, normfac [event loop / generation]");
    C.defer("recon summary: ave.error / resolution [event loop]");

    // ---- Input Spectrometer Limits (simc.f:990-1002): SPedge after slop widening, f15.4
    auto spe = [&](const char* arm, const ArmCuts& c) {
        const std::string a = arm;
        C.kv("spedge", "SPedge." + a + ".delta.min/max", 0, fmtF(c.delta.min, 15, 4));
        C.kv("spedge", "SPedge." + a + ".delta.min/max", 1, fmtF(c.delta.max, 15, 4));
        C.kv("spedge", "SPedge." + a + ".yptar.min/max", 0, fmtF(c.yptar.min, 15, 4));
        C.kv("spedge", "SPedge." + a + ".yptar.min/max", 1, fmtF(c.yptar.max, 15, 4));
        C.kv("spedge", "SPedge." + a + ".xptar.min/max", 0, fmtF(c.xptar.min, 15, 4));
        C.kv("spedge", "SPedge." + a + ".xptar.min/max", 1, fmtF(c.xptar.max, 15, 4));
    };
    spe("e", s.SPedge.e);
    spe("p", s.SPedge.p);
    // The same numbers through the run-summary JSON echo (SummaryJson.cpp, P2.2)
    {
        const Json sj = BuildSummaryJson(s, RunSummary{});
        for (const char* arm : {"e", "p"})
            for (const char* v : {"delta", "yptar", "xptar"})
                for (int k = 0; k < 2; ++k)
                    C.kv("spedge", std::string("SPedge.") + arm + "." + v + ".min/max", k,
                         fmtF(sj["spedge"][arm][v][k].get<double>(), 15, 4));
    }

    // ---- Limiting VERTEX values (simc.f:1007-1040), format 9917 f12.3 (used columns only)
    auto pair = [&](const std::string& sec, const std::string& label, double a, double b, Dep d) {
        C.row(sec, label, 0, fmtF(a, 12, 3), d);
        C.row(sec, label, 1, fmtF(b, 12, 3), d);
    };
    const GenLimits& g = s.gen;
    pair("vertex", "E arm delta", g.e.delta.min, g.e.delta.max, LO);
    pair("vertex", "E arm yptar", g.e.yptar.min * 1000., g.e.yptar.max * 1000., LO);
    pair("vertex", "E arm xptar", g.e.xptar.min * 1000., g.e.xptar.max * 1000., LO);
    pair("vertex", "P arm delta", g.p.delta.min, g.p.delta.max, LO);
    pair("vertex", "P arm yptar", g.p.yptar.min * 1000., g.p.yptar.max * 1000., LO);
    pair("vertex", "P arm xptar", g.p.xptar.min * 1000., g.p.xptar.max * 1000., LO);
    pair("vertex", "sumEgen", g.sumEgen.min, g.sumEgen.max, LO);
    pair("vertex", "Trec", s.VERTEXedge.Trec.min, s.VERTEXedge.Trec.max, LO);
    pair("vertex", "Em", s.VERTEXedge.Em.min, s.VERTEXedge.Em.max, LO);
    pair("vertex", "Pm", s.VERTEXedge.Pm.min, s.VERTEXedge.Pm.max, LO);

    // ---- Limiting ORIGINAL values (simc.f:1045-1066)
    const EdgeTrue& e = s.edge;
    pair("original", "E arm E", e.e.E.min, e.e.E.max, LO);
    pair("original", "E arm yptar", e.e.yptar.min * 1000., e.e.yptar.max * 1000., LO);
    pair("original", "E arm xptar", e.e.xptar.min * 1000., e.e.xptar.max * 1000., LO);
    pair("original", "P arm E", e.p.E.min, e.p.E.max, LO);
    pair("original", "P arm yptar", e.p.yptar.min * 1000., e.p.yptar.max * 1000., LO);
    pair("original", "P arm xptar", e.p.xptar.min * 1000., e.p.xptar.max * 1000., LO);
    pair("original", "Em", std::max(-999999.999e0, e.Em.min), std::min(999999.999e0, e.Em.max), LO);

    // ---- Limiting SPECTROMETER values (simc.f:1070-1093)
    pair("spectrometer", "E arm delta", s.SPedge.e.delta.min, s.SPedge.e.delta.max, Dep::None);
    pair("spectrometer", "E arm yptar", s.SPedge.e.yptar.min * 1000., s.SPedge.e.yptar.max * 1000., Dep::None);
    pair("spectrometer", "E arm xptar", s.SPedge.e.xptar.min * 1000., s.SPedge.e.xptar.max * 1000., Dep::None);
    pair("spectrometer", "P arm delta", s.SPedge.p.delta.min, s.SPedge.p.delta.max, Dep::None);
    pair("spectrometer", "P arm yptar", s.SPedge.p.yptar.min * 1000., s.SPedge.p.yptar.max * 1000., Dep::None);
    pair("spectrometer", "P arm xptar", s.SPedge.p.xptar.min * 1000., s.SPedge.p.xptar.max * 1000., Dep::None);

    // ---- Limiting RADIATION values (simc.f:1101-1113)
    if (sw.using_rad) {
        pair("radiation", "Egamma(1)", 0., s.rad.Egamma1_max, LO);
        pair("radiation", "Egamma(2)", 0., s.rad.Egamma2_max, LO);
        pair("radiation", "Egamma(3)", 0., s.rad.Egamma3_max, LO);
        pair("radiation", "Egamma_total", 0., s.rad.Egamma_tot_max, LO);
    }
    C.defer("all 'found' lo/hi columns of the limit tables and the slop min/max columns [event loop]");

    // ---- slops (simc.f:1116-1136), format 9918 f10.3 (used column)
    const Slops& sl = s.slop;
    C.row("slop", "E arm delta", 0, fmtF(sl.MC.e.delta.used, 10, 3));
    C.row("slop", "E arm yptar", 0, fmtF(sl.MC.e.yptar.used * 1000., 10, 3));
    C.row("slop", "E arm xptar", 0, fmtF(sl.MC.e.xptar.used * 1000., 10, 3));
    C.row("slop", "P arm delta", 0, fmtF(sl.MC.p.delta.used, 10, 3));
    C.row("slop", "P arm yptar", 0, fmtF(sl.MC.p.yptar.used * 1000., 10, 3));
    C.row("slop", "P arm xptar", 0, fmtF(sl.MC.p.xptar.used * 1000., 10, 3));
    C.row("slop", "Em", 0, fmtF(sl.total.Em.used, 10, 3), f.doing_heavy ? LO : Dep::None);

    std::printf("\nDEFERRED (%d):\n", C.res.deferred);
    for (const auto& d : C.res.deferred_items) std::printf("  - %s\n", d.c_str());
    std::printf("\nSUMMARY %s: passed %d, failed %d, deferred %d\n", histfile.c_str(), C.res.pass, C.res.fail,
                C.res.deferred);
    return C.res.fail == 0 ? 0 : 1;
}
