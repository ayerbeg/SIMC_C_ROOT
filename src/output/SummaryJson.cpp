// src/output/SummaryJson.cpp
// Fortran: simc.f:report (l.644-1142) for the .hist layout, simc.f l.462-532 for the .geni.
#include "simc/output/SummaryJson.h"
#include "simc/core/Units.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <regex>
#include <sstream>

namespace simc {

namespace {

constexpr double kHbarc = 197.327053;   // constants.inc hbarc

std::vector<std::string> splitLines(const std::string& t) {
    std::vector<std::string> v;
    std::istringstream is(t);
    for (std::string l; std::getline(is, l);) {
        if (!l.empty() && l.back() == '\r') l.pop_back();
        v.push_back(l);
    }
    return v;
}

std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t"), b = s.find_last_not_of(" \t");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}

std::string collapse(const std::string& s) {   // normalise internal white space
    std::string out;
    bool sp = false;
    for (char c : trim(s)) {
        if (c == ' ' || c == '\t') { sp = true; continue; }
        if (sp && !out.empty()) out += ' ';
        sp = false;
        out += c;
    }
    return out;
}

// A Fortran number token; "****" (field overflow) and NaN become null.
Json numTok(const std::string& t) {
    if (t.find('*') != std::string::npos) return nullptr;
    char* end = nullptr;
    std::string s = t;
    for (auto& c : s) if (c == 'D' || c == 'd') c = 'E';
    double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || std::isnan(v)) return nullptr;
    return v;
}

Json numOrInt(const std::string& t) {
    Json v = numTok(t);
    if (v.is_number() && t.find_first_of(".eEdD") == std::string::npos) return static_cast<long>(v.get<double>());
    return v;
}

std::vector<std::string> tokens(const std::string& s) {
    std::istringstream is(s);
    std::vector<std::string> v;
    for (std::string t; is >> t;) v.push_back(t);
    return v;
}

Json arr(const std::vector<std::string>& toks, size_t from, size_t n) {
    Json a = Json::array();
    for (size_t i = from; i < from + n && i < toks.size(); ++i) a.push_back(numTok(toks[i]));
    return a;
}

}  // namespace

// ------------------------------------------------------------------ stop counter keys
const std::vector<std::string>& StopCounterKeys(const std::string& spec) {
    static const std::map<std::string, std::vector<std::string>> keys = {
        // hms/mc_hms.f, struct_hms.inc (coll is not printed in the .geni)
        {"hms", {"trials", "slit_hor", "slit_vert", "slit_oct", "Q1_in", "Q1_mid", "Q1_out", "Q2_in", "Q2_mid",
                 "Q2_out", "Q3_in", "Q3_mid", "Q3_out", "D1_in", "D1_out", "hut", "dc1", "dc2", "scin", "cal",
                 "successes", "coll"}},
        // shms/mc_shms.f, struct_shms.inc (D1_men, targ_*, cer1/2, BP_*, coll are not printed in the .geni)
        // shms/mc_shms.f, struct_shms.inc. Names as in the P1.5 reader (validation/compare/RunSummary.cpp):
        // hb_*/q*_* lower case, D1_* upper case. q2_mid_dup = the value the .geni prints in the Q3 "mid"
        // column (FORTRAN-QUIRK simc.f:519-520: it is shmsSTOP_q2_mid). D1_men, targ_*, cer1/2, BP_*, coll
        // are not printed in the .geni.
        {"shms", {"trials", "hb_in", "hb_men", "hb_mex", "hb_out", "slit_hor", "slit_vert", "slit_oct", "q1_in",
                  "q1_men", "q1_mid", "q1_mex", "q1_out", "q2_in", "q2_men", "q2_mid", "q2_mex", "q2_out", "q3_in",
                  "q3_men", "q2_mid_dup", "q3_mex", "q3_out", "D1_in", "D1_flr", "D1_mid1", "D1_mid2", "D1_mid3",
                  "D1_mid4", "D1_mid5", "D1_mid6", "D1_mid7", "D1_mex", "D1_out", "hut", "dc1", "dc2", "s1", "s2",
                  "s3", "cal", "successes", "q3_mid", "D1_men", "targ_hor", "targ_vert", "targ_oct", "cer1", "cer2",
                  "BP_in", "BP_out", "coll"}},
        {"sos", {"trials", "slit_hor", "slit_vert", "slit_oct", "quad_in", "quad_mid", "quad_out", "bm01_in",
                 "bm01_out", "bm02_in", "bm02_out", "exit", "hut", "dc1", "dc2", "scin", "cal", "successes"}},
        {"hrsr", {"trials", "slit_hor", "slit_vert", "Q1_in", "Q1_mid", "Q1_out", "Q2_in", "Q2_mid", "Q2_out",
                  "D1_in", "D1_out", "Q3_in", "Q3_mid", "Q3_out", "hut", "dc1", "dc2", "s1", "s2", "cal",
                  "successes"}},
        {"hrsl", {"trials", "slit_hor", "slit_vert", "Q1_in", "Q1_mid", "Q1_out", "Q2_in", "Q2_mid", "Q2_out",
                  "D1_in", "D1_out", "Q3_in", "Q3_mid", "Q3_out", "hut", "dc1", "dc2", "s1", "s2", "cal",
                  "successes"}},
        {"calo", {"trials", "slit_hor", "slit_vert", "successes"}},
    };
    static const std::vector<std::string> none;
    auto it = keys.find(spec);
    return it == keys.end() ? none : it->second;
}

std::string SpectrometerName(int arm) {
    switch (arm) {
        case 1: return "hms";
        case 2: return "sos";
        case 3: return "hrsr";
        case 4: return "hrsl";
        case 5:
        case 6: return "shms";
        case 7:
        case 8: return "calo";
        default: return "none";
    }
}

// ------------------------------------------------------------------ C++ side
Json BuildSummaryJson(const RunSetup& s, const RunSummary& r) {
    const SimcParameters& p = s.params;
    Json j;
    j["schema"] = kSummarySchemaName;
    j["schema_version"] = kSummarySchemaVersion;
    j["run"] = {{"tag", r.run_name.empty() ? s.run_name : r.run_name},
                {"producer", "simc++"},
                {"random_seed", r.random_seed},
                {"reaction", r.reaction_title.empty() ? ReactionTitle(s.flags, p) : r.reaction_title}};

    j["counters"] = {{"ngen", r.ngen},
                     {"ntried", r.ntried},
                     {"ncontribute", r.ncontribute},
                     {"nco_no_rad_prot", r.ncontribute_no_rad_proton}};
    const double nev = r.nevent > 0 ? static_cast<double>(r.nevent) : 1.0;
    // simc.f:report — wtcontr = wtcontribute/nevent (l.920), AVERAGE.sigcc = sum_sigcc/nevent
    j["normalization"] = {{"wtcontr", r.wtcontribute / nev},
                          {"sigcc_central", r.central_sigcc},
                          {"sigcc_average", r.sum_sigcc / nev},
                          {"charge_mC", r.charge},
                          {"targetfac", r.targetfac},
                          {"luminosity_ub", r.luminosity},
                          {"luminosity_GeV2", r.luminosity * std::pow(kHbarc / 100000., 2)},
                          {"genvol", r.genvol},
                          {"normfac", r.normfac}};

    // Flags echo: the logical values of the .hist FLAGS / RADIATIVE blocks, keyed exactly as printed
    // (including the a19 truncation of using_E/P_arm_monteca and a14 of use_offshell_r).
    const ReactionFlags& f = s.flags;
    j["flags"] = {{"doing_eep", f.doing_eep},       {"doing_kaon", f.doing_kaon},
                  {"doing_pion", f.doing_pion},     {"doing_semi", f.doing_semi},
                  {"doing_rho", f.doing_rho},       {"doing_hplus", f.doing_hplus},
                  {"doing_semipi", f.doing_semipi}, {"doing_semika", f.doing_semika},
                  {"doing_pizero", f.doing_pizero}, {"doing_delta", f.doing_delta},
                  {"doing_phsp", f.doing_phsp},
                  {"doing_hyd_elast", f.doing_hyd_elast}, {"doing_deuterium", f.doing_deuterium},
                  {"doing_heavy", f.doing_heavy},   {"doing_hydpi", f.doing_hydpi},
                  {"doing_deutpi", f.doing_deutpi}, {"doing_hepi", f.doing_hepi},
                  {"doing_hydkaon", f.doing_hydkaon}, {"doing_deutkaon", f.doing_deutkaon},
                  {"doing_hekaon", f.doing_hekaon}, {"doing_hydsemi", f.doing_hydsemi},
                  {"doing_deutsemi", f.doing_deutsemi}, {"do_fermi", f.do_fermi},
                  {"doing_hydrho", f.doing_hydrho}, {"doing_deutrho", f.doing_deutrho},
                  {"doing_herho", f.doing_herho},   {"mc_smear", p.mc_smear},
                  {"using_Eloss", p.using_Eloss},
                  {"using_Coulomb", p.using_Coulomb && std::lround(p.targ.Z) != 1},   // dbase.f:547
                  {"deForest_flag", p.deForest_flag != 0},                           // printed with L edit
                  {"correct_Eloss", p.correct_Eloss && p.using_Eloss},               // dbase.f:542
                  {"correct_raster", p.correct_raster}, {"doing_decay", p.doing_decay},
                  {"using_E_arm_monteca", p.spect_mode != 1 && p.spect_mode != -1},
                  {"using_P_arm_monteca", p.spect_mode != 1 && p.spect_mode != -2},
                  {"use_benhar_sf", p.use_benhar_sf},
                  {"using_rad", p.using_rad}};
    if (p.using_rad) {
        j["flags"]["use_expon"] = p.use_expon != 0;
        j["flags"]["use_offshell_r"] = p.use_offshell_rad;
    }
    j["settings"] = {{"which_pion", p.which_pion}, {"which_kaon", p.which_kaon},
                     {"pizero_ngamma", p.pizero_ngamma}, {"electron_arm", p.electron_arm},
                     {"hadron_arm", p.hadron_arm}, {"ctau", p.ctau}};

    // Kinematics echo (simc.f:report l.791-815), .hist units
    auto arm = [](const SimcParameters::Arm& a) {
        return Json{{"angle_deg", std::abs(a.theta)}, {"momentum", a.P},
                    {"offset_x", a.offset.x},       {"offset_y", a.offset.y},
                    {"offset_z", a.offset.z},       {"offset_xptar_mr", a.offset.xptar},
                    {"offset_yptar_mr", a.offset.yptar}};
    };
    j["kinematics"] = {{"Ebeam", p.Ebeam},
                       {"dE_over_E", p.dEbeam / units::kPercent},
                       {"x_width", p.gen.xwid},
                       {"y_width", p.gen.ywid},
                       {"fr_pattern", p.targ.fr_pattern},
                       {"fr1", p.targ.fr1},
                       {"fr2", p.targ.fr2},
                       {"e", arm(p.spec.e)},
                       {"p", arm(p.spec.p)}};

    // SPedge echo, as the Fortran .hist prints it (simc.f:990-1002): after the dbase.f conversions and the
    // slop%MC widening of limits_init (init.f:201-212), i.e. RunSetup::SPedge from BuildRunSetup (P2.2).
    // A setup filled only by RunSetup::FromParameters falls back to the converted input (and says so).
    j["notes"] = Json::array();
    if (s.built) {
        auto edge = [](const ArmCuts& a) {
            return Json{{"delta", {a.delta.min, a.delta.max}},
                        {"yptar", {a.yptar.min, a.yptar.max}},
                        {"xptar", {a.xptar.min, a.xptar.max}}};
        };
        j["spedge"] = {{"e", edge(s.SPedge.e)}, {"p", edge(s.SPedge.p)}};
    } else {
        auto edge = [](const SimcParameters::ArmEdge& a) {
            const double dmin = a.delta.min <= -100.0 ? -99.99 : a.delta.min;   // dbase.f:479-480
            return Json{{"delta", {dmin, a.delta.max}},
                        {"yptar", {a.yptar.min / units::kMradPerRad, a.yptar.max / units::kMradPerRad}},
                        {"xptar", {a.xptar.min / units::kMradPerRad, a.xptar.max / units::kMradPerRad}}};
        };
        j["spedge"] = {{"e", edge(p.SPedge.e)}, {"p", edge(p.SPedge.p)}};
        if (!r.extra.contains("spedge"))
            j["notes"].push_back("spedge is the converted input, not widened by slop%MC (setup not built by BuildRunSetup)");
    }

    // Stop counters: stops.<spectrometer>.<key> (flat integers), roles in stop_roles.
    j["stops"] = Json::object();
    j["stop_roles"] = Json::object();
    for (const auto& sc : r.stops) {
        Json c = Json::object();
        for (const auto& [k, v] : sc.counts) c[k] = v;
        // mirror what the Fortran .geni prints in the SHMS Q3 "mid" column (simc.f:519-520)
        if (sc.spectrometer == "shms" && sc.counts.count("q2_mid") && !sc.counts.count("q2_mid_dup"))
            c["q2_mid_dup"] = sc.counts.at("q2_mid");
        j["stops"][sc.spectrometer] = c;
        if (!sc.role.empty()) j["stop_roles"][sc.spectrometer] = sc.role;
    }

    // Tree written by RootTreeSink (default file name, relative to the summary JSON).
    j["tree"] = {{"file", s.run_name + ".root"}, {"name", "h10"}, {"weight_branch", "Weight"}};

    for (auto it = r.extra.begin(); it != r.extra.end(); ++it) j[it.key()] = it.value();
    return j;
}

// ------------------------------------------------------------------ Fortran .geni
Json ParseFortranGeni(const std::string& geni, std::vector<std::string>* warnings) {
    Json out = Json::object();
    std::string cur;
    auto put = [&](const std::vector<std::string>& keys, const std::vector<std::string>& toks, size_t first) {
        for (size_t i = 0; i < keys.size(); ++i) {
            if (first + i >= toks.size()) break;
            if (keys[i].empty()) continue;
            out[cur][keys[i]] = std::strtol(toks[first + i].c_str(), nullptr, 10);
        }
    };
    for (const auto& line : splitLines(geni)) {
        const std::string l = collapse(line);
        if (l.empty()) continue;
        auto toks = tokens(l);
        // Number of leading label tokens = tokens that are not integers.
        size_t nlab = 0;
        // label tokens = everything before the first pure integer (^-?[0-9]+$; "1-2" is part of a label)
        auto isInt = [](const std::string& t) {
            size_t i = (t.size() > 1 && t[0] == '-') ? 1 : 0;
            return i < t.size() && t.find_first_not_of("0123456789", i) == std::string::npos;
        };
        while (nlab < toks.size() && !isInt(toks[nlab])) ++nlab;
        const std::string label = l.substr(0, l.find(toks.size() > nlab ? toks[nlab] : std::string("\x01")));
        auto starts = [&](const char* s) { return l.rfind(s, 0) == 0; };
        if (starts("HMS Trials:")) { cur = "hms"; put({"trials"}, toks, nlab); continue; }
        if (starts("SHMS Trials:")) { cur = "shms"; put({"trials"}, toks, nlab); continue; }
        if (starts("SOS Trials:")) { cur = "sos"; put({"trials"}, toks, nlab); continue; }
        if (starts("HRSr Trials:")) { cur = "hrsr"; put({"trials"}, toks, nlab); continue; }
        if (starts("HRSl Trials:")) { cur = "hrsl"; put({"trials"}, toks, nlab); continue; }
        if (starts("Calo Trials:")) { cur = "calo"; put({"trials"}, toks, nlab); continue; }
        if (cur.empty()) { if (warnings) warnings->push_back("geni: line before any spectrometer: " + l); continue; }
        std::vector<std::string> keys;
        if (cur == "hms") {
            if (starts("Slit hor/vert/corners")) keys = {"slit_hor", "slit_vert", "slit_oct"};
            else if (starts("Q1 entrance")) keys = {"Q1_in", "Q1_mid", "Q1_out"};
            else if (starts("Q2 entrance")) keys = {"Q2_in", "Q2_mid", "Q2_out"};
            else if (starts("Q3 entrance")) keys = {"Q3_in", "Q3_mid", "Q3_out"};
            else if (starts("Dipole entrance/exit")) keys = {"D1_in", "D1_out"};
            else if (starts("Events reaching hut")) keys = {"hut"};
            else if (starts("DC1, DC2, Scin, Cal")) keys = {"dc1", "dc2", "scin", "cal"};
            else if (starts("Successes")) keys = {"successes"};
        } else if (cur == "shms") {
            if (starts("HB phys")) keys = {"hb_in", "hb_men", "hb_mex", "hb_out"};
            else if (starts("Slit hor/vert/corners")) keys = {"slit_hor", "slit_vert", "slit_oct"};
            else if (starts("Q1 phys")) keys = {"q1_in", "q1_men", "q1_mid", "q1_mex", "q1_out"};
            else if (starts("Q2 phys")) keys = {"q2_in", "q2_men", "q2_mid", "q2_mex", "q2_out"};
            // FORTRAN-QUIRK simc.f:519-520: the Q3 "mid" column prints shmsSTOP_q2_mid -> "q2_mid_dup".
            else if (starts("Q3 phys")) keys = {"q3_in", "q3_men", "q2_mid_dup", "q3_mex", "q3_out"};
            else if (starts("D1 entrance/flare/mid 1-2")) keys = {"D1_in", "D1_flr", "D1_mid1", "D1_mid2"};
            else if (starts("D1 mid 3-5")) keys = {"D1_mid3", "D1_mid4", "D1_mid5"};
            else if (starts("D1 mid 6-7")) keys = {"D1_mid6", "D1_mid7", "D1_mex", "D1_out"};
            else if (starts("Events reaching hut")) keys = {"hut"};
            else if (starts("DC1, DC2, Scin, Cal")) keys = {"dc1", "dc2"};
            else if (starts("S1, S2, S3, Cal")) keys = {"s1", "s2", "s3", "cal"};
            else if (starts("Successes")) keys = {"successes"};
        } else if (cur == "sos") {
            // FORTRAN-QUIRK: simc.f prints sSTOP_slit_vert before sSTOP_slit_hor under "hor/vert".
            if (starts("Slit hor/vert/corners")) keys = {"slit_vert", "slit_hor", "slit_oct"};
            else if (starts("Quad entrance")) keys = {"quad_in", "quad_mid", "quad_out"};
            else if (starts("D1 entrance/exit")) keys = {"bm01_in", "bm01_out"};
            else if (starts("D2 entrance/exit")) keys = {"bm02_in", "bm02_out"};
            else if (starts("Vacuum exit")) keys = {"exit"};
            else if (starts("Events reaching hut")) keys = {"hut"};
            else if (starts("DC1, DC2, Scin, Cal")) keys = {"dc1", "dc2", "scin", "cal"};
            else if (starts("Successes")) keys = {"successes"};
        } else if (cur == "hrsr" || cur == "hrsl") {
            // FORTRAN-QUIRK: "Slit hor/vert" prints slit_vert first.
            if (starts("Slit hor/vert")) keys = {"slit_vert", "slit_hor"};
            else if (starts("Q1 entrance")) keys = {"Q1_in", "Q1_mid", "Q1_out"};
            else if (starts("Q2 entrance")) keys = {"Q2_in", "Q2_mid", "Q2_out"};
            else if (starts("Dipole entrance/exit")) keys = {"D1_in", "D1_out"};
            else if (starts("Q3 entrance")) keys = {"Q3_in", "Q3_mid", "Q3_out"};
            else if (starts("Events reaching hut")) keys = {"hut"};
            else if (starts("VDC1, VDC2")) keys = {"dc1", "dc2"};
            else if (starts("S1, S2, Cal")) keys = {"s1", "s2", "cal"};
        } else if (cur == "calo") {
            if (starts("Extent hor/vert")) keys = {"slit_hor", "slit_vert"};
            else if (starts("Successes")) keys = {"successes"};
        }
        if (keys.empty()) {
            if (warnings) warnings->push_back("geni: unrecognised line for " + cur + ": " + l);
            continue;
        }
        put(keys, toks, nlab);
        (void)label;
    }
    return out;
}

// ------------------------------------------------------------------ Fortran .hist
Json ParseFortranHist(const std::string& hist, const std::string& geni, const std::string& run_name,
                      std::vector<std::string>* warnings) {
    auto warn = [&](const std::string& w) { if (warnings) warnings->push_back(w); };
    Json j;
    j["schema"] = kSummarySchemaName;
    j["schema_version"] = kSummarySchemaVersion;
    j["run"] = {{"tag", run_name}, {"producer", "hist2json (Fortran .hist/.geni)"}};
    j["notes"] = Json::array();

    const auto lines = splitLines(hist);
    std::string section;
    std::string limits_table;
    std::string recon_arm;     // "e" / "p" in the RECON SUMMARY block
    std::string slop_group;    // "MC" / "total" in the slop block
    Json flags = Json::object();
    std::smatch m;

    // "label = v1 [v2] [unit]" with label possibly containing spaces/punctuation
    const std::regex re_assign(R"(^\s*(.+?)\s*=\s*(.*)$)");
    const std::regex re_flag(R"(([A-Za-z_][A-Za-z0-9_]*)=\s*([TF]|-?\d+)\b)");
    const std::regex re_table4(R"(^\s*(.+?)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$)");

    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string& raw = lines[i];
        const std::string l = collapse(raw);
        if (l.empty()) continue;

        // ---- section headers (simc.f:report)
        if (l.rfind("BEGIN Time:", 0) == 0) { j["run"]["begin_time"] = trim(l.substr(11)); continue; }
        if (l.rfind("END Time:", 0) == 0) { j["run"]["end_time"] = trim(l.substr(9)); continue; }
        if (l == "KINEMATICS:") { section = "kinematics"; continue; }
        if (l == "TARGET specs:") { section = "target"; continue; }
        if (l == "FLAGS:") { section = "flags"; continue; }
        if (l == "COUNTERS:") { section = "counters"; continue; }
        if (l == "RADIATIVE CORRECTIONS:") { section = "radiative"; continue; }
        if (l == "MISCELLANEOUS:") { section = "misc"; continue; }
        if (l.rfind("RECON SUMMARY:", 0) == 0) { section = "recon"; continue; }
        if (l == "Input Spectrometer Limits:") { section = "spedge"; continue; }
        if (l.rfind("Limiting VERTEX values", 0) == 0) { section = "limits"; limits_table = "vertex"; continue; }
        if (l.rfind("Limiting ORIGINAL values", 0) == 0) { section = "limits"; limits_table = "original"; continue; }
        if (l.rfind("Limiting SPECTROMETER values", 0) == 0) { section = "limits"; limits_table = "spectrometer"; continue; }
        if (l.rfind("Limiting RADIATION values", 0) == 0) { section = "limits"; limits_table = "radiation"; continue; }
        if (l.rfind("ACTUAL and LIMITING SLOP values", 0) == 0) { section = "slop"; continue; }
        if (l.rfind("INTEGRATED WEIGHTS", 0) == 0) { section = "weights"; continue; }
        if (l == "VALUES FOR \"CENTRAL\" EVENT:") { section = "central"; continue; }
        if (l.rfind("Central Values:", 0) == 0) { section = "radcentral"; continue; }
        if (l.find("****") != std::string::npos && section == "kinematics") {
            std::string t = l;
            auto a = t.find_first_not_of("*- ");
            auto b = t.find_last_not_of("*- ");
            t = (a == std::string::npos) ? "" : t.substr(a, b - a + 1);
            if (t.rfind("PHASE SPACE", 0) == 0) t = "PHASE SPACE";
            if (!j["run"].contains("reaction")) j["run"]["reaction"] = t; else j["run"]["reaction_detail"] = t;
            continue;
        }
        if (l.rfind("Note ", 0) == 0 || l.rfind("and may give", 0) == 0 || l.rfind("USED limits", 0) == 0 ||
            l.find("______") != std::string::npos || l.rfind("min max", 0) == 0 || l.rfind("__used__", 0) == 0 ||
            l.rfind("__ave__", 0) == 0 || l.rfind("____E arm____", 0) == 0)
            continue;

        auto toks = tokens(l);
        if (section == "kinematics" || section == "central") {
            if (!std::regex_match(l, m, re_assign)) { warn("hist: unparsed kinematics line: " + l); continue; }
            const std::string key = m[1];
            auto v = tokens(m[2]);
            auto two = [&](const std::string& k) {
                if (v.size() >= 2) {
                    j[section == "central" ? "central" : "kinematics"]["e"][k] = numTok(v[0]);
                    j[section == "central" ? "central" : "kinematics"]["p"][k] = numTok(v[1]);
                }
            };
            if (key == "Ebeam") j["kinematics"]["Ebeam"] = numTok(v.at(0));
            else if (key == "(dE/E)beam") j["kinematics"]["dE_over_E"] = numTok(v.at(0));
            else if (key == "x-width") j["kinematics"]["x_width"] = numTok(v.at(0));
            else if (key == "y-width") j["kinematics"]["y_width"] = numTok(v.at(0));
            else if (key == "fr_pattern") j["kinematics"]["fr_pattern"] = numOrInt(v.at(0));
            else if (key == "fr1") j["kinematics"]["fr1"] = numTok(v.at(0));
            else if (key == "fr2") j["kinematics"]["fr2"] = numTok(v.at(0));
            else if (key == "angle") two("angle_deg");
            else if (key == "momentum") two("momentum");
            else if (key == "x offset") two("offset_x");
            else if (key == "y offset") two("offset_y");
            else if (key == "z offset") two("offset_z");
            else if (key == "xptar offset") two("offset_xptar_mr");
            else if (key == "yptar offset") two("offset_yptar_mr");
            else if (key == "delta") two("delta");
            else if (key == "xptar") two("xptar_mr");
            else if (key == "yptar") two("yptar_mr");
            else if (key == "Q2") j["central"]["Q2_GeV2"] = numTok(v.at(0));
            else if (key == "q") j["central"]["q"] = numTok(v.at(0));
            else if (key == "nu") j["central"]["nu"] = numTok(v.at(0));
            else if (key == "recon Em") j["central"]["Em"] = numTok(v.at(0));
            else if (key == "recon Pm") j["central"]["Pm"] = numTok(v.at(0));
            else if (key == "recon W") j["central"]["W"] = numTok(v.at(0));
            else if (key == "recon MM") j["central"]["MM"] = numTok(v.at(0));
            else warn("hist: unknown kinematics key: " + key);
        } else if (section == "target") {
            // format 9911 (simc.f:report): 2(5x,a10,' = ',e12.6,1x,a5) after 2x -> fixed columns
            if (raw.find(" = ") != std::string::npos) {
                static const std::map<std::string, std::string> names = {
                    {"A", "A"}, {"Z", "Z"}, {"rho", "rho"}, {"thick", "thick_g_cm2"}, {"angle", "angle_deg"},
                    {"abundancy", "abundancy"}, {"X0", "X0"}, {"X0_cm", "X0_cm"}, {"length", "length"},
                    {"zoffset", "zoffset"}, {"xoffset", "xoffset"}, {"yoffset", "yoffset"},
                    {"Mtar_struc", "Mtar_struck"}, {"Mrec_struc", "Mrec_struck"}};
                for (int half = 0; half < 2; ++half) {
                    const size_t st = 2 + 36 * half;
                    if (raw.size() < st + 30) break;
                    const std::string k = trim(raw.substr(st + 5, 10));
                    const std::string val = trim(raw.substr(st + 18, 12));
                    const std::string unit = raw.size() > st + 31 ? trim(raw.substr(st + 31, 5)) : "";
                    std::string key;
                    if (k == "mass") key = unit == "amu" ? "mass_amu" : "mass_MeV";
                    else if (k == "Mrec") key = unit == "amu" ? "Mrec_amu" : "Mrec_MeV";
                    else if (names.count(k)) key = names.at(k);
                    else { warn("hist: unknown target key '" + k + "'"); continue; }
                    j["target"][key] = numTok(val);
                }
            } else if (toks.size() >= 4 && (toks[0] == "Coulomb" || toks[0].rfind("Eloss_", 0) == 0 ||
                                             toks[0].rfind("teff_", 0) == 0)) {
                j["target"]["losses"][toks[0]] = arr(toks, 1, 3);
            } else if (toks.size() >= 2 && toks[0].rfind("musc_", 0) == 0) {
                const std::string k = toks[0] == "musc_nsig_max" ? toks[0] : toks[0] + "_mr";
                j["target"][k] = numTok(toks[1]);
            } else {
                warn("hist: unparsed target line: " + l);
            }
        } else if (section == "flags") {
            if (toks.size() >= 2 && toks[0] == "ctau=") { j["settings"]["ctau"] = numTok(toks[1]); continue; }
            if (l.rfind("ctau=", 0) == 0) { j["settings"]["ctau"] = numTok(l.substr(5)); continue; }
            if (l.rfind("transparency=", 0) == 0) { j["settings"]["transparency"] = numTok(l.substr(13)); continue; }
            bool any = false;
            for (auto it = std::sregex_iterator(l.begin(), l.end(), re_flag); it != std::sregex_iterator(); ++it) {
                const std::string name = (*it)[1], val = (*it)[2];   // names as printed (a19)
                if (val == "T" || val == "F") flags[name] = (val == "T");
                else j["settings"][name] = std::atol(val.c_str());
                any = true;
            }
            if (!any) warn("hist: unparsed flags line: " + l);
        } else if (section == "counters") {
            if (!std::regex_match(l, m, re_assign)) { warn("hist: unparsed counters line: " + l); continue; }
            const std::string k = m[1];
            const Json v = numOrInt(trim(m[2]));
            if (k == "Ngen (request)") j["counters"]["ngen"] = v;
            else if (k == "Ntried") j["counters"]["ntried"] = v;
            else if (k == "Ncontribute") j["counters"]["ncontribute"] = v;
            else if (k == "Nco_no_rad_prot") j["counters"]["nco_no_rad_prot"] = v;
            else if (k == "-> %no_rad_prot") j["counters"]["pct_no_rad_proton"] = numTok(trim(m[2]));
            else warn("hist: unknown counter " + k);
        } else if (section == "weights") {
            if (l.rfind("MeV: wtcontr=", 0) == 0) j["normalization"]["wtcontr"] = numTok(toks.back());
            else warn("hist: unparsed weights line: " + l);
        } else if (section == "radiative") {
            bool any = false;
            static const std::regex rf(R"(([A-Za-z_][A-Za-z0-9_]*)=\s*(\S+))");
            for (auto it = std::sregex_iterator(l.begin(), l.end(), rf); it != std::sregex_iterator(); ++it) {
                const std::string name = (*it)[1], val = (*it)[2];   // names as printed (a14)
                if (val == "T" || val == "F") flags[name] = (val == "T");
                else j["radiative"][name] = numOrInt(val);
                any = true;
            }
            if (!any) warn("hist: unparsed radiative line: " + l);
        } else if (section == "radcentral") {
            if (!std::regex_match(l, m, re_assign)) { warn("hist: unparsed rad central line: " + l); continue; }
            std::string k = m[1];
            auto v = tokens(m[2]);
            auto pos = k.find('(');
            if (pos != std::string::npos) k = k.substr(0, pos);
            j["radiative"]["central"][k] = v.size() == 1 ? numTok(v[0]) : arr(v, 0, v.size());
        } else if (section == "misc") {
            if (l.rfind("Random Seed =", 0) == 0) { j["run"]["random_seed"] = numOrInt(toks.back()); continue; }
            if (l.rfind("Theory file:", 0) == 0) { j["run"]["theory_file"] = trim(l.substr(12)); continue; }
            if (l.rfind("Random State Save File:", 0) == 0) { j["run"]["random_state_file"] = trim(l.substr(23)); continue; }
            if (!std::regex_match(l, m, re_assign)) { warn("hist: unparsed misc line: " + l); continue; }
            const std::string k = m[1];
            auto v = tokens(m[2]);
            const std::string unit = v.size() > 1 ? v[1] : "";
            if (k == "CENTRAL.sigcc") j["normalization"]["sigcc_central"] = numTok(v.at(0));
            else if (k == "AVERAGE.sigcc") j["normalization"]["sigcc_average"] = numTok(v.at(0));
            else if (k == "charge") j["normalization"]["charge_mC"] = numTok(v.at(0));
            else if (k == "targetfac") j["normalization"]["targetfac"] = numTok(v.at(0));
            else if (k == "luminosity" && unit == "GeV^2") j["normalization"]["luminosity_GeV2"] = numTok(v.at(0));
            else if (k == "luminosity") j["normalization"]["luminosity_ub"] = numTok(v.at(0));
            else if (k == "genvol") j["normalization"]["genvol"] = numTok(v.at(0));
            else if (k == "normfac") j["normalization"]["normfac"] = numTok(v.at(0));
            else warn("hist: unknown misc key " + k);
        } else if (section == "recon") {
            // "Electron arm: delta = a r x10^-3" / "xptar = a r mr" ...
            std::string& arm = recon_arm;
            if (l.rfind("Electron arm:", 0) == 0) arm = "e";
            if (l.rfind("Hadron arm:", 0) == 0) arm = "p";
            if (!std::regex_match(l, m, re_assign)) { warn("hist: unparsed recon line: " + l); continue; }
            std::string k = tokens(m[1]).back();
            auto v = tokens(m[2]);
            if (v.size() >= 3)
                j["recon_summary"][arm][k] = {{"ave_error", numTok(v[0])}, {"resolution", numTok(v[1])}, {"units", v[2]}};
        } else if (section == "spedge") {
            static const std::regex re(R"(^SPedge\.([ep])\.(delta|yptar|xptar)\.min/max = (\S+) (\S+) ?(\S*)$)");
            if (std::regex_match(l, m, re)) j["spedge"][m[1].str()][m[2].str()] = {numTok(m[3]), numTok(m[4])};
            else warn("hist: unparsed SPedge line: " + l);
        } else if (section == "limits") {
            if (std::regex_match(l, m, re_table4)) {
                j["limits"][limits_table].push_back({{"label", collapse(m[1])},
                                                     {"used", {numTok(m[2]), numTok(m[3])}},
                                                     {"found", {numTok(m[4]), numTok(m[5])}},
                                                     {"units", m[6].str()}});
            } else {
                warn("hist: unparsed limits line: " + l);
            }
        } else if (section == "slop") {
            // "slop.MC E arm delta used min max unit" or continuation "E arm yptar ..."
            std::string& group = slop_group;
            std::string rest = l;
            if (rest.rfind("slop.MC ", 0) == 0) { group = "MC"; rest = rest.substr(8); }
            if (rest.rfind("slop.total ", 0) == 0) { group = "total"; rest = rest.substr(11); }
            static const std::regex re3(R"(^(.+?) (\S+) (\S+) (\S+) (\S+)$)");
            if (std::regex_match(rest, m, re3)) {
                j["slop"][group].push_back({{"label", collapse(m[1])},
                                            {"used", numTok(m[2])},
                                            {"min", numTok(m[3])},
                                            {"max", numTok(m[4])},
                                            {"units", m[5].str()}});
            } else {
                warn("hist: unparsed slop line: " + l);
            }
        } else {
            warn("hist: line outside known sections: " + l);
        }
    }
    j["flags"] = flags;
    j["stops"] = ParseFortranGeni(geni, warnings);
    j["stop_roles"] = Json::object();
    // role of each spectrometer from the arm codes of the FLAGS block
    if (j.contains("settings") && j["settings"].contains("electron_arm") && j["settings"].contains("hadron_arm")) {
        const std::string e = SpectrometerName(j["settings"]["electron_arm"].get<int>());
        const std::string h = SpectrometerName(j["settings"]["hadron_arm"].get<int>());
        for (auto it = j["stops"].begin(); it != j["stops"].end(); ++it)
            j["stop_roles"][it.key()] = (it.key() == e && it.key() == h) ? "both" : (it.key() == e ? "electron" : "hadron");
    }
    if (j["stops"].contains("shms"))
        j["notes"].push_back("shms q2_mid_dup: the Fortran .geni prints shmsSTOP_q2_mid in the Q3 mid column "
                             "(simc.f:519-520); q3_mid itself is never printed");
    j["notes"].push_back("spedge values are as printed: after widening by slop%MC (init.f:209-216)");
    return j;
}

// ------------------------------------------------------------------ text rendering
std::string FormatSummaryText(const Json& s) {
    std::ostringstream os;
    std::function<void(const Json&, const std::string&, int)> rec = [&](const Json& n, const std::string& key, int ind) {
        const std::string pad(ind * 2, ' ');
        if (n.is_object()) {
            if (!key.empty()) os << pad << key << ":\n";
            for (auto it = n.begin(); it != n.end(); ++it) rec(it.value(), it.key(), key.empty() ? ind : ind + 1);
        } else if (n.is_array() && !n.empty() && n[0].is_object()) {
            os << pad << key << ":\n";
            for (const auto& e : n) rec(e, "-", ind + 1);
        } else {
            os << pad << key << " = " << n.dump() << "\n";
        }
    };
    const Json run = s.value("run", Json::object());
    os << " SIMC run summary (" << run.value("producer", "?") << "): " << run.value("tag", "?") << "\n";
    rec(s, "", 0);
    return os.str();
}

// ------------------------------------------------------------------ schema validator (subset)
namespace {
bool typeOk(const Json& v, const std::string& t) {
    if (t == "object") return v.is_object();
    if (t == "array") return v.is_array();
    if (t == "string") return v.is_string();
    if (t == "boolean") return v.is_boolean();
    if (t == "integer") return v.is_number_integer() || (v.is_number_float() && std::floor(v.get<double>()) == v.get<double>());
    if (t == "number") return v.is_number();
    if (t == "null") return v.is_null();
    return false;
}
void validate(const Json& v, const Json& sc, const std::string& path, std::vector<std::string>& errs) {
    if (sc.contains("type")) {
        bool ok = false;
        if (sc["type"].is_array()) { for (const auto& t : sc["type"]) ok = ok || typeOk(v, t.get<std::string>()); }
        else ok = typeOk(v, sc["type"].get<std::string>());
        if (!ok) { errs.push_back(path + ": type mismatch (want " + sc["type"].dump() + ", have " + v.dump().substr(0, 60) + ")"); return; }
    }
    if (sc.contains("enum")) {
        bool ok = false;
        for (const auto& e : sc["enum"]) ok = ok || e == v;
        if (!ok) errs.push_back(path + ": value not in enum");
    }
    if (sc.contains("minimum") && v.is_number() && v.get<double>() < sc["minimum"].get<double>())
        errs.push_back(path + ": below minimum");
    if (v.is_object()) {
        if (sc.contains("required"))
            for (const auto& r : sc["required"])
                if (!v.contains(r.get<std::string>())) errs.push_back(path + ": missing required '" + r.get<std::string>() + "'");
        const Json props = sc.value("properties", Json::object());
        for (auto it = v.begin(); it != v.end(); ++it) {
            if (props.contains(it.key())) validate(it.value(), props[it.key()], path + "/" + it.key(), errs);
            else if (sc.contains("additionalProperties")) {
                const Json& ap = sc["additionalProperties"];
                if (ap.is_boolean() && !ap.get<bool>()) errs.push_back(path + ": unexpected property '" + it.key() + "'");
                else if (ap.is_object()) validate(it.value(), ap, path + "/" + it.key(), errs);
            }
        }
    }
    if (v.is_array() && sc.contains("items")) {
        for (size_t i = 0; i < v.size(); ++i) validate(v[i], sc["items"], path + "/" + std::to_string(i), errs);
        if (sc.contains("minItems") && v.size() < sc["minItems"].get<size_t>()) errs.push_back(path + ": too few items");
        if (sc.contains("maxItems") && v.size() > sc["maxItems"].get<size_t>()) errs.push_back(path + ": too many items");
    }
}
}  // namespace

std::vector<std::string> ValidateAgainstSchema(const Json& instance, const Json& schema) {
    std::vector<std::string> errs;
    validate(instance, schema, "", errs);
    return errs;
}

}  // namespace simc
