// RunSummary.cpp - readers for Fortran .hist/.geni and for the run summary JSON.
#include "RunSummary.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace simcval {

namespace {

std::string collapse(const std::string& s) {  // collapse whitespace runs, trim
    std::string o;
    bool sp = false;
    for (char c : s) {
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            sp = true;
        } else {
            if (sp && !o.empty()) o += ' ';
            o += c;
            sp = false;
        }
    }
    return o;
}

std::string lower(std::string s) {
    for (auto& c : s) c = std::tolower(static_cast<unsigned char>(c));
    return s;
}

// Parse a Fortran-printed real and return half a unit in its last digit.
Scalar parseFortranReal(const std::string& tok) {
    Scalar r;
    if (tok.find('*') != std::string::npos) return r;  // field overflow '*****'
    try {
        size_t pos = 0;
        r.value = std::stod(tok, &pos);
    } catch (...) {
        return r;
    }
    r.present = true;
    // digits after the decimal point in the mantissa, and the exponent
    std::string t = tok;
    int expo = 0;
    size_t e = t.find_first_of("EeDd");
    std::string mant = t;
    if (e != std::string::npos) {
        mant = t.substr(0, e);
        try { expo = std::stoi(t.substr(e + 1)); } catch (...) { expo = 0; }
    }
    size_t dot = mant.find('.');
    int ndec = (dot == std::string::npos) ? 0 : static_cast<int>(mant.size() - dot - 1);
    r.halfulp = 0.5 * std::pow(10.0, expo - ndec);
    return r;
}

// ---------------------------------------------------------------------------
// .geni stop-counter line table, mapped by PHYSICAL MEANING, i.e. by the
// variable actually written, not by the label (Fortran: simc.f:463-535).
//   Known label mistakes:
//   simc.f:477 SOS 'Slit hor/vert/corners' writes sSTOP_slit_vert, sSTOP_slit_hor, ...
//   simc.f:489,501 HRS 'Slit hor/vert' writes xSTOP_slit_vert, xSTOP_slit_hor
//   simc.f:519-520 SHMS 'Q3 ...mid...' writes shmsSTOP_q2_mid in the q3_mid slot
//                  (q3_mid is never printed) -> stored as "q2_mid_dup"
//   simc.f:526 SHMS 'DC1, DC2, Scin, Cal' writes only shmsSTOP_dc1, shmsSTOP_dc2
// ---------------------------------------------------------------------------
struct LineMap {
    const char* spec;
    const char* label;  // collapsed label text before the trailing integers
    std::vector<std::string> keys;
};

const std::vector<LineMap>& lineTable() {
    static const std::vector<LineMap> t = {
        // HMS  (simc.f:463-472)
        {"hms", "HMS Trials:", {"trials"}},
        {"hms", "Slit hor/vert/corners", {"slit_hor", "slit_vert", "slit_oct"}},
        {"hms", "Q1 entrance/mid/exit", {"Q1_in", "Q1_mid", "Q1_out"}},
        {"hms", "Q2 entrance/mid/exit", {"Q2_in", "Q2_mid", "Q2_out"}},
        {"hms", "Q3 entrance/mid/exit", {"Q3_in", "Q3_mid", "Q3_out"}},
        {"hms", "Dipole entrance/exit", {"D1_in", "D1_out"}},
        {"hms", "Events reaching hut", {"hut"}},
        {"hms", "DC1, DC2, Scin, Cal", {"dc1", "dc2", "scin", "cal"}},
        {"hms", "Successes", {"successes"}},
        // SOS  (simc.f:476-485)
        {"sos", "SOS Trials:", {"trials"}},
        {"sos", "Slit hor/vert/corners", {"slit_vert", "slit_hor", "slit_oct"}},
        {"sos", "Quad entrance/mid/exit", {"quad_in", "quad_mid", "quad_out"}},
        {"sos", "D1 entrance/exit", {"bm01_in", "bm01_out"}},
        {"sos", "D2 entrance/exit", {"bm02_in", "bm02_out"}},
        {"sos", "Vacuum exit", {"exit"}},
        {"sos", "Events reaching hut", {"hut"}},
        {"sos", "DC1, DC2, Scin, Cal", {"dc1", "dc2", "scin", "cal"}},
        {"sos", "Successes", {"successes"}},
        // HRSr (simc.f:488-497) and HRSl (simc.f:500-508); no Successes line
        {"hrsr", "HRSr Trials:", {"trials"}},
        {"hrsl", "HRSl Trials:", {"trials"}},
        {"hrs", "Slit hor/vert", {"slit_vert", "slit_hor"}},
        {"hrs", "Q1 entrance/mid/exit", {"Q1_in", "Q1_mid", "Q1_out"}},
        {"hrs", "Q2 entrance/mid/exit", {"Q2_in", "Q2_mid", "Q2_out"}},
        {"hrs", "Dipole entrance/exit", {"D1_in", "D1_out"}},
        {"hrs", "Q3 entrance/mid/exit", {"Q3_in", "Q3_mid", "Q3_out"}},
        {"hrs", "Events reaching hut", {"hut"}},
        {"hrs", "VDC1, VDC2", {"dc1", "dc2"}},
        {"hrs", "S1, S2, Cal", {"s1", "s2", "cal"}},
        // SHMS (simc.f:512-529)
        {"shms", "SHMS Trials:", {"trials"}},
        {"shms", "HB phys entrance/mag entr/mag exit/phys exit", {"hb_in", "hb_men", "hb_mex", "hb_out"}},
        {"shms", "Slit hor/vert/corners", {"slit_hor", "slit_vert", "slit_oct"}},
        {"shms", "Q1 phys entrance/mag entr/mid/mag exit/phys exit", {"q1_in", "q1_men", "q1_mid", "q1_mex", "q1_out"}},
        {"shms", "Q2 phys entrance/mag entr/mid/mag exit/phys exit", {"q2_in", "q2_men", "q2_mid", "q2_mex", "q2_out"}},
        {"shms", "Q3 phys entrance/mag entr/mid/mag exit/phys exit", {"q3_in", "q3_men", "q2_mid_dup", "q3_mex", "q3_out"}},
        {"shms", "D1 entrance/flare/mid 1-2", {"D1_in", "D1_flr", "D1_mid1", "D1_mid2"}},
        {"shms", "D1 mid 3-5", {"D1_mid3", "D1_mid4", "D1_mid5"}},
        {"shms", "D1 mid 6-7/mag exit/phys exit", {"D1_mid6", "D1_mid7", "D1_mex", "D1_out"}},
        {"shms", "Events reaching hut", {"hut"}},
        {"shms", "DC1, DC2, Scin, Cal", {"dc1", "dc2"}},
        {"shms", "S1, S2, S3, Cal", {"s1", "s2", "s3", "cal"}},
        {"shms", "Successes", {"successes"}},
        // Calo (simc.f:532-535)
        {"calo", "Calo Trials:", {"trials"}},
        {"calo", "Extent hor/vert", {"slit_hor", "slit_vert"}},
        {"calo", "Successes", {"successes"}},
    };
    return t;
}

bool isInteger(const std::string& t) {
    if (t.empty()) return false;
    size_t i = (t[0] == '-' || t[0] == '+') ? 1 : 0;
    if (i >= t.size()) return false;
    for (; i < t.size(); ++i)
        if (!std::isdigit(static_cast<unsigned char>(t[i]))) return false;
    return true;
}

// ---------------------------------------------------------------------------
// JSON helpers: tolerant key lookup (the schema of docs/summary.schema.json is
// owned by P1.4; aliases keep this reader working if names differ slightly).
// ---------------------------------------------------------------------------
const json* findKey(const json& j, const std::vector<std::string>& aliases, int depth = 3) {
    if (!j.is_object()) return nullptr;
    for (const auto& a : aliases) {
        auto it = j.find(a);
        if (it != j.end()) return &(*it);
    }
    if (depth <= 0) return nullptr;
    for (auto it = j.begin(); it != j.end(); ++it) {
        if (it.key() == "stops" || it.key() == "stop_counters") continue;
        if (const json* r = findKey(*it, aliases, depth - 1)) return r;
    }
    return nullptr;
}

bool jsonNumber(const json* p, double& v) {
    if (!p) return false;
    if (p->is_number()) { v = p->get<double>(); return true; }
    if (p->is_object() && p->contains("value") && (*p)["value"].is_number()) {
        v = (*p)["value"].get<double>();
        return true;
    }
    if (p->is_string()) {
        try { v = std::stod(p->get<std::string>()); return true; } catch (...) {}
    }
    return false;
}

}  // namespace

// ---------------------------------------------------------------------------
void parseFortranHist(const std::string& path, RunSummary& s) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    s.source = path;
    s.format = "fortran-hist";
    s.producer = "fortran";
    // 'label = value [unit]' lines of simc.f:report
    static const std::regex kv(R"(^\s*([A-Za-z][A-Za-z0-9_.:() %\-]*?)\s*=\s*([-+0-9.EeDd*]+)\s*(\S*)\s*$)");
    static const std::regex flag(R"(([A-Za-z_]+)\s*=\s*([TF])\b)");
    std::string line;
    int nlumi = 0;
    bool haveNgen = false, haveNtried = false, haveNcontr = false;
    // P12.3: pass_cuts window (simc.f:report "Input Spectrometer Limits" and "slop.MC" lines) and electron_arm
    static const std::regex spe(R"(^\s*SPedge\.([ep])\.(delta|yptar|xptar)\.min/max\s*=\s*(\S+)\s+(\S+))");
    static const std::regex slp(R"(^\s*(?:slop\.MC)?\s+([EP]) arm (delta|yptar|xptar)\s+(\S+))");
    static const std::regex earm(R"(electron_arm=\s*(-?\d+))");
    bool inSlop = false;
    int nSpedge = 0, nSlop = 0;
    s.fortran_compat = 1;
    auto vidx = [](const std::string& a, const std::string& v) {
        return (a == "p" || a == "P" ? 3 : 0) + (v == "delta" ? 0 : (v == "yptar" ? 1 : 2));
    };
    while (std::getline(in, line)) {
        {
            std::smatch q;
            if (std::regex_search(line, q, spe)) {
                const int k = vidx(q[1], q[2]);
                try {
                    s.spedge_min[k] = std::stod(q[3]);
                    s.spedge_max[k] = std::stod(q[4]);
                    ++nSpedge;
                } catch (...) {}
            }
            if (line.find("slop.MC") != std::string::npos) inSlop = true;
            if (inSlop && std::regex_search(line, q, slp)) {
                try {
                    s.slop_used[vidx(q[1], q[2])] = std::stod(q[3]);
                    ++nSlop;
                } catch (...) {}
                if (nSlop >= 6) inSlop = false;
            }
            if (s.electron_arm == 0 && std::regex_search(line, q, earm)) s.electron_arm = std::stoi(q[1]);
        }
        // flags (FLAGS: block, 'name= T')
        for (auto it = std::sregex_iterator(line.begin(), line.end(), flag); it != std::sregex_iterator(); ++it) {
            std::string n = (*it)[1];
            if (n.rfind("doing_", 0) == 0 || n.rfind("using_", 0) == 0 || n == "mc_smear" ||
                n.rfind("correct_", 0) == 0)
                s.flags[n] = ((*it)[2] == "T");
        }
        std::smatch m;
        if (!std::regex_match(line, m, kv)) continue;
        std::string label = collapse(m[1]);
        std::string val = m[2];
        std::string unit = m[3];
        auto setI = [&](int64_t& dst, bool& have) {
            try { dst = std::stoll(val); have = true; } catch (...) {}
        };
        if (label == "Ngen (request)") setI(s.ngen, haveNgen);
        else if (label == "Ntried") setI(s.ntried, haveNtried);
        else if (label == "Ncontribute") setI(s.ncontribute, haveNcontr);
        else if (label == "Nco_no_rad_prot") { bool h; setI(s.nco_no_rad_prot, h); }
        else if (label == "Random Seed") setI(s.random_seed, s.has_seed);
        else if (label == "MeV: wtcontr") s.scalars["wtcontr"] = parseFortranReal(val);
        else if (label == "CENTRAL.sigcc") s.scalars["sigcc_central"] = parseFortranReal(val);
        else if (label == "AVERAGE.sigcc") s.scalars["sigcc_average"] = parseFortranReal(val);
        else if (label == "charge") s.scalars["charge_mC"] = parseFortranReal(val);
        else if (label == "targetfac") s.scalars["targetfac"] = parseFortranReal(val);
        else if (label == "luminosity") {
            // simc.f:961 (ub^-1) then simc.f:962 (GeV^2)
            if (unit == "GeV^2" || nlumi == 1) s.scalars["luminosity_GeV2"] = parseFortranReal(val);
            else s.scalars["luminosity_ub"] = parseFortranReal(val);
            ++nlumi;
        } else if (label == "genvol") s.scalars["genvol"] = parseFortranReal(val);
        else if (label == "normfac") s.scalars["normfac"] = parseFortranReal(val);
    }
    s.has_counters = haveNgen && haveNtried && haveNcontr;
    if (!s.has_counters) s.warnings.push_back("COUNTERS block incomplete in " + path);
    s.has_spedge = nSpedge == 6;
    s.has_slop = nSlop == 6;
}

// ---------------------------------------------------------------------------
void parseFortranGeni(const std::string& path, RunSummary& s) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line, cur;
    const auto& tab = lineTable();
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        std::vector<std::string> tok;
        std::string t;
        while (ss >> t) tok.push_back(t);
        if (tok.empty()) continue;
        size_t k = tok.size();
        while (k > 0 && isInteger(tok[k - 1])) --k;
        if (k == tok.size()) {
            // Header of the distributions part (not produced by simc.f in .geni,
            // but stop reading at any non-counter text block).
            continue;
        }
        std::string label;
        for (size_t i = 0; i < k; ++i) label += (i ? " " : "") + tok[i];
        std::vector<int64_t> vals;
        for (size_t i = k; i < tok.size(); ++i) vals.push_back(std::stoll(tok[i]));

        // spectrometer header?
        for (const auto& lm : tab)
            if (label == lm.label && std::string(lm.label).find("Trials:") != std::string::npos) cur = lm.spec;
        if (cur.empty()) {
            s.warnings.push_back(".geni: counter line before any 'Trials:' header: " + line);
            continue;
        }
        const LineMap* hit = nullptr;
        std::string fam = (cur == "hrsr" || cur == "hrsl") ? "hrs" : cur;
        for (const auto& lm : tab) {
            if (label != lm.label) continue;
            if (lm.spec == cur || lm.spec == fam) { hit = &lm; break; }
        }
        if (!hit) {
            // unknown line: keep it under a label-derived key
            std::string base = lower(label);
            for (auto& c : base) if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
            for (size_t i = 0; i < vals.size(); ++i) {
                std::string key = base + "_" + std::to_string(i);
                s.stops[cur][key] = vals[i];
                s.stopOrder[cur].push_back(key);
            }
            s.warnings.push_back(".geni: unmapped line in " + cur + ": '" + label + "'");
            continue;
        }
        if (hit->keys.size() != vals.size())
            s.warnings.push_back(".geni: '" + label + "' has " + std::to_string(vals.size()) +
                                 " values, expected " + std::to_string(hit->keys.size()));
        for (size_t i = 0; i < std::min(vals.size(), hit->keys.size()); ++i) {
            s.stops[cur][hit->keys[i]] = vals[i];
            s.stopOrder[cur].push_back(hit->keys[i]);
        }
    }
}

// ---------------------------------------------------------------------------
void parseSummaryJson(const std::string& path, RunSummary& s) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    json j;
    try {
        in >> j;
    } catch (const std::exception& e) {
        throw std::runtime_error("JSON parse error in " + path + ": " + e.what());
    }
    s.source = path;
    s.format = "summary-json";
    if (const json* p = findKey(j, {"producer", "program"}); p && p->is_string()) s.producer = *p;
    if (const json* p = findKey(j, {"tag", "run_tag", "base"}); p && p->is_string()) s.tag = *p;
    // values converted from a Fortran .hist carry the .hist print precision
    const bool fromText = lower(s.producer).find("fortran") != std::string::npos ||
                          lower(s.producer).find("hist2json") != std::string::npos;

    double v;
    bool ok = true;
    auto geti = [&](const std::vector<std::string>& al, int64_t& dst) {
        if (jsonNumber(findKey(j, al), v)) { dst = static_cast<int64_t>(std::llround(v)); return true; }
        return false;
    };
    ok &= geti({"ngen", "Ngen", "ngen_request"}, s.ngen);
    ok &= geti({"ntried", "Ntried"}, s.ntried);
    ok &= geti({"ncontribute", "Ncontribute"}, s.ncontribute);
    geti({"nco_no_rad_prot", "Nco_no_rad_prot", "ncontribute_no_rad_proton"}, s.nco_no_rad_prot);
    s.has_seed = geti({"random_seed", "seed"}, s.random_seed);
    s.has_counters = ok;
    if (!ok) s.warnings.push_back("summary JSON lacks ngen/ntried/ncontribute: " + path);

    const std::vector<std::pair<std::string, std::vector<std::string>>> sc = {
        {"wtcontr", {"wtcontr", "wtcontribute"}},
        {"sigcc_central", {"sigcc_central", "central_sigcc", "CENTRAL.sigcc"}},
        {"sigcc_average", {"sigcc_average", "average_sigcc", "ave_sigcc", "AVERAGE.sigcc"}},
        {"charge_mC", {"charge_mC", "charge"}},
        {"targetfac", {"targetfac"}},
        {"luminosity_ub", {"luminosity_ub", "luminosity_ubinv", "luminosity"}},
        {"luminosity_GeV2", {"luminosity_GeV2", "luminosity_gev2"}},
        {"genvol", {"genvol"}},
        {"normfac", {"normfac"}},
    };
    for (const auto& [key, al] : sc) {
        if (jsonNumber(findKey(j, al), v)) {
            Scalar x;
            x.value = v;
            x.present = true;
            if (fromText && v != 0.0) {
                // e16.8 for wtcontr (simc.f:920), e16.6 otherwise (simc.f:9915)
                int nd = (key == "wtcontr") ? 8 : 6;
                int ex = static_cast<int>(std::floor(std::log10(std::fabs(v)))) + 1;
                x.halfulp = 0.5 * std::pow(10.0, ex - nd);
            }
            s.scalars[key] = x;
        }
    }
    if (const json* f = findKey(j, {"flags"}); f && f->is_object()) {
        for (auto it = f->begin(); it != f->end(); ++it) {
            if (it->is_boolean()) s.flags[it.key()] = it->get<bool>();
            else if (it->is_number_integer()) s.flags[it.key()] = it->get<int>() != 0;
        }
    }

    const json* st = nullptr;
    for (const char* k : {"stops", "stop_counters"})
        if (j.contains(k)) st = &j[k];
    if (!st) st = findKey(j, {"stops", "stop_counters"}, 1);
    if (st && st->is_object()) {
        for (auto it = st->begin(); it != st->end(); ++it) {
            if (!it->is_object()) continue;
            std::string spec = lower(it.key());
            for (auto c = it->begin(); c != it->end(); ++c) {
                if (!jsonNumber(&(*c), v)) continue;
                s.stops[spec][c.key()] = static_cast<int64_t>(std::llround(v));
                s.stopOrder[spec].push_back(c.key());
            }
        }
    }

    // P12.3: pass_cuts window ("spedge" [min, max], "slop".MC[].used), settings.electron_arm, run.fortran_compat
    if (j.contains("spedge") && j["spedge"].is_object()) {
        int n = 0;
        const char* arm[2] = {"e", "p"};
        const char* var[3] = {"delta", "yptar", "xptar"};
        for (int a = 0; a < 2; ++a)
            for (int k = 0; k < 3; ++k) {
                const json& sp = j["spedge"];
                if (sp.contains(arm[a]) && sp[arm[a]].contains(var[k]) && sp[arm[a]][var[k]].is_array() &&
                    sp[arm[a]][var[k]].size() == 2) {
                    s.spedge_min[3 * a + k] = sp[arm[a]][var[k]][0].get<double>();
                    s.spedge_max[3 * a + k] = sp[arm[a]][var[k]][1].get<double>();
                    ++n;
                }
            }
        s.has_spedge = n == 6;
    }
    if (j.contains("slop") && j["slop"].is_object() && j["slop"].contains("MC") && j["slop"]["MC"].is_array()) {
        int n = 0;
        static const char* labels[6] = {"E arm delta", "E arm yptar", "E arm xptar",
                                        "P arm delta", "P arm yptar", "P arm xptar"};
        for (const auto& e : j["slop"]["MC"]) {
            if (!e.is_object() || !e.contains("label") || !e.contains("used")) continue;
            for (int k = 0; k < 6; ++k)
                if (e["label"] == labels[k]) {
                    s.slop_used[k] = e["used"].get<double>();
                    ++n;
                }
        }
        s.has_slop = n == 6;
    }
    if (j.contains("settings") && j["settings"].is_object() && j["settings"].contains("electron_arm") &&
        j["settings"]["electron_arm"].is_number())
        s.electron_arm = j["settings"]["electron_arm"].get<int>();
    if (j.contains("run") && j["run"].is_object() && j["run"].contains("fortran_compat") &&
        j["run"]["fortran_compat"].is_boolean())
        s.fortran_compat = j["run"]["fortran_compat"].get<bool>() ? 1 : 0;
    else if (fromText)
        s.fortran_compat = 1;   // converted from a Fortran .hist

    if (const json* t = findKey(j, {"tree"}); t && t->is_object()) {
        if (t->contains("file") && (*t)["file"].is_string()) s.treeFile = (*t)["file"];
        if (t->contains("name") && (*t)["name"].is_string()) s.treeName = (*t)["name"];
        if (t->contains("weight_branch") && (*t)["weight_branch"].is_string()) s.weightBranch = (*t)["weight_branch"];
    } else if (const json* t2 = findKey(j, {"root_file", "tree_file"}); t2 && t2->is_string()) {
        s.treeFile = *t2;
    }
    if (!s.treeFile.empty() && fs::path(s.treeFile).is_relative())
        s.treeFile = (fs::path(path).parent_path() / s.treeFile).string();
}

// ---------------------------------------------------------------------------
RunSummary loadRun(const std::string& spec, bool preferJson) {
    RunSummary s;
    fs::path p(spec);
    if (!fs::exists(p)) throw std::runtime_error("run spec not found: " + spec);

    auto stemOf = [](const fs::path& f) {
        std::string n = f.filename().string();
        for (const char* suf : {".summary.json", ".json", ".hist", ".geni", ".root"}) {
            std::string sf = suf;
            if (n.size() > sf.size() && n.compare(n.size() - sf.size(), sf.size(), sf) == 0)
                return (f.parent_path() / n.substr(0, n.size() - sf.size())).string();
        }
        return f.string();
    };

    fs::path jsonFile, histFile;
    std::string stem;
    if (fs::is_directory(p)) {
        std::vector<fs::path> js, hs;
        for (auto& e : fs::directory_iterator(p)) {
            std::string n = e.path().filename().string();
            if (n.size() > 13 && n.compare(n.size() - 13, 13, ".summary.json") == 0) js.push_back(e.path());
            if (e.path().extension() == ".hist") hs.push_back(e.path());
        }
        if (js.size() > 1 || hs.size() > 1)
            throw std::runtime_error("directory " + spec + " has several runs; give a file");
        if (!js.empty()) jsonFile = js[0];
        if (!hs.empty()) histFile = hs[0];
        stem = stemOf(!jsonFile.empty() ? jsonFile : histFile);
    } else {
        stem = stemOf(p);
        std::string ext = p.extension().string();
        if (ext == ".json") jsonFile = p;
        else if (ext == ".hist") histFile = p;
    }
    if (jsonFile.empty() && fs::exists(stem + ".summary.json")) jsonFile = stem + ".summary.json";
    if (histFile.empty() && fs::exists(stem + ".hist")) histFile = stem + ".hist";

    bool useJson = !jsonFile.empty() && (preferJson || histFile.empty());
    if (p.extension() == ".hist") useJson = false;
    if (p.extension() == ".json") useJson = true;
    if (useJson) {
        parseSummaryJson(jsonFile.string(), s);
    } else if (!histFile.empty()) {
        parseFortranHist(histFile.string(), s);
        std::string geni = stemOf(histFile) + ".geni";
        if (fs::exists(geni)) parseFortranGeni(geni, s);
        else s.warnings.push_back("no .geni next to " + histFile.string() + " (no stop counters)");
    } else {
        s.source = spec;
        s.format = "tree-only";
        s.warnings.push_back("no summary JSON or .hist found for " + spec + "; run-level tests skipped");
    }
    if (s.tag.empty()) s.tag = fs::path(stem).filename().string();
    if (s.treeFile.empty()) {
        if (p.extension() == ".root") s.treeFile = p.string();
        else if (fs::exists(stem + ".root")) s.treeFile = stem + ".root";
    }
    return s;
}

// ---------------------------------------------------------------------------
void writeSummaryJson(const RunSummary& s, const std::string& path) {
    json j;
    j["schema"] = "simc_run_summary";
    j["schema_version"] = 1;
    j["run"] = {{"tag", s.tag}, {"producer", s.format == "fortran-hist" ? "fortran-hist(simc_compare)" : s.producer}};
    if (s.has_seed) j["run"]["random_seed"] = s.random_seed;
    j["counters"] = {{"ngen", s.ngen}, {"ntried", s.ntried}, {"ncontribute", s.ncontribute},
                     {"nco_no_rad_prot", s.nco_no_rad_prot}};
    json n = json::object();
    for (const auto& [k, v] : s.scalars)
        if (v.present) n[k] = v.value;
    j["normalization"] = n;
    json f = json::object();
    for (const auto& [k, v] : s.flags) f[k] = v;
    j["flags"] = f;
    json st = json::object();
    for (const auto& [spec, keys] : s.stopOrder) {
        json o = json::object();
        for (const auto& k : keys) o[k] = s.stops.at(spec).at(k);
        st[spec] = o;
    }
    j["stops"] = st;
    j["tree"] = {{"file", fs::path(s.treeFile).filename().string()}, {"name", s.treeName},
                 {"weight_branch", s.weightBranch}};
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write " + path);
    o << j.dump(2) << "\n";
}

}  // namespace simcval
