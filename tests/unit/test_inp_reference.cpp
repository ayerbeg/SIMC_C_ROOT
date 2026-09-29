// tests/unit/test_inp_reference.cpp — .inp conversion against the Fortran reference runs.
//   all_infiles <dir>            : every <dir>/*.inp converts (with its extra file) without errors
//   hist <file.inp> <file.hist>  : converted values == values echoed by the Fortran in the .hist
//   notreg <file.inp> <file.out> <n> : "not registered" list == the Fortran log (n names, same order)
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "TestUtil.h"
#include "simc/input/InpReader.h"

using namespace simc;
using simctest::check;

namespace {

std::vector<std::string> readLines(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::string> v;
    for (std::string l; std::getline(in, l);) v.push_back(l);
    return v;
}

// Absolute tolerance implied by the way a number was printed (half a unit in the last digit, +1%).
double printTolerance(const std::string& tok) {
    std::string mant = tok;
    int expo = 0;
    auto e = tok.find_first_of("eEdD");
    if (e != std::string::npos) {
        mant = tok.substr(0, e);
        expo = std::atoi(tok.c_str() + e + 1);
    }
    auto dot = mant.find('.');
    int dec = dot == std::string::npos ? 0 : static_cast<int>(mant.size() - dot - 1);
    return 0.505 * std::pow(10.0, expo - dec);
}

struct Hist {
    std::vector<std::string> lines;
    // First match of `re` (ECMAScript); returns captured groups as strings.
    bool find(const std::string& re, std::vector<std::string>& groups, size_t from = 0) const {
        std::regex r(re);
        for (size_t i = from; i < lines.size(); ++i) {
            std::smatch m;
            if (std::regex_search(lines[i], m, r)) {
                groups.clear();
                for (size_t g = 1; g < m.size(); ++g) groups.push_back(m[g]);
                return true;
            }
        }
        return false;
    }
};

std::vector<std::string> g_verified;

void num(const Hist& h, const std::string& label, const std::string& re, const std::vector<double>& expect) {
    std::vector<std::string> g;
    if (!h.find(re, g)) {
        check(false, "label not found in .hist: " + label);
        return;
    }
    check(g.size() == expect.size(), label + ": number of values");
    bool all = true;
    for (size_t i = 0; i < g.size() && i < expect.size(); ++i) {
        const double v = std::strtod(g[i].c_str(), nullptr);
        const double tol = printTolerance(g[i]);
        const bool ok = std::abs(v - expect[i]) <= tol;
        all = all && ok;
        check(ok, label + "[" + std::to_string(i) + "]: .hist " + g[i] + " vs converted " + std::to_string(expect[i]));
    }
    if (all) g_verified.push_back(label);
}

void flag(const Hist& h, const std::string& name, bool expect, const std::string& label = "") {
    std::vector<std::string> g;
    if (!h.find("\\b" + name + "=\\s*([TF])\\b", g)) {
        check(false, "flag not found in .hist: " + name);
        return;
    }
    const bool ok = (g[0] == "T") == expect;
    check(ok, name + ": .hist " + g[0] + " vs converted " + (expect ? "T" : "F"));
    if (ok) g_verified.push_back(label.empty() ? name : label);
}

void ival(const Hist& h, const std::string& name, int expect, const std::string& label = "") {
    std::vector<std::string> g;
    if (!h.find("\\b" + name + "=\\s*(-?\\d+)\\b", g)) {
        check(false, "int not found in .hist: " + name);
        return;
    }
    const bool ok = std::atoi(g[0].c_str()) == expect;
    check(ok, name + ": .hist " + g[0] + " vs converted " + std::to_string(expect));
    if (ok) g_verified.push_back(label.empty() ? name : label);
}

int histCheck(const std::string& inp, const std::string& histfile) {
    if (!simctest::exists(inp) || !simctest::exists(histfile)) {
        std::printf("SKIP: reference files missing\n");
        return simctest::kSkip;
    }
    SimcParameters p;
    InpReader r(p);
    check(r.loadRun(inp), "load " + inp);
    Hist h{readLines(histfile)};

    // ---- KINEMATICS (simc.f:report)
    num(h, "Ebeam", "^\\s*Ebeam =\\s*(\\S+)", {p.Ebeam});
    num(h, "dEbeam", "\\(dE/E\\)beam =\\s*(\\S+)", {p.dEbeam / 100.0});
    num(h, "gen%xwid", "x-width =\\s*(\\S+)", {p.gen.xwid});
    num(h, "gen%ywid", "y-width =\\s*(\\S+)", {p.gen.ywid});
    num(h, "targ%fr_pattern", "fr_pattern =\\s*(\\S+)", {double(p.targ.fr_pattern)});
    num(h, "targ%fr1", "^\\s*fr1 =\\s*(\\S+)", {p.targ.fr1});
    num(h, "targ%fr2", "^\\s*fr2 =\\s*(\\S+)", {p.targ.fr2});
    num(h, "spec%e/p%theta", "^\\s*angle =\\s*(\\S+)\\s+(\\S+)\\s+deg\\s*$",
        {std::abs(p.spec.e.theta), std::abs(p.spec.p.theta)});
    num(h, "spec%e/p%P", "^\\s*momentum =\\s*(\\S+)\\s+(\\S+)", {p.spec.e.P, p.spec.p.P});
    num(h, "spec%e/p%offset%x", "^\\s*x offset =\\s*(\\S+)\\s+(\\S+)", {p.spec.e.offset.x, p.spec.p.offset.x});
    num(h, "spec%e/p%offset%y", "^\\s*y offset =\\s*(\\S+)\\s+(\\S+)", {p.spec.e.offset.y, p.spec.p.offset.y});
    num(h, "spec%e/p%offset%z", "^\\s*z offset =\\s*(\\S+)\\s+(\\S+)", {p.spec.e.offset.z, p.spec.p.offset.z});
    num(h, "spec%e/p%offset%xptar", "xptar offset =\\s*(\\S+)\\s+(\\S+)",
        {p.spec.e.offset.xptar, p.spec.p.offset.xptar});
    num(h, "spec%e/p%offset%yptar", "yptar offset =\\s*(\\S+)\\s+(\\S+)",
        {p.spec.e.offset.yptar, p.spec.p.offset.yptar});
    // ---- TARGET
    num(h, "targ%A,targ%Z", "^\\s*A =\\s*(\\S+)\\s+Z =\\s*(\\S+)", {p.targ.A, p.targ.Z});
    num(h, "targ%mass_amu", "^\\s*mass =\\s*(\\S+)\\s+amu", {p.targ.mass_amu});
    num(h, "targ%rho,targ%thick", "rho =\\s*(\\S+) g/cm3\\s+thick =\\s*(\\S+) g/cm2",
        {p.targ.rho, p.targ.thick / 1000.0});
    const double tangle = (p.targ.Z < 2.4) ? 0.0 : p.targ.angle;   // dbase.f:456-460
    num(h, "targ%angle,targ%abundancy", "^\\s*angle =\\s*(\\S+)\\s+deg\\s+abundancy =\\s*(\\S+)",
        {tangle, p.targ.abundancy});
    num(h, "targ%zoffset", "zoffset =\\s*(\\S+)", {p.targ.zoffset});
    num(h, "targ%xoffset,targ%yoffset", "xoffset =\\s*(\\S+)\\s+cm\\s+yoffset =\\s*(\\S+)",
        {p.targ.xoffset, p.targ.yoffset});
    // ---- FLAGS (doing_phsp is off in both reference runs; no cross-flag rewriting needed)
    flag(h, "doing_kaon", p.doing_kaon);
    flag(h, "doing_pion", p.doing_pion);
    flag(h, "doing_semi", p.doing_semi);
    flag(h, "doing_rho", p.doing_rho);
    flag(h, "doing_hplus", p.doing_hplus);
    flag(h, "doing_pizero", p.doing_pizero);
    flag(h, "doing_delta", p.doing_delta);
    flag(h, "doing_phsp", p.doing_phsp);
    flag(h, "do_fermi", p.do_fermi);
    flag(h, "mc_smear", p.mc_smear);
    flag(h, "using_Eloss", p.using_Eloss);
    flag(h, "using_Coulomb", p.using_Coulomb && std::lround(p.targ.Z) != 1, "using_Coulomb (Z=1 rule dbase.f:547)");
    flag(h, "deForest_flag", p.deForest_flag != 0, "deForest_flag (printed with L edit)");
    flag(h, "correct_Eloss", p.correct_Eloss && p.using_Eloss, "correct_Eloss (dbase.f:542)");
    flag(h, "correct_raster", p.correct_raster);
    flag(h, "doing_decay", p.doing_decay);
    flag(h, "use_benhar_sf", p.use_benhar_sf);
    flag(h, "using_E_arm_monteca", p.spect_mode != 1 && p.spect_mode != -1, "spect_mode (E arm, dbase.f:531)");
    flag(h, "using_P_arm_monteca", p.spect_mode != 1 && p.spect_mode != -2, "spect_mode (P arm, dbase.f:532)");
    {
        std::vector<std::string> g;
        if (h.find("\\buse_first_cer=\\s*([TF])", g)) flag(h, "use_first_cer", p.use_first_cer != 0);
    }
    ival(h, "which_pion", p.which_pion);
    ival(h, "which_kaon", p.which_kaon);
    ival(h, "pizero_ngamma", p.pizero_ngamma);
    ival(h, "electron_arm", p.electron_arm);
    ival(h, "hadron_arm", p.hadron_arm);
    num(h, "ctau", "ctau=\\s*(\\S+)", {p.ctau});
    num(h, "ngen", "Ngen \\(request\\) =\\s*(\\S+)", {double(p.ngen)});
    // ---- RADIATIVE CORRECTIONS
    flag(h, "use_expon", p.use_expon != 0, "use_expon (printed with L edit)");
    flag(h, "using_rad", p.using_rad);
    flag(h, "use_offshell_r", p.use_offshell_rad, "use_offshell_rad");
    ival(h, "rad_flag", p.rad_flag);
    int extrad = p.extrad_flag;
    if (extrad == 0) extrad = (p.rad_flag == 0) ? 3 : 1;   // init.f:621-626
    ival(h, "extrad_flag", extrad);
    ival(h, "one_tail", p.one_tail);
    ival(h, "intcor_mode", p.intcor_mode);
    num(h, "dE_edge_test", "dE_edge_test=\\s*(\\S+)", {p.dE_edge_test});
    num(h, "EXPER%charge", "^\\s*charge =\\s*(\\S+)\\s+mC", {p.EXPER.charge});
    // ---- SPedge: printed after widening by slop%MC%*%used (init.f:209-216)
    {
        size_t k = 0;
        for (; k < h.lines.size(); ++k)
            if (h.lines[k].find("slop.MC") != std::string::npos) break;
        check(k + 5 < h.lines.size(), "slop.MC table present");
        if (k + 5 < h.lines.size()) {
            double slop[6];
            for (int i = 0; i < 6; ++i) {
                std::smatch m;
                std::regex re("arm (delta|yptar|xptar)\\s+(\\S+)");
                std::regex_search(h.lines[k + i], m, re);
                slop[i] = std::strtod(m[2].str().c_str(), nullptr);
            }
            auto edge = [&](const std::string& arm, const std::string& var, const ParamRange& rg, double s) {
                const bool ang = var != "delta";
                double lo = rg.min, hi = rg.max;
                if (!ang && lo <= -100.0) lo = -99.99;           // dbase.f:479-480
                const double f = ang ? 1000.0 : 1.0;            // mrad -> rad (dbase.f:481-488)
                num(h, "SPedge%" + arm + "%" + var + "%min/max",
                    "SPedge\\." + arm + "\\." + var + "\\.min/max =\\s*(\\S+)\\s+(\\S+)",
                    {lo / f - s / f, hi / f + s / f});
            };
            edge("e", "delta", p.SPedge.e.delta, slop[0]);
            edge("e", "yptar", p.SPedge.e.yptar, slop[1]);
            edge("e", "xptar", p.SPedge.e.xptar, slop[2]);
            edge("p", "delta", p.SPedge.p.delta, slop[3]);
            edge("p", "yptar", p.SPedge.p.yptar, slop[4]);
            edge("p", "xptar", p.SPedge.p.xptar, slop[5]);
        }
    }
    std::printf("verified against .hist (%zu):", g_verified.size());
    for (const auto& v : g_verified) std::printf(" %s;", v.c_str());
    std::printf("\n");
    return simctest::finish("inp_vs_hist");
}

int notreg(const std::string& inp, const std::string& out, int n_expected) {
    if (!simctest::exists(inp) || !simctest::exists(out)) {
        std::printf("SKIP: reference files missing\n");
        return simctest::kSkip;
    }
    std::vector<std::string> fortran;
    std::regex re("^(\\S+) not registered, registering as int\\s*$");
    for (const auto& l : readLines(out)) {
        std::smatch m;
        if (std::regex_match(l, m, re)) fortran.push_back(m[1]);
    }
    SimcParameters p;
    InpReader r(p);
    check(r.loadRun(inp), "load " + inp);
    const auto& mine = r.notRegistered();
    check(static_cast<int>(fortran.size()) == n_expected,
          "Fortran log has " + std::to_string(n_expected) + " names (got " + std::to_string(fortran.size()) + ")");
    check(mine.size() == fortran.size(),
          "same count: C++ " + std::to_string(mine.size()) + " vs Fortran " + std::to_string(fortran.size()));
    for (size_t i = 0; i < std::min(mine.size(), fortran.size()); ++i)
        check(mine[i] == fortran[i], "name " + std::to_string(i) + ": C++ '" + mine[i] + "' vs Fortran '" + fortran[i] + "'");
    return simctest::finish("inp_not_registered");
}

int allInfiles(const std::string& dir) {
    DIR* d = opendir(dir.c_str());
    if (!d) {
        std::printf("SKIP: %s not found\n", dir.c_str());
        return simctest::kSkip;
    }
    std::vector<std::string> files;
    while (dirent* e = readdir(d)) {
        std::string n = e->d_name;
        if (n.size() > 4 && n.compare(n.size() - 4, 4, ".inp") == 0) files.push_back(n);
    }
    closedir(d);
    std::sort(files.begin(), files.end());
    check(!files.empty(), "infiles found");
    for (const auto& f : files) {
        SimcParameters p;
        InpReader r(p);
        const bool ok = r.loadRun(dir + "/" + f);
        check(ok, "convert " + f);
        int errors = 0;
        for (const auto& m : r.messages())
            if (m.text.find("Tried to fill") != std::string::npos || m.text.find("Error interpreting") != std::string::npos ||
                m.text.find("Failed to open") != std::string::npos || m.text.find("can't be a number") != std::string::npos ||
                m.text.find("Unknown block type") != std::string::npos)
                ++errors, std::printf("  %s:%d: %s\n", m.file.c_str(), m.line, m.text.c_str());
        check(errors == 0, f + ": no CTP parse errors");
        // Round trip through JSON must be lossless.
        SimcParameters q;
        FromJson(ToJson(p), q);
        check(ToJson(q) == ToJson(p), f + ": JSON round trip");
        std::printf("%-40s files=%zu set=%zu not_registered=%zu", f.c_str(), r.filesLoaded().size(),
                    r.provenance().size(), r.notRegistered().size());
        if (!r.notRegistered().empty()) {
            std::printf(" [");
            for (size_t i = 0; i < r.notRegistered().size(); ++i) std::printf("%s%s", i ? " " : "", r.notRegistered()[i].c_str());
            std::printf("]");
        }
        std::printf("\n");
    }
    return simctest::finish("inp_convert_all_infiles");
}

}  // namespace

int main(int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "";
    if (mode == "all_infiles" && argc > 2) return allInfiles(argv[2]);
    if (mode == "hist" && argc > 3) return histCheck(argv[2], argv[3]);
    if (mode == "notreg" && argc > 4) return notreg(argv[2], argv[3], std::atoi(argv[4]));
    std::printf("usage: test_inp_reference all_infiles <dir> | hist <inp> <hist> | notreg <inp> <out> <n>\n");
    return 2;
}
