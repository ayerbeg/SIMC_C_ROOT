// tests/l1/test_l1_driver.cpp
//
// L1 test of simc::driver::MonteCarlo (task P7.1, simc.f:montecarlo) against the P4.F1 captures of the unmodified
// SIMC program (reference/functions/event_capture/<cfg>/, reference/runs/p4f1_<cfg>_cap/):
//
//   replay <cap_dir> <run_dir> <cfg> <data_dir> [run_name]   (data_dir: the theory files of D/A(e,e'p) captures)
//     * RunSetup from the run's .inp (InpReader + BuildRunSetup), checked against cap_final (spec, Mh2, SPedge);
//     * every recorded montecarlo call, in order, with ONE persistent main/recon (as in SIMC, so the stale values
//       of failed tries are reproduced): orig and main%target from cap_generate, the phase-2 random numbers of
//       cap_rng in a ScriptedRandom; after a success, recon is updated with the recorded complete_recon_ev
//       outputs (cap_recon) — complete_recon_ev overwrites recon%p%P/E for phsp;
//     * compared bitwise: mc_ok, fail arm, HMS/SHMS stop codes, SHMS stop_id, grnd/gauss1 counts, main%SP,
//       main%RECON, main%FP, ntup%resfac, and the recon arm fields (rec_*);
//     * the 46 .bin columns (results_write.f, electron_arm = 5 layout, H elastic) of every written event, built
//       from the C++ montecarlo outputs (recon arm fields, main%FP) plus the recorded fields of the other stages
//       (cap_events: orig, main%target, sigcc, weight, radphot; cap_recon: the complete_recon_ev quantities),
//       compared bit for bit with the archived .bin.
// Exit: 0 pass, 1 fail, 77 reference missing.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "simc/core/Json.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/MonteCarlo.h"
#include "simc/input/InpReader.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/setup/SetupBuilder.h"
#include "simc/spectrometers/hms/kernel/McHms.h"
#include "simc/spectrometers/shms/ShmsSpectrometer.h"
#include "simc/spectrometers/shms/kernel/McShms.h"

using namespace simc;

namespace {

// ------------------------------------------------------------------------------------------------ tables
struct Tab {
    std::vector<std::string> names;
    std::map<std::string, int> idx;
    std::vector<std::vector<double>> rows;
    int col(const std::string& n) const {
        auto it = idx.find(n);
        if (it == idx.end()) throw std::runtime_error("missing column " + n);
        return it->second;
    }
};

bool exists(const std::string& p) {
    std::ifstream f(p);
    return static_cast<bool>(f);
}

/// Load a whitespace table (plain or .xz through `xz -dc`). Header = first non-# line.
bool loadTab(const std::string& base, Tab& t) {
    std::string cmd;
    FILE* f = nullptr;
    bool piped = false;
    if (exists(base)) {
        f = std::fopen(base.c_str(), "r");
    } else if (exists(base + ".xz")) {
        cmd = "xz -dc '" + base + ".xz'";
        f = popen(cmd.c_str(), "r");
        piped = true;
    }
    if (!f) return false;
    std::string line;
    char buf[1 << 16];
    bool header = false;
    std::string acc;
    while (std::fgets(buf, sizeof buf, f)) {
        acc += buf;
        if (acc.empty() || acc.back() != '\n') continue;   // long line: keep reading
        line.swap(acc);
        acc.clear();
        std::size_t p = line.find_first_not_of(" \t");
        if (p == std::string::npos || line[p] == '#') continue;
        std::istringstream is(line);
        if (!header) {
            std::string n;
            while (is >> n) {
                t.idx[n] = static_cast<int>(t.names.size());
                t.names.push_back(n);
            }
            header = true;
            continue;
        }
        std::vector<double> r;
        r.reserve(t.names.size());
        std::string tok;
        while (is >> tok) r.push_back(std::strtod(tok.c_str(), nullptr));   // handles E+000, nan, inf
        if (r.size() != t.names.size()) {
            std::printf("  %s: row with %zu values, header %zu\n", base.c_str(), r.size(), t.names.size());
            if (piped) pclose(f); else std::fclose(f);
            return false;
        }
        t.rows.push_back(std::move(r));
    }
    if (piped) pclose(f); else std::fclose(f);
    return header;
}

/// Read a Fortran sequential unformatted file: record payloads.
bool readRecords(const std::string& path, std::vector<std::string>& recs) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    while (true) {
        std::int32_t n1, n2;
        if (!f.read(reinterpret_cast<char*>(&n1), 4)) break;
        std::string d(static_cast<std::size_t>(n1), '\0');
        f.read(&d[0], n1);
        f.read(reinterpret_cast<char*>(&n2), 4);
        if (!f || n1 != n2) return false;
        recs.push_back(std::move(d));
    }
    return true;
}

bool same(double a, double b) {   // bitwise, NaN == NaN, -0 != +0
    std::uint64_t x, y;
    std::memcpy(&x, &a, 8);
    std::memcpy(&y, &b, 8);
    if (std::isnan(a) && std::isnan(b)) return true;
    return x == y;
}

struct Counter {
    long rows{0}, diff{0};
    double maxabs{0};
    long firstRow{-1};
    double fRef{0}, fCpp{0};
};

void check(std::map<std::string, Counter>& c, const std::string& name, double cpp, double ref, long ntried) {
    Counter& k = c[name];
    ++k.rows;
    if (!same(cpp, ref)) {
        ++k.diff;
        k.maxabs = std::max(k.maxabs, std::fabs(cpp - ref));
        if (k.firstRow < 0) {
            k.firstRow = ntried;
            k.fRef = ref;
            k.fCpp = cpp;
        }
    }
}

bool report(const char* title, const std::map<std::string, Counter>& c) {
    bool ok = true;
    long tot = 0, bad = 0;
    std::printf("  %s\n  %-14s %8s %8s %12s %10s  %s\n", title, "column", "rows", "differ", "max|abs|", "1st try",
                "ref / cpp");
    for (const auto& kv : c) {
        tot += kv.second.rows;
        bad += kv.second.diff;
        if (kv.second.diff) {
            ok = false;
            std::printf("  %-14s %8ld %8ld %12.3e %10ld  %.17g / %.17g\n", kv.first.c_str(), kv.second.rows,
                        kv.second.diff, kv.second.maxabs, kv.second.firstRow, kv.second.fRef, kv.second.fCpp);
        }
    }
    std::printf("  %zu columns, %ld values compared, %ld differ%s\n", c.size(), tot, bad,
                ok ? " (all bitwise)" : "");
    return ok;
}

int stopIndexHms(const std::string& s) {
    if (s.empty()) return 0;
    for (int i = 1; i <= hms::kNumCounters; ++i)
        if (s == hms::counterName(i)) return i;
    return -1;
}
int stopIndexShms(const std::string& s) {
    if (s.empty()) return 0;
    for (int i = 1; i <= shms::kNumSlots; ++i)
        if (s == shms::counterName(i)) return i;
    return -1;
}

// ------------------------------------------------------------------------------------------------ replay
int replay(const std::string& capdir, const std::string& rundir, const std::string& cfg, const std::string& fdir,
           const std::string& runName) {
    const std::string inp = rundir + "/" + runName + ".inp";
    const std::string bin = rundir + "/" + runName + ".bin";
    if (!exists(inp) || !exists(capdir + "/cap_final.txt")) {
        std::printf("replay %s: SKIP (reference missing)\n", cfg.c_str());
        return 77;
    }
    // ---- setup
    SimcParameters p;
    InpReader rd(p);
    if (!rd.loadRun(inp)) {
        std::printf("cannot read %s\n", inp.c_str());
        return 1;
    }
    SetupOptions opt;
    opt.fortran_compat = true;   // D-031: Fortran references -> exact Fortran behaviour
    opt.run_name = runName;
    opt.data_dir = fdir;
    RunSetup setup = BuildRunSetup(p, opt);
    if (!setup.ok()) {
        std::printf("setup stopped: %s\n", setup.stops.front().c_str());
        return 1;
    }
    material::HallCCryoTarget target;
    std::ostringstream tlog;
    target.setLog(&tlog);
    target.configure(Json{{"Mh", setup.Mh}}, p);
    target.init(setup);
    driver::MonteCarlo mc(setup, target);
    std::ostringstream mlog;
    mc.setLog(&mlog);

    Tab fin, gen, rng, mcr, rec, ev;
    if (!loadTab(capdir + "/cap_final.txt", fin) || !loadTab(capdir + "/cap_generate.txt", gen) ||
        !loadTab(capdir + "/cap_rng.txt", rng) || !loadTab(capdir + "/cap_montecarlo.txt", mcr) ||
        !loadTab(capdir + "/cap_recon.txt", rec) || !loadTab(capdir + "/cap_events.txt", ev)) {
        std::printf("replay %s: SKIP (cannot read capture files; xz needed for .xz)\n", cfg.c_str());
        return 77;
    }
    std::printf("replay %s: %zu generate, %zu montecarlo calls, %zu rng values, %zu recon, %zu events\n", cfg.c_str(),
                gen.rows.size(), mcr.rows.size(), rng.rows.size(), rec.rows.size(), ev.rows.size());

    // ---- setup vs cap_final (the state montecarlo reads)
    std::map<std::string, Counter> cs;
    {
        const auto& F = fin.rows.at(0);
        auto f = [&](const char* n) { return F[fin.col(n)]; };
        check(cs, "spec_e_P", setup.spec.e.P, f("spec_e_P"), 0);
        check(cs, "spec_e_theta", setup.spec.e.theta, f("spec_e_theta"), 0);
        check(cs, "spec_e_phi", setup.spec.e.phi, f("spec_e_phi"), 0);
        check(cs, "spec_e_cos_th", setup.spec.e.cos_th, f("spec_e_cos_th"), 0);
        check(cs, "spec_e_sin_th", setup.spec.e.sin_th, f("spec_e_sin_th"), 0);
        check(cs, "spec_p_P", setup.spec.p.P, f("spec_p_P"), 0);
        check(cs, "spec_p_theta", setup.spec.p.theta, f("spec_p_theta"), 0);
        check(cs, "spec_p_phi", setup.spec.p.phi, f("spec_p_phi"), 0);
        check(cs, "spec_p_cos_th", setup.spec.p.cos_th, f("spec_p_cos_th"), 0);
        check(cs, "spec_p_sin_th", setup.spec.p.sin_th, f("spec_p_sin_th"), 0);
        check(cs, "Mh2", setup.Mh2, f("Mh2"), 0);
        check(cs, "Ebeam", setup.Ebeam, f("Ebeam"), 0);
    }
    bool okSetup = report("setup vs cap_final:", cs);

    const bool pion = setup.flags.doing_pion;   // P9.2: pion captures (decay in the hadron arm)
    const bool hadronIsHms = p.hadron_arm == 1;
    // ---- indices
    std::map<long, std::size_t> genRow, recRow;
    for (std::size_t i = 0; i < gen.rows.size(); ++i) genRow[std::lround(gen.rows[i][gen.col("ntried")])] = i;
    for (std::size_t i = 0; i < rec.rows.size(); ++i) recRow[std::lround(rec.rows[i][rec.col("ntried")])] = i;
    std::map<long, std::pair<std::vector<double>, std::vector<double>>> stream;   // ntried -> (gauss, unif)
    {
        const int it = rng.col("ntried"), ip = rng.col("phase"), ik = rng.col("kind"), iv = rng.col("value");
        for (const auto& r : rng.rows) {
            const long ph = std::lround(r[ip]);
            // montecarlo phase 2; P9.2 pion captures: also 7/8 = the decay draws of project/transp (in call order)
            if (ph != 2 && !(pion && (ph == 7 || ph == 8))) continue;
            auto& s = stream[std::lround(r[it])];
            (std::lround(r[ik]) == 2 ? s.first : s.second).push_back(r[iv]);
        }
    }

    // ---- replay
    ScriptedRandom srng(ScriptedRandom::GaussMode::Scripted);
    MainEvent main;
    Event recon;
    auto* shmsArm = dynamic_cast<shms::ShmsSpectrometer*>(&mc.electronArm());
    if (!shmsArm) shmsArm = dynamic_cast<shms::ShmsSpectrometer*>(&mc.hadronArm());
    std::map<std::string, Counter> c;
    std::map<long, std::pair<MainEvent, Event>> atSuccess;   // ntried -> state after montecarlo (for the .bin)
    long nOk = 0, nFailP = 0, nFailE = 0, nRngErr = 0;
    const char* SP[] = {"SP_e_delta", "SP_e_yptar", "SP_e_xptar", "SP_e_z", "SP_p_delta", "SP_p_yptar", "SP_p_xptar", "SP_p_z"};
    const char* RE[] = {"RE_e_delta", "RE_e_yptar", "RE_e_xptar", "RE_e_z", "RE_p_delta", "RE_p_yptar", "RE_p_xptar", "RE_p_z"};
    for (const auto& M : mcr.rows) {
        auto m = [&](const char* n) { return M[mcr.col(n)]; };
        const long nt = std::lround(m("ntried"));
        const auto& G = gen.rows.at(genRow.at(nt));
        auto g = [&](const char* n) { return G[gen.col(n)]; };
        Event orig;
        orig.Ein = g("o_Ein");
        orig.e.xptar = g("o_e_xptar");
        orig.e.yptar = g("o_e_yptar");
        orig.e.E = g("o_e_E");
        orig.e.P = g("o_e_P");
        orig.e.delta = g("o_e_delta");
        orig.p.xptar = g("o_p_xptar");
        orig.p.yptar = g("o_p_yptar");
        orig.p.E = g("o_p_E");
        orig.p.P = g("o_p_P");
        orig.p.delta = g("o_p_delta");
        main.target.x = g("t_x");
        main.target.y = g("t_y");
        main.target.z = g("t_z");
        main.target.rasterx = g("t_rasterx");
        main.target.rastery = g("t_rastery");
        main.target.Eloss = {g("t_Eloss_1"), g("t_Eloss_2"), g("t_Eloss_3")};
        main.target.teff = {g("t_teff_1"), g("t_teff_2"), g("t_teff_3")};
        main.target.Coulomb = g("t_Coulomb");
        main.ntup.resfac = 0.0;   // simc.f:178
        main.ntup.decdist = 0.0;  // simc.f:177
        main.ntup.Mh2_final = setup.Mh2;
        auto s = stream[nt];
        srng.load(s.first, s.second);

        const bool ok = mc.run(orig, main, recon, srng);
        const auto& I = mc.info();
        if (srng.error()) ++nRngErr;
        ok ? ++nOk : (I.fail_arm == 1 ? ++nFailP : ++nFailE);

        check(c, "mc_ok", ok ? 1. : 0., m("mc_ok"), nt);
        // fail_arm of the capture: 1 = the HMS failed, 2 = the SHMS (the C++ 1 = hadron arm, 2 = electron arm)
        const double failSpec = (hadronIsHms || I.fail_arm == 0) ? I.fail_arm : 3 - I.fail_arm;
        check(c, "fail_arm", ok ? 0. : failSpec, ok ? 0. : m("fail_arm"), nt);
        if (hadronIsHms) {
            check(c, "hms_stop", I.p_called ? stopIndexHms(I.p_stop) : -9., m("hms_stop"), nt);
            check(c, "shms_stop", I.e_called ? stopIndexShms(I.e_stop) : -9., m("shms_stop"), nt);
        } else {   // hadron in the SHMS, electron in the HMS (pion K1 captures)
            check(c, "shms_stop", I.p_called ? stopIndexShms(I.p_stop) : -9., m("shms_stop"), nt);
            check(c, "hms_stop", I.e_called ? stopIndexHms(I.e_stop) : -9., m("hms_stop"), nt);
        }
        if (pion) {
            check(c, "decdist", main.ntup.decdist, m("decdist"), nt);
            check(c, "Mh2_final", main.ntup.Mh2_final, m("Mh2_final"), nt);
        }
        check(c, "shms_stop_id", shmsArm ? shmsArm->kernel().stopId() : 0., m("shms_stop_id"), nt);
        // n_unif of the capture counts phase 2 only; the decay draws (phases 7/8) are n_dec_unif (pion captures)
        check(c, "n_unif", static_cast<double>(srng.nUniform()), m("n_unif") + (pion ? m("n_dec_unif") : 0.), nt);
        check(c, "n_gauss", static_cast<double>(srng.nGauss()), m("n_gauss"), nt);
        const double sp[8] = {main.SP.e.delta, main.SP.e.yptar, main.SP.e.xptar, main.SP.e.z,
                              main.SP.p.delta, main.SP.p.yptar, main.SP.p.xptar, main.SP.p.z};
        const double re[8] = {main.RECON.e.delta, main.RECON.e.yptar, main.RECON.e.xptar, main.RECON.e.z,
                              main.RECON.p.delta, main.RECON.p.yptar, main.RECON.p.xptar, main.RECON.p.z};
        for (int k = 0; k < 8; ++k) {
            check(c, SP[k], sp[k], m(SP[k]), nt);
            check(c, RE[k], re[k], m(RE[k]), nt);
        }
        const double fp[10] = {main.FP.e.x, main.FP.e.dx, main.FP.e.y, main.FP.e.dy, main.FP.e.path,
                               main.FP.p.x, main.FP.p.dx, main.FP.p.y, main.FP.p.dy, main.FP.p.path};
        const char* FPn[] = {"FP_e_x", "FP_e_dx", "FP_e_y", "FP_e_dy", "FP_e_path",
                             "FP_p_x", "FP_p_dx", "FP_p_y", "FP_p_dy", "FP_p_path"};
        for (int k = 0; k < 10; ++k) check(c, FPn[k], fp[k], m(FPn[k]), nt);
        check(c, "resfac", main.ntup.resfac, m("resfac"), nt);
        const ArmFull* arms[2] = {&recon.e, &recon.p};
        const char* an[2] = {"rec_e_", "rec_p_"};
        for (int a = 0; a < 2; ++a) {
            const ArmFull& A = *arms[a];
            const double v[8] = {A.delta, A.xptar, A.yptar, A.z, A.theta, A.phi, A.E, A.P};
            const char* f[8] = {"delta", "xptar", "yptar", "z", "theta", "phi", "E", "P"};
            for (int k = 0; k < 8; ++k) {
                const std::string n = std::string(an[a]) + f[k];
                check(c, n, v[k], M[mcr.col(n)], nt);
            }
        }
        if (ok) {
            atSuccess[nt] = {main, recon};
            // complete_recon_ev (not ported here) may modify recon arm fields (phsp: p%P, p%E, event.f:1132-1136):
            // continue with the recorded outputs so that later stale values match.
            auto itR = recRow.find(nt);
            if (itR != recRow.end()) {
                const auto& R = rec.rows[itR->second];
                auto r = [&](const std::string& n) { return R[rec.col(n)]; };
                ArmFull* arms2[2] = {&recon.e, &recon.p};
                const char* pre[2] = {"r_e_", "r_p_"};
                for (int a = 0; a < 2; ++a) {
                    ArmFull& A = *arms2[a];
                    const std::string q = pre[a];
                    A.delta = r(q + "delta");
                    A.xptar = r(q + "xptar");
                    A.yptar = r(q + "yptar");
                    A.z = r(q + "z");
                    A.theta = r(q + "theta");
                    A.phi = r(q + "phi");
                    A.E = r(q + "E");
                    A.P = r(q + "P");
                }
            }
        }
    }
    std::printf("  montecarlo calls %zu: success %ld, hadron-arm fail %ld, electron-arm fail %ld, rng queue errors %ld\n",
                mcr.rows.size(), nOk, nFailP, nFailE, nRngErr);
    bool okMc = report("montecarlo replay (bitwise):", c) && nRngErr == 0;

    // ---- .bin
    bool okBin = true;
    std::vector<std::string> recs;
    if (pion) {
        std::printf("  .bin: the 53 pion columns are checked by the full-run replay (l1_simcrun_replay_pion_*)\n");
    } else if (!readRecords(bin, recs) || recs.empty()) {
        std::printf("  .bin: cannot read %s\n", bin.c_str());
        okBin = false;
    } else {
        std::int32_t nsz;
        std::memcpy(&nsz, recs[0].data(), 4);
        const std::size_t first = 1 + static_cast<std::size_t>(nsz);
        const std::size_t nev = (recs.size() - first) / static_cast<std::size_t>(nsz);
        std::map<std::string, Counter> cb;
        const double Mp = constants::Mp;
        const double Ebeam = setup.Ebeam;
        std::size_t k = 0;
        for (const auto& E : ev.rows) {
            auto e = [&](const char* n) { return E[ev.col(n)]; };
            const long nt = std::lround(e("ntried"));
            if (k >= nev) break;
            auto itS = atSuccess.find(nt);
            if (itS == atSuccess.end()) {
                std::printf("  .bin: event of try %ld has no successful montecarlo replay\n", nt);
                okBin = false;
                ++k;
                continue;
            }
            const MainEvent& Mn = itS->second.first;
            const Event& Rc = itS->second.second;   // recon as left by montecarlo
            const auto& R = rec.rows.at(recRow.at(nt));
            auto r = [&](const char* n) { return R[rec.col(n)]; };
            double ntu[47];
            // results_write.f:44-55 (doing_hyd_elast)
            const double sh = std::sin(Rc.e.theta / 2.);
            const double poftheta = Mp * Ebeam / (2 * Ebeam * (sh * sh) + Mp);
            const double corrsing = Rc.e.P - poftheta;
            const double uqy = r("r_uq_y"), uqz = r("r_uq_z");
            const double Pm_Heepz = -(r("r_Pmy") * uqy + r("r_Pmz") * uqz) / std::sqrt(uqy * uqy + uqz * uqz);
            const double Pm_Heepy = (r("r_Pmz") * uqy - r("r_Pmy") * uqz) / std::sqrt(uqy * uqy + uqz * uqz);
            const double Pm_Heepx = -r("r_Pmx");
            // electron_arm = 5 (left): l.89-113
            ntu[1] = Rc.p.delta;
            ntu[2] = Rc.p.yptar;
            ntu[3] = Rc.p.xptar;
            ntu[4] = Rc.p.z;
            ntu[5] = Mn.FP.p.x;
            ntu[6] = Mn.FP.p.dx;
            ntu[7] = Mn.FP.p.y;
            ntu[8] = Mn.FP.p.dy;
            ntu[9] = e("o_p_delta");
            ntu[10] = e("o_p_yptar");
            ntu[11] = e("o_p_xptar");
            ntu[12] = e("m_t_z") * setup.spec.p.sin_th;
            ntu[13] = Rc.e.delta;
            ntu[14] = Rc.e.yptar;
            ntu[15] = Rc.e.xptar;
            ntu[16] = Rc.e.z;
            ntu[17] = Mn.FP.e.x;
            ntu[18] = Mn.FP.e.dx;
            ntu[19] = Mn.FP.e.y;
            ntu[20] = Mn.FP.e.dy;
            ntu[21] = e("o_e_delta");
            ntu[22] = e("o_e_yptar");
            ntu[23] = e("o_e_xptar");
            ntu[24] = -e("m_t_z") * setup.spec.e.sin_th;
            ntu[25] = r("r_q") / 1000.;
            ntu[26] = r("r_nu") / 1000.;
            ntu[27] = r("r_Q2") / 1.e6;
            ntu[28] = r("r_W") / 1000.;
            ntu[29] = r("r_epsilon");
            ntu[30] = r("r_Em") / 1000.;
            ntu[31] = r("r_Pm") / 1000.;
            ntu[32] = r("r_theta_pq");
            ntu[33] = r("r_phi_pq");
            ntu[34] = corrsing / 1000.;
            ntu[35] = Pm_Heepx / 1000.;
            ntu[36] = Pm_Heepy / 1000.;
            ntu[37] = Pm_Heepz / 1000.;
            ntu[38] = r("r_PmPar") / 1000.;
            ntu[39] = r("r_PmPer") / 1000.;
            ntu[40] = r("r_PmOop") / 1000.;
            ntu[41] = -e("m_t_rastery");
            ntu[42] = e("radphot") / 1000.;
            ntu[43] = e("m_sigcc");
            ntu[44] = e("m_weight");
            ntu[45] = Rc.e.theta;
            ntu[46] = Rc.p.theta;
            for (int i = 1; i <= 46 && i <= nsz; ++i) {
                const std::string& d = recs[first + k * static_cast<std::size_t>(nsz) + static_cast<std::size_t>(i - 1)];
                double v;
                std::memcpy(&v, d.data(), 8);
                char nm[16];
                std::snprintf(nm, sizeof nm, "ntu%02d", i);
                check(cb, nm, ntu[i], v, nt);
            }
            ++k;
        }
        std::printf("  .bin %s: NtupleSize %d, %zu events in file, %zu compared\n", bin.c_str(), nsz, nev, k);
        if (k != nev || nsz != 46) okBin = false;
        okBin = report(".bin columns from the C++ montecarlo + recorded fields (bitwise):", cb) && okBin;
    }
    if (!mlog.str().empty()) std::printf("  messages:\n%s", mlog.str().c_str());
    const bool ok = okSetup && okMc && okBin;
    std::printf("replay %s: %s\n", cfg.c_str(), ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 6 || std::string(argv[1]) != "replay") {
        std::fprintf(stderr, "usage: %s replay <cap_dir> <run_dir> <cfg> <fortran_dir>\n", argv[0]);
        return 2;
    }
    try {
        // optional argv[6]: the run name (default p4f1_<cfg>_cap; P8.2: p8f1_d_rad_cap, p8f1_fe56_rad_cap)
        return replay(argv[2], argv[3], argv[4], argv[5],
                      argc > 6 ? std::string(argv[6]) : "p4f1_" + std::string(argv[4]) + "_cap");
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
}
