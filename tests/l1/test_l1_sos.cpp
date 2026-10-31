// tests/l1/test_l1_sos.cpp
//
// Tests of the ISpectrometer "sos" (task P14.1: sos/mc_sos.f, mc_sos_hut.f, mc_sos_recon.f).
//   mc <table_dir> [rel_tol] : a P14.F1 ray table in the drv_mc.f format of reference/functions/mc_hms (README there):
//                              input.txt = spectr(=2) p_spec th_spec dpp x y z dxdz dydz m2 fry ms_flag wcs_flag
//                              decay_flag [using_coll] [ctau] u1.. [g1..]; output.txt = the inputs + the output columns.
//                              Every output column of the reference that this adapter knows is compared (real*8:
//                              |rel| <= rel_tol; real*4 focal plane, codes, RNG counts: exact), plus the stop reasons
//                              row by row and the COMMON /spectrometer2/ counters summed over the table.
//                              stop_code = index of the sSTOP_* counter in struct_sos.inc order (2 slit_hor ... 17 cal).
//   selftest                 : no Fortran table needed: maps, RNG consumption per stop (no smearing / scripted zero
//                              gaussians), counter bookkeeping on a flat ray scan, recon sanity, the mc_sos_recon zero
//                              guard, mock sieve, decay, the plug-in adapter (raster fry), the default arm choice.
// Exit: 0 pass, 1 fail, 77 table missing.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/MonteCarlo.h"
#include "simc/optics/CosyEngine.h"
#include "simc/spectrometers/sos/SosSpectrometer.h"
#include "simc/spectrometers/sos/kernel/McSos.h"

using namespace simc;
using simc::sos::McSos;
using simc::sos::McSosArgs;
using simc::sos::McSosDetail;
using simc::sos::SosSpectrometer;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

constexpr double kUnset = -999999.;   // drv_mc.f preset of x_fp..dy_fp, resmult
constexpr double kMe2 = 0.51099906 * 0.51099906;
constexpr double kMp2 = 938.27231 * 938.27231;
constexpr double kMpi = 139.57018;
constexpr double kMmu = 105.6583755;   // constants.inc:25 (project.f:77 decay to a muon)

// ---------------------------------------------------------------------------------------------- mc (P14.F1 tables)
int runMc(const std::string& dir, double tol) {
    simcl1::Table ref;
    std::string err;
    if (!ref.load(dir + "/output.txt", err)) {
        std::printf("mc_sos (%s): SKIP (%s)\n", dir.c_str(), err.c_str());
        return 77;
    }
    SosSpectrometer plugin;
    McSos& sos = plugin.kernel();
    sos.setFortranCompat(true);   // D-031: Fortran tables -> exact Fortran behaviour
    std::ostringstream log;
    sos.setLog(&log);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<int> codes;

    // the output columns this adapter can produce, compared if present in the reference
    std::vector<ColumnSpec> known;
    for (const char* n : {"ok", "stop_code", "reached_hut"}) known.push_back({n, 0.0, 0.0});
    for (const char* n : {"x_fp", "dx_fp", "y_fp", "dy_fp"}) known.push_back({n, 0.0, 0.0});   // real*4: exact
    for (const char* n : {"dpp_out", "dxdz_out", "dydz_out", "y_out", "pathlen", "m2_out"}) known.push_back({n, tol, 0.0});
    known.push_back({"resmult", 0.0, 0.0});
    for (const char* n : {"xs_end", "ys_end", "dxdzs_end", "dydzs_end", "dpps_end"}) known.push_back({n, tol, 0.0});
    for (const char* n : {"decdist", "Mh2_final"}) known.push_back({n, tol, 0.0});
    for (const char* n : {"n_gauss", "n_unif", "rng_err", "nsig_last", "row"}) known.push_back({n, 0.0, 0.0});
    std::vector<ColumnSpec> c;
    for (const auto& k : known)
        if (ref.col(k.name) >= 0) c.push_back(k);
    std::printf("  comparing %zu output columns\n", c.size());

    int rc = simcl1::run("mc_sos (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        std::vector<double> g, u;
        for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
        for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
        rng.load(g, u);
        if (std::lround(r["spectr"]) != 2) throw std::runtime_error("not an SOS row");
        McSosArgs a;
        a.p_spec = r["p_spec"];
        a.th_spec = r["th_spec"];
        a.dpp = r["dpp"];
        a.x = r["x"];
        a.y = r["y"];
        a.z = r["z"];
        a.dxdz = r["dxdz"];
        a.dydz = r["dydz"];
        a.m2 = r["m2"];
        a.fry = r["fry"];
        a.ms_flag = std::lround(r["ms_flag"]) != 0;
        a.wcs_flag = std::lround(r["wcs_flag"]) != 0;
        a.decay_flag = std::lround(r["decay_flag"]) != 0;
        if (r.has("ctau")) a.ctau = r["ctau"];
        a.decdist = r.has("decdist_in") ? r["decdist_in"] : 0.;
        a.Mh2_final = r.has("Mh2_final_in") ? r["Mh2_final_in"] : a.m2;
        a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
        a.resmult = kUnset;
        a.pathlen = 0.;
        const McSosDetail det = sos.mcSos(a, rng);
        codes.push_back(det.stop_code);
        o["ok"] = a.ok_spec ? 1. : 0.;
        o["stop_code"] = det.stop_code;
        o["reached_hut"] = det.reached_hut ? 1. : 0.;
        o["x_fp"] = a.x_fp;
        o["dx_fp"] = a.dx_fp;
        o["y_fp"] = a.y_fp;
        o["dy_fp"] = a.dy_fp;
        o["dpp_out"] = a.dpp;
        o["dxdz_out"] = a.dxdz;
        o["dydz_out"] = a.dydz;
        o["y_out"] = a.y;
        o["pathlen"] = a.pathlen;
        o["m2_out"] = a.m2;
        o["resmult"] = a.resmult;
        const optics::Track& t = sos.track();
        o["xs_end"] = t.xs;
        o["ys_end"] = t.ys;
        o["dxdzs_end"] = t.dxdzs;
        o["dydzs_end"] = t.dydzs;
        o["dpps_end"] = t.dpps;
        o["decdist"] = a.decdist;
        o["Mh2_final"] = a.Mh2_final;
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        o["nsig_last"] = rng.lastNsig();
        o["row"] = static_cast<double>(r.number());
    });
    if (rc == 77) return rc;

    // Per-reason agreement and the COMMON /spectrometer2/ counters summed over the table.
    bool okc = true;
    const int ic = ref.col("stop_code"), ih = ref.col("reached_hut");
    if (ic >= 0) {
        std::map<int, long> nref, ncpp, nagree;
        long hutRef = 0;
        for (std::size_t i = 0; i < ref.rows.size(); ++i) {
            const int cr = static_cast<int>(std::lround(ref.rows[i][ic]));
            ++nref[cr];
            if (ih >= 0) hutRef += std::lround(ref.rows[i][ih]);
            if (i < codes.size()) {
                ++ncpp[codes[i]];
                if (codes[i] == cr) ++nagree[cr];
            }
        }
        std::printf("  stop reason agreement (row by row):\n  %4s %-10s %8s %8s %8s\n", "code", "name", "ref", "cpp",
                    "agree");
        for (const auto& kv : nref) {
            const char* nm = kv.first == 0 ? "success" : sos::counterName(kv.first);
            std::printf("  %4d %-10s %8ld %8ld %8ld\n", kv.first, nm, kv.second, ncpp[kv.first], nagree[kv.first]);
            if (nagree[kv.first] != kv.second || ncpp[kv.first] != kv.second) okc = false;
        }
        const sos::SosCounters& cn = sos.counters();
        std::printf("  counters: trials %ld (rows %zu), hut %ld (ref %ld), successes %ld (ref %ld)\n", cn[sos::kTrials],
                    ref.rows.size(), cn[sos::kHut], hutRef, cn[sos::kSuccesses], nref[0]);
        if (cn[sos::kTrials] != static_cast<long>(ref.rows.size()) || (ih >= 0 && cn[sos::kHut] != hutRef) ||
            cn[sos::kSuccesses] != nref[0])
            okc = false;
        for (int k = sos::kSlitHor; k <= sos::kCal; ++k) {
            if (k == sos::kHut) continue;
            if (cn[k] != nref[k]) {
                std::printf("  counter %s = %ld, reference rows with that stop code = %ld\n", sos::counterName(k), cn[k],
                            nref[k]);
                okc = false;
            }
        }
    }
    if (!log.str().empty()) {
        std::printf("  unexpected messages:\n%s", log.str().c_str());
        okc = false;
    }
    std::printf("  stop reasons / counters: %s\n", okc ? "PASS" : "FAIL");
    return (rc == 0 && okc) ? 0 : 1;
}

// ---------------------------------------------------------------------------------------------- selftest
int nfail = 0;
void check(bool ok, const std::string& what) {
    std::printf("  %-96s %s\n", what.c_str(), ok ? "ok" : "FAIL");
    if (!ok) ++nfail;
}

McSosArgs ray(double p_spec, double dpp, double x, double y, double dxdz, double dydz, double m2, bool smear = false) {
    McSosArgs a;
    a.p_spec = p_spec;
    a.th_spec = 0.5;
    a.dpp = dpp;
    a.x = x;
    a.y = y;
    a.z = 0.;
    a.dxdz = dxdz;
    a.dydz = dydz;
    a.m2 = m2;
    a.fry = x;
    a.ms_flag = smear;
    a.wcs_flag = smear;
    a.decay_flag = false;
    a.Mh2_final = m2;
    a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
    a.resmult = kUnset;
    return a;
}

int runSelftest() {
    // 1. maps
    {
        McSos s;
        std::ostringstream log;
        s.setLog(&log);
        bool loaded = true;
        try {
            s.loadMaps();
        } catch (const std::exception& e) {
            loaded = false;
            std::printf("  %s\n", e.what());
        }
        check(loaded, "SOS COSY maps load (" + McSos::defaultMatrixDir() + "/sos)");
        if (!loaded) return 1;
        optics::CosyEngine& e = s.engine();
        check(e.nClasses(2) == 10, "sos/forward_cosy.dat: 10 classes (mc_sos.f:158)");
        check(e.adrift(2, 1) && e.driftdist(2, 1) == 150., "class 1 is a 150 cm drift (mc_sos.f:202-203)");
        check(e.reconTerms(2).size() == 461, "sos/recon_cosy.dat: 461 terms (got " +
                                                 std::to_string(e.reconTerms(2).size()) + ")");
        check(log.str().empty(), "no transp_init length warnings");
    }

    // 2. RNG consumption: no smearing -> 1 grnd for a ray reaching the hut, 0 gauss1; nothing before the hut
    McSos s;
    std::ostringstream log;
    s.setLog(&log);
    s.setFortranCompat(true);
    ScriptedRandom rng;
    {
        McSosArgs a = ray(1000., 1.0, 0.1, -0.05, 0.002, -0.003, kMe2);
        rng.load({}, {0.5});
        const McSosDetail d = s.mcSos(a, rng);
        check(a.ok_spec && d.stop_code == 0 && d.reached_hut && d.fit_done, "near-central electron passes");
        check(rng.nUniform() == 1 && rng.nGauss() == 0 && !rng.error(), "no smearing: 1 grnd (resmult), 0 gauss1");
        check(a.resmult == 1.0, "u = 0.5 -> resmult 1");
        check(std::fabs(a.dpp - 1.0) < 0.02 && std::fabs(a.dxdz - 0.002) < 5e-4 && std::fabs(a.dydz + 0.003) < 5e-4,
              "recon close to the generated ray (delta " + std::to_string(a.dpp) + ", xptar " + std::to_string(a.dxdz) +
                  ", yptar " + std::to_string(a.dydz) + ")");
        std::printf("    pathlen = %.6f cm, x_fp %.5f dx_fp %.6f y_fp %.5f dy_fp %.6f\n", a.pathlen, a.x_fp, a.dx_fp,
                    a.y_fp, a.dy_fp);
        check(a.pathlen > 1000. && a.pathlen < 1300., "pathlen to the back of the calorimeter");
        // the same ray with ms = wcs = on and all gaussians 0: same result bitwise, exactly 180 gauss1(99)
        McSosArgs b = ray(1000., 1.0, 0.1, -0.05, 0.002, -0.003, kMe2, true);
        rng.load(std::vector<double>(400, 0.0), {0.1});
        const McSosDetail db = s.mcSos(b, rng);
        check(b.ok_spec && db.stop_code == 0, "smeared (zero gaussians) passes");
        check(rng.nGauss() == 180 && rng.nUniform() == 1 && rng.lastNsig() == 99.0,
              "ms=wcs=on success: 180 gauss1(99) + 1 grnd (got " + std::to_string(rng.nGauss()) + ")");
        check(b.resmult == 2.0, "u = 0.1 -> resmult 2");
        check(b.dpp == a.dpp && b.dxdz == a.dxdz && b.dydz == a.dydz && b.y == a.y && b.x_fp == a.x_fp &&
                  b.dy_fp == a.dy_fp && b.pathlen == a.pathlen,
              "zero gaussians reproduce the unsmeared ray bitwise");
    }
    // a ray stopped at the slit: no random numbers
    {
        McSosArgs a = ray(1000., 0., 0., 0., 0.0, 0.2, kMp2);
        rng.load({}, {});
        const McSosDetail d = s.mcSos(a, rng);
        check(!a.ok_spec && d.stop_code == sos::kSlitHor && rng.nUniform() == 0 && rng.nGauss() == 0,
              "dydz = 0.2 stops at slit_hor, no random numbers");
        check(a.x_fp == kUnset && a.resmult == kUnset && a.dpp == 0. && a.dydz == 0.2,
              "outputs not written for a stopped ray");
        McSosArgs b = ray(1000., 0., 0., 0., 0.1, 0.0, kMp2);
        const McSosDetail db = s.mcSos(b, rng);
        check(db.stop_code == sos::kSlitVert, "dxdz = 0.1 stops at slit_vert");
        McSosArgs c = ray(1000., 0., 0., 0., 0.033, 0.05, kMp2);
        const McSosDetail dc = s.mcSos(c, rng);
        check(dc.stop_code == sos::kSlitOct, "corner ray stops at slit_oct (got " + std::to_string(dc.stop_code) + ")");
    }

    // 3. flat scan: counter bookkeeping, every stop code consistent with the counters, recon sanity
    {
        McSos k;
        k.setLog(&log);
        k.setFortranCompat(true);
        Mt64Random mt(20261030);
        std::map<int, long> tally;
        double s2d = 0, s2x = 0, s2y = 0;
        long nres = 0;
        const long N = 200000;
        for (long i = 0; i < N; ++i) {
            const bool point = (i % 2) == 0;   // half of the rays from a point target (recon sanity)
            const double dpp = -30. + 60. * mt.uniform();
            const double x = point ? 0. : -1. + 2. * mt.uniform();
            const double y = point ? 0. : -2. + 4. * mt.uniform();
            const double dxdz = -0.07 + 0.14 * mt.uniform();
            const double dydz = -0.09 + 0.18 * mt.uniform();
            McSosArgs a = ray(1200., dpp, x, y, dxdz, dydz, (i % 3) ? kMp2 : kMe2);
            a.fry = point ? 0. : x;
            const McSosDetail d = k.mcSos(a, mt);
            ++tally[d.stop_code];
            if (d.stop_code == 0 && point) {
                s2d += (a.dpp - dpp) * (a.dpp - dpp);
                s2x += (a.dxdz - dxdz) * (a.dxdz - dxdz);
                s2y += (a.dydz - dydz) * (a.dydz - dydz);
                ++nres;
            }
        }
        const sos::SosCounters& c = k.counters();
        std::printf("    %-10s %8s %8s\n", "counter", "value", "rays");
        long sumStops = 0;
        bool same = true;
        for (int i = 1; i <= sos::kNumCounters; ++i) {
            const long rays = i == sos::kSuccesses ? tally[0] : (i == sos::kTrials || i == sos::kHut) ? -1 : tally[i];
            std::printf("    %-10s %8ld %8ld\n", sos::counterName(i), c[i], rays);
            if (rays >= 0 && rays != c[i]) same = false;
            if (i != sos::kTrials && i != sos::kHut && i != sos::kSuccesses) sumStops += c[i];
        }
        check(same, "every stop code matches the counter it incremented");
        check(c[sos::kTrials] == N && c[sos::kTrials] == sumStops + c[sos::kSuccesses], "trials = stops + successes");
        check(c[sos::kHut] == c[sos::kDc1] + c[sos::kDc2] + c[sos::kScin] + c[sos::kCal] + c[sos::kSuccesses],
              "hut = dc1 + dc2 + scin + cal + successes");
        check(c[sos::kCal] == 0, "sSTOP_cal never incremented (mc_sos_hut.f:519-525 commented out)");
        int nhit = 0;
        for (int i : {sos::kSlitHor, sos::kSlitVert, sos::kSlitOct, sos::kQuadIn, sos::kQuadMid, sos::kQuadOut,
                      sos::kBm01In, sos::kBm01Out, sos::kBm02In, sos::kBm02Out, sos::kExit, sos::kDc1, sos::kDc2,
                      sos::kScin})
            if (c[i] > 0) ++nhit;
        std::printf("    %d of 14 stop reasons reached by the flat scan\n", nhit);
        check(c[sos::kSuccesses] > 1000 && c[sos::kSlitHor] > 0 && c[sos::kSlitVert] > 0 && c[sos::kSlitOct] > 0,
              "successes and the three slit stops are populated");
        const double rd = std::sqrt(s2d / nres), rx = std::sqrt(s2x / nres), ry = std::sqrt(s2y / nres);
        std::printf("    point-target recon rms (no smearing, %ld rays): delta %.4g %%, xptar %.4g, yptar %.4g\n", nres,
                    rd, rx, ry);
        check(nres > 1000 && rd < 0.1 && rx < 2e-3 && ry < 2e-3, "recon rms small (forward and recon maps agree)");
        check(log.str().empty(), "no messages");
    }

    // 4. mc_sos_recon zero guard: EVERY hut value |v| <= 1e-30 -> 1e-30 (mc_sos_recon.f:120-122)
    {
        optics::CosyEngine& e = s.engine();
        optics::Track t0, t1;
        t1.xs = 1e-28;      // xs/100 = 1e-30
        t1.dxdzs = 1e-30;
        t1.ys = 1e-28;
        t1.dydzs = 1e-30;
        double a0, b0, c0, d0, a1, b1, c1, d1;
        e.recon(2, t0, 0.0, a0, b0, c0, d0);
        e.recon(2, t1, 1e-28, a1, b1, c1, d1);
        check(a0 == a1 && b0 == b1 && c0 == c1 && d0 == d1, "recon(2) of a zero track == recon of the 1e-30 track");
        check(std::isfinite(a0) && std::isfinite(d0), "recon(2) of a zero track is finite");
    }

    // 5. mock sieve (mc_sos.f:122-137): 2 grnd per call before the transport
    {
        McSos k;
        k.setLog(&log);
        k.setUseSieve(true);
        McSosArgs a = ray(1000., 0.5, 0.0, 0.0, 0.004, 0.003, kMe2);
        rng.load({}, {0.25, 0.5, 0.5});
        const McSosDetail d = k.mcSos(a, rng);
        check(d.stop_code == 0 && rng.nUniform() == 3, "use_sieve: 2 grnd + 1 (hut) for a success");
    }

    // 6. decay in flight (pions, short ctau): decays happen, m2 changes to the muon mass, decdist accumulated
    {
        McSos k;
        k.setLog(&log);
        Mt64Random mt(7);
        long ndec = 0, nok = 0;
        for (int i = 0; i < 20000; ++i) {
            McSosArgs a = ray(800., -10. + 20. * mt.uniform(), 0., 0., -0.04 + 0.08 * mt.uniform(),
                              -0.06 + 0.12 * mt.uniform(), kMpi * kMpi);
            a.decay_flag = true;
            a.ctau = 780.45;
            k.mcSos(a, mt);
            if (a.ok_spec) ++nok;
            if (a.m2 != kMpi * kMpi) {
                ++ndec;
                if (std::fabs(std::sqrt(a.m2) - kMmu) > 0.01 || a.Mh2_final != a.m2) {
                    check(false, "decay product is a muon with Mh2_final = m2");
                    break;
                }
            }
        }
        std::printf("    pions: %ld decays, %ld successes of 20000\n", ndec, nok);
        check(ndec > 100 && nok > 100, "pion decays in flight occur in the SOS");
    }

    // 7. the plug-in: registry, counters, raster fry (simc.f:1472) as the recon argument
    {
        auto sp = Registry<ISpectrometer>::create("sos");
        SimcParameters p;
        sp->configure(Json::object(), p);
        RunSetup setup;
        sp->init(StageContext(setup));
        SpecInput in;
        in.p_spec = 1000.;
        in.th_spec = 0.5;
        in.delta = 1.0;
        in.x = 0.1;
        in.y = -0.05;
        in.dxdz = 0.002;
        in.dydz = -0.003;
        in.m2 = kMe2;
        in.fry = 0.1;          // xtar_init: NOT used by the SOS
        in.raster_fry = 0.2;   // simc.f fry
        Mt64Random mt(3);
        const SpecResult r = sp->transport(in, mt);
        McSos k;
        k.setLog(nullptr);
        McSosArgs a = ray(1000., 1.0, 0.1, -0.05, 0.002, -0.003, kMe2);
        a.fry = 0.2;
        Mt64Random mt2(3);
        k.mcSos(a, mt2);
        check(sp->name() == "sos" && r.ok && r.resolution_factor.has_value(), "plug-in \"sos\": transport ok, resfac set");
        check(r.delta_recon == a.dpp && r.dxdz_recon == a.dxdz && r.dydz_recon == a.dydz && r.y_recon == a.y &&
                  r.x_fp == a.x_fp && r.pathlen == a.pathlen,
              "plug-in passes SpecInput::raster_fry as the recon fry (kernel with fry = 0.2, bitwise)");
        const auto& m = sp->stopCounters();
        check(m.size() == 18 && m.at("trials") == 1 && m.at("successes") == 1 && m.count("bm02_out") == 1,
              "stopCounters: 18 sSTOP_* keys");
    }

    // 8. the default arm of code 2 is the SOS plug-in
    {
        bool ok = false;
        try {
            ok = driver::MonteCarlo::DefaultArm(2).type == "sos";
        } catch (const std::exception& e) {
            std::printf("  %s\n", e.what());
        }
        check(ok, "MonteCarlo::DefaultArm(2) = \"sos\"");
        bool np = false;
        try {
            driver::MonteCarlo::DefaultArm(3);
        } catch (const NotPorted&) {
            np = true;
        }
        check(np, "MonteCarlo::DefaultArm(3) (HRS-R) still NotPorted (P14.2)");
    }

    std::printf("selftest: %s (%d failures)\n", nfail ? "FAIL" : "PASS", nfail);
    return nfail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s mc <table_dir> [rel_tol] | selftest\n", argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    try {
        if (mode == "selftest") return runSelftest();
        if (mode == "mc" && argc > 2) return runMc(argv[2], argc > 3 ? std::atof(argv[3]) : 1e-12);
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
