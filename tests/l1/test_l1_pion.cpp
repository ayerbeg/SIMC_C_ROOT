// tests/l1/test_l1_pion.cpp
//
// L1 tests of the exclusive-pion physics (task P9.1) against the P9.F1 Fortran tables (reference/functions/, see
// reference/P9.F1_report.md). fortran_compat mode (D-031); the tolerance argument is |rel| (1e-12 in pion.cmake); the
// runner reports the bitwise-equal rows per column, NaN/Inf rows must match (NaN == NaN).
//
//   files <data/pion> <fortran dir>                     md5 of the copies (MD5SUMS) + byte equality with the Fortran
//   sig_blok|sig_param04|sig_param_3000 <dir> <tol>     physics_pion.f:404-573
//   sig_param_2021 <dir> <tol>                          physics_pion.f:738-805 (+ COMMON /excmn/)
//   exclfit <dir> <tol>                                 physics_pion.f:807-854
//   sigmaid <dir> <data/pion> <tol>                     physics_pion.f:577-736 (+ the nearest-bin indices)
//   sig_multipole <dir> <data/pion> <tol>               physics_pion.f:197-398, rows in order (SAVEd state)
//   transform_to_cm <dir> <tol>                         jacobians.f:1-282 on the pion table (reactions::transformToCm)
//   peepi <dir> <pion_state dir> <data/pion> <tol>      physics_pion.f:1-194 with the per-configuration state
//   pizero_decay <dir> <tol>                            pizero_decay.f with ScriptedRandom (u1, u2; u3 detects extras)
//   plugin <peepi dir> <pion_state dir> <data/pion>     the IPhysicsModel "pion" from the Registry on the peepi rows
//   fortran_diff_D038_F10 <peepi dir> <pion_state dir> <data/pion>   stale ntup%sigcm2, both modes
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

#include "L1Runner.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/physics/pion/PionCrossSections.h"
#include "simc/physics/pion/PionModel.h"
#include "simc/physics/pion/PizeroDecay.h"
#include "simc/reactions/EventChain.h"
#include "simc/stages/IPhysicsModel.h"

using namespace simc;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

int iv(const Row& r, const char* c) { return static_cast<int>(std::lround(r[c])); }
bool same(double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0; }

std::vector<ColumnSpec> specs(const std::vector<std::string>& names, double tol) {
    std::vector<ColumnSpec> v;
    for (const auto& n : names) v.push_back(ColumnSpec{n, tol, 0.0});
    return v;
}

// ------------------------------------------------------------------------------------------------- md5 (RFC 1321)
std::string md5hex(const std::string& msg) {
    static const std::uint32_t K[64] = {
        0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
        0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
        0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
        0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
        0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
        0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
        0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
        0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391};
    static const int R[64] = {7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9,  14, 20, 5, 9,
                              14, 20, 5, 9,  14, 20, 5, 9,  14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
                              4,  11, 16, 23, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
    std::string m = msg;
    const std::uint64_t bits = static_cast<std::uint64_t>(msg.size()) * 8u;
    m += static_cast<char>(0x80);
    while (m.size() % 64 != 56) m += '\0';
    for (int i = 0; i < 8; ++i) m += static_cast<char>((bits >> (8 * i)) & 0xff);
    std::uint32_t h0 = 0x67452301, h1 = 0xefcdab89, h2 = 0x98badcfe, h3 = 0x10325476;
    for (std::size_t off = 0; off < m.size(); off += 64) {
        std::uint32_t w[16];
        for (int i = 0; i < 16; ++i) {
            const unsigned char* p = reinterpret_cast<const unsigned char*>(m.data() + off + 4 * i);
            w[i] = p[0] | (p[1] << 8) | (p[2] << 16) | (static_cast<std::uint32_t>(p[3]) << 24);
        }
        std::uint32_t a = h0, b = h1, c = h2, d = h3;
        for (int i = 0; i < 64; ++i) {
            std::uint32_t f;
            int g;
            if (i < 16) { f = (b & c) | (~b & d); g = i; }
            else if (i < 32) { f = (d & b) | (~d & c); g = (5 * i + 1) % 16; }
            else if (i < 48) { f = b ^ c ^ d; g = (3 * i + 5) % 16; }
            else { f = c ^ (b | ~d); g = (7 * i) % 16; }
            const std::uint32_t tmp = d;
            d = c;
            c = b;
            const std::uint32_t x = a + f + K[i] + w[g];
            b = b + ((x << R[i]) | (x >> (32 - R[i])));
            a = tmp;
        }
        h0 += a; h1 += b; h2 += c; h3 += d;
    }
    char out[33];
    const std::uint32_t hs[4] = {h0, h1, h2, h3};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) std::snprintf(out + 8 * i + 2 * j, 3, "%02x", (hs[i] >> (8 * j)) & 0xff);
    return std::string(out, 32);
}

bool slurp(const std::string& p, std::string& out) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return false;
    std::ostringstream s;
    s << f.rdbuf();
    out = s.str();
    return true;
}

int runFiles(const std::string& dataDir, const std::string& fortranDir) {
    std::ifstream sums(dataDir + "/MD5SUMS");
    if (!sums) return 77;
    int bad = 0, n = 0;
    std::string md5, name;
    while (sums >> md5 >> name) {
        ++n;
        std::string a, b;
        if (!slurp(dataDir + "/" + name, a)) {
            std::printf("  missing copy %s\n", name.c_str());
            ++bad;
            continue;
        }
        const std::string got = md5hex(a);
        const bool okmd5 = got == md5;
        std::string cmp = "Fortran file missing (not compared)";
        bool okF = true;
        if (slurp(fortranDir + "/" + name, b)) {
            okF = a == b;
            cmp = okF ? "identical to the Fortran file" : "DIFFERS from the Fortran file";
        }
        std::printf("  %-16s md5 %s %s, %s\n", name.c_str(), got.c_str(), okmd5 ? "ok" : "MISMATCH", cmp.c_str());
        if (!okmd5 || !okF) ++bad;
    }
    std::printf("pion data files: %d listed, %d bad -> %s\n", n, bad, bad || n != 5 ? "FAIL" : "PASS");
    return (bad || n != 5) ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------- simple fits
int runSimpleFit(const std::string& which, const std::string& dir, double tol) {
    using Fn = double (*)(double, double, double, double, double, double, double, int);
    Fn fn = which == "sig_blok" ? pion::sigBlok : which == "sig_param04" ? pion::sigParam04 : pion::sigParam3000;
    return simcl1::run(which + " (physics_pion.f)", dir, specs({"sig"}, tol), [&](const Row& r, Outputs& o) {
        o["sig"] = fn(r["thetacm"], r["phicm"], r["t"], r["q2_gev"], r["s_gev"], r["eps"], r["mtar_gev"],
                      iv(r, "which_pion"));
    });
}

int runParam2021(const std::string& dir, double tol) {
    return simcl1::run("sig_param_2021 (physics_pion.f:738-805)", dir,
                       specs({"sig", "excmn_sigl", "excmn_sigt", "excmn_sigttv", "excmn_sigltv"}, tol),
                       [&](const Row& r, Outputs& o) {
                           pion::Excmn c;
                           o["sig"] = pion::sigParam2021(r["thcm"], r["phicm"], r["t"], r["q2"], r["wsq"], r["eps"],
                                                         iv(r, "which_pion"), iv(r, "doing_pizero") != 0, &c);
                           o["excmn_sigl"] = c.sigl;
                           o["excmn_sigt"] = c.sigt;
                           o["excmn_sigttv"] = c.sigttv;
                           o["excmn_sigltv"] = c.sigltv;
                       });
}

int runExclfit(const std::string& dir, double tol) {
    return simcl1::run("exclfit (physics_pion.f:807-854)", dir, specs({"sigl", "sigt", "sigtt", "siglt", "sig"}, tol),
                       [&](const Row& r, Outputs& o) {
                           std::array<double, 17> p{};
                           for (int k = 1; k <= 17; ++k) p[static_cast<std::size_t>(k - 1)] = r["p" + std::to_string(k)];
                           const pion::ExclfitResult x = pion::exclfit(r["t"], r["thetacm"], r["phicm"], r["q2_gev"],
                                                                       r["s_gev"], r["eps"], p, r["fpifact"]);
                           o["sigl"] = x.sigl;
                           o["sigt"] = x.sigt;
                           o["sigtt"] = x.sigtt;
                           o["siglt"] = x.siglt;
                           o["sig"] = x.sig;
                       });
}

int runSigmaid(const std::string& dir, const std::string& dataDir, double tol) {
    pion::MaidTables maid;
    std::string err;
    if (!maid.load(dataDir, &err)) {
        std::printf("sigmaid: %s\n", err.c_str());
        return 1;
    }
    return simcl1::run("sigmaid (physics_pion.f:577-736)", dir, specs({"sig0", "hx_iq", "hx_iw", "hx_ith"}, tol),
                       [&](const Row& r, Outputs& o) {
                           o["sig0"] = pion::sigmaid(iv(r, "ipi"), r["q2"], r["w"], r["e0"], r["costh"], r["phi"], maid).sig0;
                           const pion::MaidBin b = pion::sigmaidBin(r["q2"], r["w"], r["costh"]);
                           o["hx_iq"] = b.iq;
                           o["hx_iw"] = b.iw;
                           o["hx_ith"] = b.ith;
                       });
}

int runMultipole(const std::string& dir, const std::string& dataDir, double tol) {
    pion::SigMultipole mp;
    std::string err;
    if (!mp.load(dataDir + "/MULTIPOL.file", &err)) {
        std::printf("sig_multipole: %s\n", err.c_str());
        return 1;
    }
    std::ostringstream log;
    const int rc = simcl1::run("sig_multipole (physics_pion.f:197-398, rows in order)", dir, specs({"sig", "nu_out"}, tol),
                               [&](const Row& r, Outputs& o) {
                                   double nu = r["nu"];
                                   o["sig"] = mp(nu, r["efer"], r["qstar"], r["tfcos"], r["ppicm"], r["s_gev"],
                                                 r["thetacm"], r["phicm"], r["eps"], iv(r, "final_state"), r["q2_gev"],
                                                 r["mtar_gev"], r["pfer"], r["mh"], &log);
                                   o["nu_out"] = nu;
                               });
    std::size_t warn = 0;
    for (std::size_t p = log.str().find("CHEESY"); p != std::string::npos; p = log.str().find("CHEESY", p + 1)) ++warn;
    std::printf("  'LOOK OUT CHEESY POOF BREATH!' warnings: %zu\n", warn);
    return rc;
}

int runTransform(const std::string& dir, double tol) {
    const std::vector<std::string> cols = {"gstar",   "bstar",   "bstarx",   "bstary",   "bstarz",  "nustar",   "qstar",
                                           "qstarx",  "qstary",  "qstarz",   "ehadcm",   "phadcm",  "phadcmx",  "phadcmy",
                                           "phadcmz", "ebeamcm", "pbeamcm",  "pbeamcmx", "pbeamcmy", "pbeamcmz", "etarcm",
                                           "ptarcm",  "ptarcmx", "ptarcmy",  "ptarcmz",  "thetacm", "phicm",    "phiqn",
                                           "jacobian", "jac_old"};
    return simcl1::run("transform_to_cm pion (jacobians.f:1-282, reactions::transformToCm)", dir, specs(cols, tol),
                       [&](const Row& r, Outputs& o) {
                           Event v;
                           MainEvent m;
                           v.Ein = r["Ein"];
                           v.nu = r["nu"];
                           v.q = r["q"];
                           v.p.P = r["p_P"];
                           v.p.E = r["p_E"];
                           v.uq = {r["uq_x"], r["uq_y"], r["uq_z"]};
                           v.up = {r["up_x"], r["up_y"], r["up_z"]};
                           m.phi_pq = r["phi_pq"];
                           const reactions::CmTransform c = reactions::transformToCm(
                               v, m, r["pfer"], r["pferx"], r["pfery"], r["pferz"], r["efer"]);
                           const double vals[] = {c.gstar,   c.bstar,    c.bstarx,   c.bstary,   c.bstarz,  c.nustar,
                                                  c.qstar,   c.qstarx,   c.qstary,   c.qstarz,   c.ehadcm,  c.phadcm,
                                                  c.phadcmx, c.phadcmy,  c.phadcmz,  c.ebeamcm,  c.pbeamcm, c.pbeamcmx,
                                                  c.pbeamcmy, c.pbeamcmz, c.etarcm,  c.ptarcm,   c.ptarcmx, c.ptarcmy,
                                                  c.ptarcmz, c.thetacm,  c.phicm,    c.phiqn,    c.jacobian, c.jac_old};
                           for (std::size_t k = 0; k < cols.size(); ++k) o[cols[k]] = vals[k];
                       });
}

// ------------------------------------------------------------------------------------------------- peepi
/// Per-configuration state from pion_state/output.txt (the Fortran start-up of each configuration).
struct CfgState {
    int which_pion{0};
    bool doing_pizero{false};
    double Mtar_struck{0};
};
bool loadStates(const std::string& stateDir, std::map<int, CfgState>& st) {
    simcl1::Table t;
    std::string err;
    if (!t.load(stateDir + "/output.txt", err)) {
        std::printf("pion_state: %s\n", err.c_str());
        return false;
    }
    for (std::size_t i = 0; i < t.rows.size(); ++i) {
        Row r(t, i);
        CfgState c;
        c.which_pion = iv(r, "which_pion");
        c.doing_pizero = iv(r, "doing_pizero") != 0;
        c.Mtar_struck = r["targ_Mtar_struck"];
        st[iv(r, "cfg")] = c;
    }
    return !st.empty();
}

void readEvent(const Row& r, const std::string& pre, Event& v) {
    auto g = [&](const char* n) { return r[pre + n]; };
    v.Ein = g("Ein"); v.Em = g("Em"); v.Pm = g("Pm"); v.Emiss = g("Emiss"); v.Pmiss = g("Pmiss");
    v.Pmx = g("Pmx"); v.Pmy = g("Pmy"); v.Pmz = g("Pmz"); v.PmPar = g("PmPar"); v.PmPer = g("PmPer");
    v.PmOop = g("PmOop"); v.nu = g("nu"); v.q = g("q"); v.Q2 = g("Q2"); v.Trec = g("Trec"); v.W = g("W");
    v.Mrec = g("Mrec"); v.epsilon = g("epsilon"); v.theta_pq = g("theta_pq"); v.theta_tarq = g("theta_tarq");
    v.phi_pq = g("phi_pq"); v.phi_targ = g("phi_targ"); v.beta = g("beta"); v.phi_s = g("phi_s"); v.phi_c = g("phi_c");
    v.zhad = g("zhad"); v.pt2 = g("pt2"); v.xbj = g("xbj");
    auto arm = [&](const std::string& a, ArmFull& x) {
        x.delta = r[pre + a + "_delta"]; x.xptar = r[pre + a + "_xptar"]; x.yptar = r[pre + a + "_yptar"];
        x.z = r[pre + a + "_z"]; x.theta = r[pre + a + "_theta"]; x.phi = r[pre + a + "_phi"];
        x.E = r[pre + a + "_E"]; x.P = r[pre + a + "_P"];
    };
    arm("e", v.e);
    arm("p", v.p);
    v.ue = {g("ue_x"), g("ue_y"), g("ue_z")};
    v.up = {g("up_x"), g("up_y"), g("up_z")};
    v.uq = {g("uq_x"), g("uq_y"), g("uq_z")};
}

void readMain(const Row& r, MainEvent& m) {
    m.epsilon = r["in_m_epsilon"];
    m.t = r["in_m_t"];
    m.phi_pq = r["in_m_phi_pq"];
    m.ntup.sigcm1 = r["in_ntup_sigcm1"];
    m.ntup.sigcm2 = r["in_ntup_sigcm2"];
    m.ntup.sigcm = r["in_ntup_sigcm"];
    m.pfermi = FermiMomentum{r["pfer"], r["pferx"], r["pfery"], r["pferz"], r["efer"]};
}

void putPeepi(double sig, const MainEvent& m, Outputs& o) {
    o["peepi"] = sig;
    o["m_thetacm"] = m.thetacm;
    o["m_phicm"] = m.phicm;
    o["m_pcm"] = m.pcm;
    o["m_davejac"] = m.davejac;
    o["m_johnjac"] = m.johnjac;
    o["m_wcm"] = m.wcm;
    o["ntup_sigcm1"] = m.ntup.sigcm1;
    o["ntup_sigcm2"] = m.ntup.sigcm2;
    o["ntup_sigcm"] = m.ntup.sigcm;
}

const std::vector<std::string> kPeepiCols = {"peepi",   "m_thetacm", "m_phicm",     "m_pcm",       "m_davejac",
                                             "m_johnjac", "m_wcm",   "ntup_sigcm1", "ntup_sigcm2", "ntup_sigcm"};

int runPeepi(const std::string& dir, const std::string& stateDir, const std::string& dataDir, double tol) {
    std::map<int, CfgState> st;
    if (!loadStates(stateDir, st)) return 77;
    pion::MaidTables maid;
    std::string err;
    if (!maid.load(dataDir, &err)) {
        std::printf("peepi: %s\n", err.c_str());
        return 1;
    }
    std::map<int, std::size_t> perCfg;
    const int rc = simcl1::run("peepi (physics_pion.f:1-194), fortran_compat", dir, specs(kPeepiCols, tol),
                               [&](const Row& r, Outputs& o) {
                                   const CfgState& c = st.at(iv(r, "cfg"));
                                   ++perCfg[iv(r, "cfg")];
                                   pion::PeepiConfig pc{c.which_pion, c.doing_pizero, c.Mtar_struck, true};
                                   Event v;
                                   MainEvent m;
                                   readEvent(r, "v_", v);
                                   readMain(r, m);
                                   const double sig = pion::peepi(v, m, pc, maid);
                                   putPeepi(sig, m, o);
                               });
    std::printf("  rows per cfg:");
    for (const auto& kv : perCfg) std::printf(" %d:%zu", kv.first, kv.second);
    std::printf("\n");
    return rc;
}

// ------------------------------------------------------------------------------------------------- pizero_decay
int runPizero(const std::string& dir, double tol) {
    std::ostringstream log;
    const int rc = simcl1::run(
        "pizero_decay (pizero_decay.f:1-94), ScriptedRandom", dir,
        specs({"success", "gamma1_1", "gamma1_2", "gamma1_3", "gamma1_4", "gamma2_1", "gamma2_2", "gamma2_3", "gamma2_4",
               "n_unif", "rng_err"},
              tol),
        [&](const Row& r, Outputs& o) {
            ScriptedRandom rng;
            rng.load({}, {r["u1"], r["u2"], r["u3"]});
            Event v;
            v.p.P = r["p_P"];
            v.up = {r["up_x"], r["up_y"], r["up_z"]};
            const pion::PizeroDecayResult d = pion::pizeroDecay(v, r["Mh"], rng, &log);
            o["success"] = d.success ? 1.0 : 0.0;
            for (int k = 0; k < 4; ++k) {
                o["gamma1_" + std::to_string(k + 1)] = d.gamma1[static_cast<std::size_t>(k)];
                o["gamma2_" + std::to_string(k + 1)] = d.gamma2[static_cast<std::size_t>(k)];
            }
            o["n_unif"] = static_cast<double>(rng.nUniform());
            o["rng_err"] = rng.error() ? 1.0 : 0.0;
        });
    std::size_t warn = 0;
    for (std::size_t p = log.str().find("pi0 decay"); p != std::string::npos; p = log.str().find("pi0 decay", p + 1)) ++warn;
    std::printf("  conservation warnings printed: %zu\n", warn);
    return rc;
}

// ------------------------------------------------------------------------------------------------- plug-in
std::unique_ptr<IPhysicsModel> makeModel(const CfgState& c, const std::string& dataDir, bool compat) {
    SimcParameters p;
    p.doing_pion = true;
    p.which_pion = c.which_pion;
    p.doing_pizero = c.doing_pizero;
    RunSetup s = RunSetup::FromParameters(p);
    s.flags.doing_pizero = c.doing_pizero;
    s.targ.Mtar_struck = c.Mtar_struck;
    s.fortran_compat = compat;
    auto m = Registry<IPhysicsModel>::create("pion");
    m->configure(Json{{"data_dir", dataDir}, {"verbose", false}}, p);
    m->init(s);
    return m;
}

int runPlugin(const std::string& dir, const std::string& stateDir, const std::string& dataDir) {
    std::map<int, CfgState> st;
    if (!loadStates(stateDir, st)) return 77;
    simcl1::Table in, ref;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !ref.load(dir + "/output.txt", err)) return 77;
    int bad = 0;
    // a non-pion set-up is refused
    try {
        SimcParameters p;
        RunSetup s = RunSetup::FromParameters(p);
        auto m = Registry<IPhysicsModel>::create("pion");
        m->configure(Json{{"data_dir", dataDir}}, p);
        m->init(s);
        std::printf("  plugin: init without targ%%Mtar_struck did not throw\n");
        ++bad;
    } catch (const std::invalid_argument& e) {
        std::printf("  plugin: non-pion set-up refused: %s\n", e.what());
    }
    std::map<int, std::unique_ptr<IPhysicsModel>> models;
    std::size_t rows = 0, eq = 0;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        Row r(in, i), f(ref, i);
        const int cfg = iv(r, "cfg");
        auto it = models.find(cfg);
        if (it == models.end()) it = models.emplace(cfg, makeModel(st.at(cfg), dataDir, true)).first;
        IPhysicsModel& model = *it->second;
        Event v;
        MainEvent m;
        readEvent(r, "v_", v);
        readMain(r, m);
        const MainEvent before = m;
        const double c1 = model.sigcc(v, m);   // const: main must not change
        const bool unchanged = std::memcmp(&before.thetacm, &m.thetacm, sizeof(double)) == 0 &&
                               std::memcmp(&before.ntup.sigcm, &m.ntup.sigcm, sizeof(double)) == 0;
        const double c2 = model.sigccMain(v, m);
        const double rr = model.sigccRecon(v, m);
        Outputs o;
        putPeepi(c2, m, o);
        bool ok = same(c1, c2) && unchanged && rr == 1.0;
        for (const auto& c : kPeepiCols) ok = ok && same(o[c], f[c]);
        ++rows;
        if (ok) {
            ++eq;
        } else if (bad++ < 5) {
            std::printf("  plugin row %zu (cfg %d): sigcc %.17g sigccMain %.17g ref %.17g unchanged %d recon %g\n", i + 1, cfg,
                        c1, c2, f["peepi"], unchanged ? 1 : 0, rr);
        }
    }
    std::printf("plugin \"pion\" (Registry): %zu rows, %zu bitwise (all outputs) -> %s\n", rows, eq,
                bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------- D-038 F10
int runFortranDiffF10(const std::string& dir, const std::string& stateDir, const std::string& dataDir) {
    std::map<int, CfgState> st;
    if (!loadStates(stateDir, st)) return 77;
    simcl1::Table in, ref;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !ref.load(dir + "/output.txt", err)) return 77;
    pion::MaidTables maid;
    if (!maid.load(dataDir, &err)) {
        std::printf("%s\n", err.c_str());
        return 1;
    }
    std::size_t highW = 0, staleNonzero = 0, lowW = 0, bad = 0;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        Row r(in, i), f(ref, i);
        const CfgState& c = st.at(iv(r, "cfg"));
        Event v;
        readEvent(r, "v_", v);
        MainEvent mc, md;
        readMain(r, mc);
        readMain(r, md);
        const double sc = pion::peepi(v, mc, pion::PeepiConfig{c.which_pion, c.doing_pizero, c.Mtar_struck, true}, maid);
        const double sd = pion::peepi(v, md, pion::PeepiConfig{c.which_pion, c.doing_pizero, c.Mtar_struck, false}, maid);
        // the cross section and every other output are the same in both modes
        bool ok = same(sc, sd) && same(mc.ntup.sigcm, md.ntup.sigcm) && same(mc.ntup.sigcm1, md.ntup.sigcm1) &&
                  same(mc.wcm, md.wcm) && same(sc, f["peepi"]);
        if (mc.wcm < 2000) {
            ++lowW;
            ok = ok && same(mc.ntup.sigcm2, md.ntup.sigcm2) && same(mc.ntup.sigcm2, f["ntup_sigcm2"]);
        } else {
            ++highW;
            // compat: the Fortran keeps the incoming (stale) value; corrected: 0
            ok = ok && same(mc.ntup.sigcm2, r["in_ntup_sigcm2"]) && same(mc.ntup.sigcm2, f["ntup_sigcm2"]) &&
                 md.ntup.sigcm2 == 0.0;
            if (r["in_ntup_sigcm2"] != 0.0) ++staleNonzero;
        }
        if (!ok && bad++ < 5)
            std::printf("  row %zu: W %.6f compat sigcm2 %.17g corrected %.17g ref %.17g\n", i + 1, mc.wcm, mc.ntup.sigcm2,
                        md.ntup.sigcm2, f["ntup_sigcm2"]);
    }
    std::printf("D-038 F10 (physics_pion.f:133-154 stale ntup%%sigcm2): %zu rows W >= 2 GeV (%zu with a non-zero stale "
                "value: compat keeps it, corrected 0), %zu rows W < 2 GeV (same in both modes); peepi/sigcm/sigcm1 "
                "identical in both modes; %zu bad -> %s\n",
                highW, staleNonzero, lowW, bad, (bad || !staleNonzero) ? "FAIL" : "PASS");
    return (bad || !staleNonzero) ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mode> ... (see the file header)\n", argv[0]);
        return 2;
    }
    const std::string what = argv[1];
    try {
        if (what == "files" && argc > 3) return runFiles(argv[2], argv[3]);
        if ((what == "sig_blok" || what == "sig_param04" || what == "sig_param_3000") && argc > 3)
            return runSimpleFit(what, argv[2], std::atof(argv[3]));
        if (what == "sig_param_2021" && argc > 3) return runParam2021(argv[2], std::atof(argv[3]));
        if (what == "exclfit" && argc > 3) return runExclfit(argv[2], std::atof(argv[3]));
        if (what == "sigmaid" && argc > 4) return runSigmaid(argv[2], argv[3], std::atof(argv[4]));
        if (what == "sig_multipole" && argc > 4) return runMultipole(argv[2], argv[3], std::atof(argv[4]));
        if (what == "transform_to_cm" && argc > 3) return runTransform(argv[2], std::atof(argv[3]));
        if (what == "peepi" && argc > 5) return runPeepi(argv[2], argv[3], argv[4], std::atof(argv[5]));
        if (what == "pizero_decay" && argc > 3) return runPizero(argv[2], std::atof(argv[3]));
        if (what == "plugin" && argc > 4) return runPlugin(argv[2], argv[3], argv[4]);
        if (what == "fortran_diff_D038_F10" && argc > 4) return runFortranDiffF10(argv[2], argv[3], argv[4]);
    } catch (const std::exception& e) {
        std::printf("ERROR: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "bad arguments for mode %s\n", what.c_str());
    return 2;
}
