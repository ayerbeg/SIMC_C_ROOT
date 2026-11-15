// src/optics/ForwardChains.cpp — see include/simc/optics/ForwardChains.h.
// Aperture-free forward transport target -> focal plane for the HMS and SHMS (the drift/map/rotation sequence of
// mc_hms.f and mc_shms.f with every aperture test removed), and pure COSY chains. Used by the L1 tests and the optics
// diagnostics to check the maps independently of the acceptance; a callback sees the track at every plane.
#include "simc/optics/ForwardChains.h"

namespace simc {
namespace optics {

namespace {

struct Stepper {
    Track& t;
    double& pathlen;
    const ChainCallback& cb;
    bool stopped{false};
    bool rec(const char* name, double xt, double yt) {
        if (stopped) return false;
        ChainStep s{name, t, pathlen, xt, yt};
        if (cb && !cb(s)) stopped = true;
        return !stopped;
    }
    bool rec(const char* name) { return rec(name, t.xs, t.ys); }
};

}  // namespace

// Fortran: hms/mc_hms.f:201-413 (constants mc_hms.f:84-94), apertures removed.
bool hmsChain(CosyEngine& eng, Track& t, double& pathlen, const ChainCallback& cb, DecayState* decay) {
    constexpr int spectr = 1;
    constexpr double z_off = +40.17;
    constexpr double z_entr = 126.2e0 + z_off;   // nominally 1.262 m
    constexpr double z_exit = z_entr + 6.3e0;    // 6.3 cm thick
    constexpr double z_dip1 = 64.77e0;
    constexpr double z_dip2 = z_dip1 + 297.18e0;
    constexpr double z_dip3 = z_dip2 + 115.57e0;
    Stepper s{t, pathlen, cb};
    double xt, yt;

    project(t, z_entr, pathlen, decay);                              // front of fixed slit
    if (!s.rec("coll_ent")) return false;
    project(t, z_exit - z_entr, pathlen, decay);                     // back of fixed slit
    if (!s.rec("coll_ext")) return false;
    project(t, eng.driftdist(spectr, 1) - z_exit, pathlen, decay);   // Q1 IN mag bound
    if (!s.rec("q1_in")) return false;
    eng.transp(spectr, 2, t, 125.233e0, pathlen, decay);
    if (!s.rec("q1_mid")) return false;
    eng.transp(spectr, 3, t, 62.617e0, pathlen, decay);
    if (!s.rec("q1_out")) return false;
    project(t, eng.driftdist(spectr, 4), pathlen, decay);
    if (!s.rec("q2_in")) return false;
    eng.transp(spectr, 5, t, 143.90e0, pathlen, decay);
    if (!s.rec("q2_mid")) return false;
    eng.transp(spectr, 6, t, 71.95e0, pathlen, decay);
    if (!s.rec("q2_out")) return false;
    project(t, eng.driftdist(spectr, 7), pathlen, decay);
    if (!s.rec("q3_in")) return false;
    eng.transp(spectr, 8, t, 143.8e0, pathlen, decay);
    if (!s.rec("q3_mid")) return false;
    eng.transp(spectr, 9, t, 71.9e0, pathlen, decay);
    if (!s.rec("q3_out")) return false;
    project(t, eng.driftdist(spectr, 10), pathlen, decay);            // D1 IN
    xt = t.xs;
    yt = t.ys;
    rotateHaxis(t, -6.0e0, xt, yt, eng.fortranCompat());
    if (!s.rec("d1_in", xt, yt)) return false;
    eng.transp(spectr, 11, t, 526.053e0, pathlen, decay);             // D1 OUT
    xt = t.xs;
    yt = t.ys;
    rotateHaxis(t, 6.0e0, xt, yt, eng.fortranCompat());
    if (!s.rec("d1_out", xt, yt)) return false;
    project(t, z_dip1, pathlen, decay);
    if (!s.rec("pipe1")) return false;
    project(t, z_dip2 - z_dip1, pathlen, decay);
    if (!s.rec("pipe2")) return false;
    project(t, z_dip3 - z_dip2, pathlen, decay);
    if (!s.rec("pipe3")) return false;
    project(t, eng.driftdist(spectr, 12) - z_dip3, pathlen, decay);   // to the focal plane (mc_hms.f:413)
    return s.rec("fp");
}

// Fortran: shms/mc_shms.f:322-1047 (use_coll, no skip_hb), apertures removed; zd_* PARAMETERs of mc_shms.f.
bool shmsChain(CosyEngine& eng, Track& t, double& pathlen, const ChainCallback& cb, DecayState* decay) {
    constexpr int spectr = 5;
    constexpr double zd_hbin = 118.39, zd_hbmen = 17.61, zd_hbmex = 80.0, zd_hbout = 17.61;
    constexpr double z_entr = 25.189, z_thick = 6.35;
    constexpr double zd_q1in = 58.39, zd_q1men = 28.35, zd_q1mid = 93.65, zd_q1mex = 93.65, zd_q1out = 28.35;
    constexpr double zd_q2in = 25.55, zd_q2men = 39.1, zd_q2mid = 79.35, zd_q2mex = 79.35, zd_q2out = 39.1;
    constexpr double zd_q3in = 28.10, zd_q3men = 39.1, zd_q3mid = 79.35, zd_q3mex = 79.35, zd_q3out = 39.1;
    constexpr double zd_q3d1trans = 18.00, zd_d1flare = 30.10;
    constexpr double zd_d1men = 39.47;   // FORTRAN-DIFF-CANDIDATE(D-015 F5): COSY class 22 drift is 39.4715
    constexpr double zd_d1mid = 36.406263;
    constexpr double zd_d1mex = 36.406263;
    constexpr double zd_d1out = 60.68;   // COSY class 31: 60.68058
    constexpr double zd_fp = 307.95;     // COSY class 32: 307.9514
    Stepper s{t, pathlen, cb};
    double xt, yt;
    auto vrot = [&](double ang, double yoff) {
        xt = t.xs;
        yt = t.ys;
        rotateVaxis(t, ang, xt, yt, eng.fortranCompat());
        yt = yt + yoff;
    };
    auto hrot = [&](double ang, double xoff) {
        xt = t.xs;
        yt = t.ys;
        rotateHaxis(t, ang, xt, yt, eng.fortranCompat());
        xt = xt + xoff;
    };
    // HB
    project(t, zd_hbin, pathlen, decay);
    vrot(1.5, 1.51);
    if (!s.rec("hb_in", xt, yt)) return false;
    project(t, zd_hbmen, pathlen, decay);
    vrot(1.5, 0.98);
    if (!s.rec("hb_men", xt, yt)) return false;
    eng.transp(spectr, 3, t, zd_hbmex, pathlen, decay);
    vrot(-1.5, 0.98);
    if (!s.rec("hb_mex", xt, yt)) return false;
    project(t, zd_hbout, pathlen, decay);
    vrot(-1.5, 1.51);
    if (!s.rec("hb_out", xt, yt)) return false;
    // collimator (use_coll)
    project(t, z_entr, pathlen, decay);
    if (!s.rec("coll_ent")) return false;
    project(t, z_thick, pathlen, decay);
    if (!s.rec("coll_ext")) return false;
    // Q1
    project(t, zd_q1in - z_entr - z_thick, pathlen, decay);
    if (!s.rec("q1_in")) return false;
    project(t, zd_q1men, pathlen, decay);
    if (!s.rec("q1_men")) return false;
    eng.transp(spectr, 7, t, zd_q1mid, pathlen, decay);
    if (!s.rec("q1_mid")) return false;
    eng.transp(spectr, 8, t, zd_q1mex, pathlen, decay);
    if (!s.rec("q1_mex")) return false;
    project(t, zd_q1out, pathlen, decay);
    if (!s.rec("q1_out")) return false;
    // Q2
    project(t, zd_q2in, pathlen, decay);
    if (!s.rec("q2_in")) return false;
    project(t, zd_q2men, pathlen, decay);
    if (!s.rec("q2_men")) return false;
    eng.transp(spectr, 12, t, zd_q2mid, pathlen, decay);
    if (!s.rec("q2_mid")) return false;
    eng.transp(spectr, 13, t, zd_q2mex, pathlen, decay);
    if (!s.rec("q2_mex")) return false;
    project(t, zd_q2out, pathlen, decay);
    if (!s.rec("q2_out")) return false;
    // Q3
    project(t, zd_q3in, pathlen, decay);
    if (!s.rec("q3_in")) return false;
    project(t, zd_q3men, pathlen, decay);
    if (!s.rec("q3_men")) return false;
    eng.transp(spectr, 17, t, zd_q3mid, pathlen, decay);
    if (!s.rec("q3_mid")) return false;
    eng.transp(spectr, 18, t, zd_q3mex, pathlen, decay);
    if (!s.rec("q3_mex")) return false;
    project(t, zd_q3out, pathlen, decay);
    if (!s.rec("q3_out")) return false;
    // D1
    project(t, zd_q3d1trans, pathlen, decay);
    if (!s.rec("d1_in")) return false;
    project(t, zd_d1flare, pathlen, decay);
    hrot(9.200, -3.5);
    if (!s.rec("d1_flr", xt, yt)) return false;
    project(t, zd_d1men, pathlen, decay);
    hrot(9.200, 2.82);
    if (!s.rec("d1_men", xt, yt)) return false;
    struct Mid { int cls; double ang, xoff; const char* name; };
    static const Mid mids[] = {{23, 6.9, 8.05, "d1_m1"},   {24, 4.6, 11.75, "d1_m2"},  {25, 2.3, 13.96, "d1_m3"},
                               {26, 0.0, 14.70, "d1_m4"},  {27, -2.3, 13.96, "d1_m5"}, {28, -4.6, 11.75, "d1_m6"},
                               {29, -6.9, 8.05, "d1_m7"}};
    for (const Mid& m : mids) {
        eng.transp(spectr, m.cls, t, zd_d1mid, pathlen, decay);
        hrot(m.ang, m.xoff);
        if (!s.rec(m.name, xt, yt)) return false;
    }
    eng.transp(spectr, 30, t, zd_d1mex, pathlen, decay);
    hrot(-9.2, 2.82);
    if (!s.rec("d1_mex", xt, yt)) return false;
    project(t, zd_d1out, pathlen, decay);
    hrot(-9.20, -6.88);
    if (!s.rec("d1_out", xt, yt)) return false;
    project(t, zd_fp, pathlen, decay);
    return s.rec("fp");
}

// zd passed to transp by mc_hms.f / mc_shms.f for magnet classes (reference driver drv_chain.f:zdm).
double mcMagnetZd(int spectr, int k) {
    if (spectr == 1) {
        switch (k) {
            case 2: return 125.233e0;
            case 3: return 62.617e0;
            case 5: return 143.90e0;
            case 6: return 71.95e0;
            case 8: return 143.8e0;
            case 9: return 71.9e0;
            case 11: return 526.053e0;
            default: return 0.;
        }
    }
    if (spectr == 5) {
        if (k == 3) return 80.0;
        if (k == 7 || k == 8) return 93.65;
        if (k == 12 || k == 13 || k == 17 || k == 18) return 79.35;
        if (k >= 23 && k <= 30) return 36.406263;
    }
    return 0.;
}

void cosyChain(CosyEngine& eng, int spectr, Track& t, double& pathlen) {
    const int n = eng.nClasses(spectr);
    for (int k = 1; k <= n; ++k) {
        const double zd = eng.adrift(spectr, k) ? eng.driftdist(spectr, k) : mcMagnetZd(spectr, k);
        eng.transp(spectr, k, t, zd, pathlen, nullptr);
    }
}

std::vector<ChainStep> collectHmsChain(CosyEngine& eng, Track t) {
    std::vector<ChainStep> v;
    double path = 0;
    hmsChain(eng, t, path, [&](const ChainStep& s) {
        v.push_back(s);
        return true;
    });
    return v;
}

std::vector<ChainStep> collectShmsChain(CosyEngine& eng, Track t) {
    std::vector<ChainStep> v;
    double path = 0;
    shmsChain(eng, t, path, [&](const ChainStep& s) {
        v.push_back(s);
        return true;
    });
    return v;
}

}  // namespace optics
}  // namespace simc
