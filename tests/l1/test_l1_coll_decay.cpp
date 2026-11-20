// tests/l1/test_l1_coll_decay.cpp
//
// L1 tests of the HMS/SHMS collimator (pions/muons) and decay-in-flight paths (task P6.4) against the Fortran harness
// tables of P6.F2 (reference/functions/README.md, section P6.F2; reference/P6.F2_report.md). All comparisons are
// BITWISE (tolerance 0), in fortran_compat mode (D-031):
//   mc <table_dir> hms|shms        : mc_hms_coll, mc_hms_decay, mc_shms_coll, mc_shms_decay: every output column of
//                                    drv_mccd.f (arguments, /track/, decdist, Mh2_final, save_coll, RNG consumption,
//                                    collimator columns n_absorb_calls..eloss_sum, decay columns n_decay_checks..
//                                    dec_step_len, every aperture-counter delta d_*), plus every project/transp step of
//                                    steps.txt.xz (entry and exit state) and every pion_coll_absorb/enerloss_new call
//                                    of coll_calls.txt.
//   tp <table_dir>                 : transp_decay_hms, transp_decay_shms, project_decay: one project/transp call.
//   pca <table_dir>                : pion_coll_absorb.
//   replay <cap_calls.txt> [corrected] : every mc_hms/mc_shms call recorded in a P6.F2 capture run
//                                    (fortran_harness/spectrometers_coll_decay/work/p6f2_pi_*_cap/cap_calls.txt), in call
//                                    order through ONE McHms and ONE McShms (one Fortran process), with the call's own
//                                    RANLUX values: outputs, counters, RNG consumption. "corrected" runs the default mode
//                                    and reports the effect of the D-037 corrections instead (no pass/fail on the diff).
// Exit: 0 pass, 1 fail, 77 table missing.
#include <cmath>
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
#include "simc/core/random/ScriptedRandom.h"
#include "simc/optics/CosyEngine.h"
#include "simc/optics/TransportTrace.h"
#include "simc/spectrometers/hms/kernel/McHms.h"
#include "simc/spectrometers/shms/kernel/McShms.h"

using namespace simc;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

constexpr double kUnset = -999999.;

bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0 || (std::isnan(a) && std::isnan(b)); }

// ------------------------------------------------------------------------------------------------ one arm
// The arguments common to mc_hms and mc_shms, and the kernel call with all diagnostics.
struct CallIn {
    double p_spec{0}, th_spec{0}, dpp{0}, x{0}, y{0}, z{0}, dxdz{0}, dydz{0}, m2{0}, fry{0};
    bool ms{false}, wcs{false}, decay{false}, coll{false};
    double ctau{0}, decdist{0}, Mh2_final{0}, resmult{0}, pathlen{0};
    int spectr{1};
};
struct CallOut {
    bool ok{false};
    int stop_code{0}, stop_id{-1};
    bool reached_hut{false};
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0}, dpp{0}, dxdz{0}, dydz{0}, y{0}, pathlen{0}, m2{0}, resmult{0};
    optics::Track t;
    double decdist{0}, Mh2_final{0}, save_coll{0};
    std::vector<long> dcnt;   ///< counter deltas, index 1..N (HMS 22, SHMS 52 with slot 51 = stop_id after)
};

class Arm {
public:
    explicit Arm(bool isHms, bool compat) : hms_(isHms) {
        auto eng = std::make_shared<optics::CosyEngine>();
        eng->setLog(&log_);
        if (hms_) {
            h_ = std::make_unique<hms::McHms>();
            h_->setEngine(eng);
            h_->setFortranCompat(compat);
            h_->setLog(&log_);
            h_->setTransportTrace(&tt_);
            h_->setCollTrace(&ct_);
        } else {
            s_ = std::make_unique<shms::McShms>();
            s_->setEngine(eng);
            s_->setFortranCompat(compat);
            s_->setLog(&log_);
            s_->setTransportTrace(&tt_);
            s_->setCollTrace(&ct_);
        }
        tt_.steps = &steps_;
        ct_.calls = &calls_;
    }
    bool isHms() const { return hms_; }
    void setSave(double v) {
        if (hms_) h_->setSaveCollSuccess(v);
        else s_->setSaveCollSuccess(v);
    }
    CallOut call(const CallIn& in, IRandom& rng) {
        CallOut o;
        if (hms_) {
            hms::McHmsArgs a;
            a.p_spec = in.p_spec;
            a.th_spec = in.th_spec;
            a.dpp = in.dpp;
            a.x = in.x;
            a.y = in.y;
            a.z = in.z;
            a.dxdz = in.dxdz;
            a.dydz = in.dydz;
            a.m2 = in.m2;
            a.fry = in.fry;
            a.ms_flag = in.ms;
            a.wcs_flag = in.wcs;
            a.decay_flag = in.decay;
            a.using_HMScoll = in.coll;
            a.ctau = in.ctau;
            a.decdist = in.decdist;
            a.Mh2_final = in.Mh2_final;
            a.resmult = in.resmult;
            a.pathlen = in.pathlen;
            a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
            const hms::HmsCounters before = h_->counters();
            const hms::McHmsDetail det = h_->mcHms(a, rng);
            const hms::HmsCounters& after = h_->counters();
            o.dcnt.assign(hms::kNumCounters + 1, 0);
            for (int k = 1; k <= hms::kNumCounters; ++k) o.dcnt[k] = after[k] - before[k];
            o.ok = a.ok_spec;
            o.stop_code = det.stop_code;
            o.stop_id = -1;
            o.reached_hut = det.reached_hut;
            o.x_fp = a.x_fp;
            o.dx_fp = a.dx_fp;
            o.y_fp = a.y_fp;
            o.dy_fp = a.dy_fp;
            o.dpp = a.dpp;
            o.dxdz = a.dxdz;
            o.dydz = a.dydz;
            o.y = a.y;
            o.pathlen = a.pathlen;
            o.m2 = a.m2;
            o.resmult = a.resmult;
            o.t = h_->track();
            o.decdist = a.decdist;
            o.Mh2_final = a.Mh2_final;
            o.save_coll = h_->saveCollSuccess();
        } else {
            shms::McShmsArgs a;
            a.spectr = in.spectr;
            a.p_spec = in.p_spec;
            a.th_spec = in.th_spec;
            a.dpp = in.dpp;
            a.x = in.x;
            a.y = in.y;
            a.z = in.z;
            a.dxdz = in.dxdz;
            a.dydz = in.dydz;
            a.m2 = in.m2;
            a.fry = in.fry;
            a.ms_flag = in.ms;
            a.wcs_flag = in.wcs;
            a.decay_flag = in.decay;
            a.using_SHMScoll = in.coll;
            a.ctau = in.ctau;
            a.decdist = in.decdist;
            a.Mh2_final = in.Mh2_final;
            a.resmult = in.resmult;
            a.pathlen = in.pathlen;
            a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
            const shms::ShmsCounters before = s_->counters();
            const shms::McShmsDetail det = s_->mcShms(a, rng);
            const shms::ShmsCounters& after = s_->counters();
            o.dcnt.assign(shms::kNumSlots + 1, 0);
            for (int k = 1; k <= shms::kNumSlots; ++k)
                o.dcnt[k] = (k == shms::kStopIdSlot) ? after[k] : after[k] - before[k];
            o.ok = a.ok_spec;
            o.stop_code = det.stop_code;
            o.stop_id = det.stop_id;
            o.reached_hut = det.reached_hut;
            o.x_fp = a.x_fp;
            o.dx_fp = a.dx_fp;
            o.y_fp = a.y_fp;
            o.dy_fp = a.dy_fp;
            o.dpp = a.dpp;
            o.dxdz = a.dxdz;
            o.dydz = a.dydz;
            o.y = a.y;
            o.pathlen = a.pathlen;
            o.m2 = a.m2;
            o.resmult = a.resmult;
            o.t = s_->track();
            o.decdist = a.decdist;
            o.Mh2_final = a.Mh2_final;
            o.save_coll = s_->saveCollSuccess();
        }
        return o;
    }
    const optics::TransportTrace& tt() const { return tt_; }
    const optics::CollTrace& ct() const { return ct_; }
    const std::vector<optics::StepRecord>& steps() const { return steps_; }
    const std::vector<optics::CollCall>& calls() const { return calls_; }
    std::string takeLog() {
        std::string s = log_.str();
        log_.str("");
        return s;
    }

private:
    bool hms_;
    std::ostringstream log_;
    std::unique_ptr<hms::McHms> h_;
    std::unique_ptr<shms::McShms> s_;
    optics::TransportTrace tt_;
    optics::CollTrace ct_;
    std::vector<optics::StepRecord> steps_;
    std::vector<optics::CollCall> calls_;
};

const char* const kHmsNames[23] = {"",      "trials", "slit_hor", "slit_vert", "slit_oct", "Q1_in",  "Q1_mid", "Q1_out",
                                   "Q2_in", "Q2_mid", "Q2_out",   "Q3_in",     "Q3_mid",   "Q3_out", "D1_in",  "D1_out",
                                   "hut",   "dc1",    "dc2",      "scin",      "cal",      "successes", "coll"};
const char* const kShmsNames[53] = {
    "",        "trials",  "successes", "targ_hor", "targ_vert", "targ_oct", "slit_hor", "slit_vert", "slit_oct",
    "HB_in",   "HB_men",  "HB_mex",    "HB_out",   "Q1_in",     "Q1_men",   "Q1_mid",   "Q1_mex",    "Q1_out",
    "Q2_in",   "Q2_men",  "Q2_mid",    "Q2_mex",   "Q2_out",    "Q3_in",    "Q3_men",   "Q3_mid",    "Q3_mex",
    "Q3_out",  "D1_in",   "D1_flr",    "D1_men",   "D1_mid1",   "D1_mid2",  "D1_mid3",  "D1_mid4",   "D1_mid5",
    "D1_mid6", "D1_mid7", "D1_mex",    "D1_out",   "BP_in",     "BP_out",   "hut",      "dc1",       "dc2",
    "cer1",    "cer2",    "s1",        "s2",       "s3",        "cal",      "stop_id",  "coll"};

void loadQueues(const Row& r, ScriptedRandom& rng) {
    std::vector<double> g, u;
    for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
    for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
    rng.load(g, u);
}

void fillDecay(Outputs& o, const optics::TransportTrace& tt, bool withIndices) {
    const optics::DecayRecord& d = tt.dec;
    const bool z = !tt.decayed;
    auto v = [&](double x) { return z ? 0.0 : x; };
    o["decayed"] = tt.decayed ? 1. : 0.;
    o["dec_kind"] = v(d.kind);
    if (withIndices) {
        o["dec_check"] = v(static_cast<double>(tt.dec_check));
        o["dec_step"] = v(static_cast<double>(tt.dec_step));
        o["dec_class"] = v(d.cls);
    }
    o["dec_x"] = v(d.x);
    o["dec_y"] = v(d.y);
    o["dec_dxdz_before"] = v(d.dxdz_before);
    o["dec_dydz_before"] = v(d.dydz_before);
    o["dec_dpp_before"] = v(d.dpp_before);
    o["dec_p_before"] = v(d.p_before);
    o["dec_pathlen_at"] = v(d.pathlen_at);
    o["dec_decdist_at"] = v(d.decdist_at);
    o["dec_m2_before"] = v(d.m2_before);
    o["dec_gamma"] = v(d.gamma);
    o["dec_bx"] = v(d.bx);
    o["dec_by"] = v(d.by);
    o["dec_bz"] = v(d.bz);
    o["dec_er"] = v(d.er);
    o["dec_pxr"] = v(d.pxr);
    o["dec_pyr"] = v(d.pyr);
    o["dec_pzr"] = v(d.pzr);
    o["dec_ef"] = v(d.ef);
    o["dec_pxf"] = v(d.pxf);
    o["dec_pyf"] = v(d.pyf);
    o["dec_pzf"] = v(d.pzf);
    o["dec_pf"] = v(d.pf);
    o["dec_dxdz_after"] = v(d.dxdz_after);
    o["dec_dydz_after"] = v(d.dydz_after);
    if (withIndices) {
        o["dec_dpp_after"] = v(d.dpp_after);
        o["dec_m2_after"] = v(d.m2_after);
        o["dec_p_after"] = v(d.p_after);
        o["dec_step_len"] = v(d.step_len);
    }
}
const char* const kDecCols[] = {"decayed", "dec_kind", "dec_x", "dec_y", "dec_dxdz_before", "dec_dydz_before",
                                "dec_dpp_before", "dec_p_before", "dec_pathlen_at", "dec_decdist_at", "dec_m2_before",
                                "dec_gamma", "dec_bx", "dec_by", "dec_bz", "dec_er", "dec_pxr", "dec_pyr", "dec_pzr",
                                "dec_ef", "dec_pxf", "dec_pyf", "dec_pzf", "dec_pf", "dec_dxdz_after",
                                "dec_dydz_after"};
const char* const kDecColsMc[] = {"dec_check", "dec_step", "dec_class", "dec_dpp_after", "dec_m2_after", "dec_p_after",
                                  "dec_step_len"};

// ------------------------------------------------------------------------------------------------ steps / coll calls
struct StepRef {
    std::vector<double> v;   // the 28 columns after "row"
};

// steps.txt.xz: row step kind class zd decay_check dflag_in dflag_out n_unif nu_before, then entry and exit
// xs ys dxdzs dydzs dpps m2 ph pathlen decdist.
bool loadSteps(const std::string& path, std::map<long, std::vector<std::vector<double>>>& out) {
    std::ifstream probe(path);
    if (!probe) return false;
    const std::string cmd = "xz -dc '" + path + "'";
    FILE* f = popen(cmd.c_str(), "r");
    if (!f) return false;
    char buf[4096];
    bool header = true;
    while (std::fgets(buf, sizeof buf, f)) {
        if (header) {
            header = false;
            continue;
        }
        std::istringstream ss(buf);
        std::vector<double> v;
        double x;
        while (ss >> x) v.push_back(x);
        if (v.size() < 2) continue;
        out[std::lround(v[0])].push_back(v);
    }
    return pclose(f) == 0;
}

// coll_calls.txt: row type idx a1..a6 ('#' comments).
bool loadCollCalls(const std::string& path, std::map<long, std::vector<std::vector<double>>>& out) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    bool header = true;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (header) {
            header = false;
            continue;
        }
        std::istringstream ss(line);
        std::vector<double> v;
        double x;
        while (ss >> x) v.push_back(x);
        if (v.size() < 9) continue;
        out[std::lround(v[0])].push_back(v);
    }
    return true;
}

struct TraceCheck {
    long rows{0}, stepsRef{0}, stepsBad{0}, rowsStepCountBad{0}, callsRef{0}, callsBad{0}, rowsCallCountBad{0};
    long firstBadRow{0};
};

void checkSteps(long row, const std::vector<optics::StepRecord>& cs, const std::vector<std::vector<double>>* ref,
                TraceCheck& tc) {
    ++tc.rows;
    const std::size_t n = ref ? ref->size() : 0;
    tc.stepsRef += static_cast<long>(n);
    if (cs.size() != n) {
        ++tc.rowsStepCountBad;
        if (!tc.firstBadRow) tc.firstBadRow = row;
    }
    for (std::size_t i = 0; i < std::min(n, cs.size()); ++i) {
        const std::vector<double>& r = (*ref)[i];
        const optics::StepRecord& s = cs[i];
        const double c[] = {static_cast<double>(s.kind), static_cast<double>(s.cls), s.zd, s.decay_check ? 1. : 0.,
                            s.dflag_in ? 1. : 0., s.dflag_out ? 1. : 0.};
        bool ok = true;
        for (int k = 0; k < 6; ++k) ok = ok && sameBits(c[k], r[2 + k]);
        const optics::StepState* st[2] = {&s.in, &s.out};
        for (int h = 0; h < 2; ++h) {
            const optics::StepState& q = *st[h];
            const double qv[9] = {q.xs, q.ys, q.dxdzs, q.dydzs, q.dpps, q.m2, q.ph, q.pathlen, q.decdist};
            for (int k = 0; k < 9; ++k) ok = ok && sameBits(qv[k], r[10 + 9 * h + k]);
        }
        if (!ok) {
            ++tc.stepsBad;
            if (!tc.firstBadRow) tc.firstBadRow = row;
        }
    }
}

// coll_calls.txt lists, per row, the pion_coll_absorb calls (type 1) and the enerloss_new calls (type 2), each type in
// call order (idx); the two sequences are compared separately.
void checkCalls(long row, const std::vector<optics::CollCall>& cs, const std::vector<std::vector<double>>* ref,
                TraceCheck& tc) {
    for (int type = 1; type <= 2; ++type) {
        std::vector<const optics::CollCall*> c;
        std::vector<const std::vector<double>*> r;
        for (const auto& x : cs)
            if (x.type == type) c.push_back(&x);
        if (ref)
            for (const auto& x : *ref)
                if (std::lround(x[1]) == type) r.push_back(&x);
        tc.callsRef += static_cast<long>(r.size());
        if (c.size() != r.size()) {
            ++tc.rowsCallCountBad;
            if (!tc.firstBadRow) tc.firstBadRow = row;
        }
        for (std::size_t i = 0; i < std::min(c.size(), r.size()); ++i) {
            const int na = type == 1 ? 5 : 4;   // type 2: the gauss1 value (a5) is not recorded by the C++
            bool ok = true;
            for (int k = 0; k < na; ++k) ok = ok && sameBits(c[i]->a[k], (*r[i])[3 + k]);
            if (!ok) {
                ++tc.callsBad;
                if (!tc.firstBadRow) tc.firstBadRow = row;
            }
        }
    }
}

// ------------------------------------------------------------------------------------------------ mc tables
int runMc(const std::string& dir, bool isHms) {
    Arm arm(isHms, /*compat=*/true);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::map<long, std::vector<std::vector<double>>> stepsRef, callsRef;
    const bool haveSteps = loadSteps(dir + "/steps.txt.xz", stepsRef);
    const bool haveCalls = loadCollCalls(dir + "/coll_calls.txt", callsRef);
    TraceCheck tc;
    long msgRows = 0;

    std::vector<ColumnSpec> c;
    auto add = [&](const char* n) { c.push_back({n, 0.0, 0.0}); };
    for (const char* n : {"ok", "stop_code", "stop_id", "reached_hut", "x_fp", "dx_fp", "y_fp", "dy_fp", "dpp_out",
                          "dxdz_out", "dydz_out", "y_out", "pathlen", "m2_out", "resmult", "xs_end", "ys_end",
                          "dxdzs_end", "dydzs_end", "dpps_end", "decdist", "Mh2_final", "save_Hcollsuccess",
                          "save_Scollsuccess", "n_gauss", "n_unif", "rng_err", "nsig_last", "n_absorb_calls",
                          "absorbed", "ab_p_last", "ab_thick_last", "ab_trans_last", "ab_u_last", "n_eloss_calls",
                          "eloss_sum", "n_decay_checks", "n_steps"})
        add(n);
    for (const char* n : kDecCols) add(n);
    for (const char* n : kDecColsMc) add(n);
    const int ncnt = isHms ? 22 : 52;
    std::vector<std::string> dnames;
    for (int k = 1; k <= ncnt; ++k) {
        if (!isHms && k == shms::kStopIdSlot) continue;
        dnames.push_back(std::string("d_") + (isHms ? kHmsNames[k] : kShmsNames[k]));
        add(dnames.back().c_str());
    }
    add("row");
    // (the ColumnSpec names are copied: dnames may reallocate, so rebuild the d_ entries from their own strings)
    {
        std::size_t j = 0;
        for (auto& s : c)
            if (s.name.rfind("d_", 0) == 0) s.name = dnames[j++];
    }

    int rc = simcl1::run(std::string("coll/decay mc ") + (isHms ? "hms" : "shms") + " (" + dir + ")", dir, c,
                         [&](const Row& r, Outputs& o) {
        loadQueues(r, rng);
        CallIn in;
        in.spectr = static_cast<int>(std::lround(r["spectr"]));
        in.p_spec = r["p_spec"];
        in.th_spec = r["th_spec"];
        in.dpp = r["dpp"];
        in.x = r["x"];
        in.y = r["y"];
        in.z = r["z"];
        in.dxdz = r["dxdz"];
        in.dydz = r["dydz"];
        in.m2 = r["m2"];
        in.fry = r["fry"];
        in.ms = std::lround(r["ms_flag"]) != 0;
        in.wcs = std::lround(r["wcs_flag"]) != 0;
        in.decay = std::lround(r["decay_flag"]) != 0;
        in.coll = std::lround(r["using_coll"]) != 0;
        in.ctau = r["ctau"];
        in.decdist = r["decdist_in"];
        in.Mh2_final = r["Mh2_final_in"];
        in.resmult = r["resmult_in"];
        in.pathlen = 0.;
        arm.setSave(isHms ? r["save_Hcoll_in"] : r["save_Scoll_in"]);
        const CallOut x = arm.call(in, rng);
        if (!arm.takeLog().empty()) ++msgRows;
        o["ok"] = x.ok ? 1. : 0.;
        o["stop_code"] = x.stop_code;
        o["stop_id"] = x.stop_id;
        o["reached_hut"] = x.reached_hut ? 1. : 0.;
        o["x_fp"] = x.x_fp;
        o["dx_fp"] = x.dx_fp;
        o["y_fp"] = x.y_fp;
        o["dy_fp"] = x.dy_fp;
        o["dpp_out"] = x.dpp;
        o["dxdz_out"] = x.dxdz;
        o["dydz_out"] = x.dydz;
        o["y_out"] = x.y;
        o["pathlen"] = x.pathlen;
        o["m2_out"] = x.m2;
        o["resmult"] = x.resmult;
        o["xs_end"] = x.t.xs;
        o["ys_end"] = x.t.ys;
        o["dxdzs_end"] = x.t.dxdzs;
        o["dydzs_end"] = x.t.dydzs;
        o["dpps_end"] = x.t.dpps;
        o["decdist"] = x.decdist;
        o["Mh2_final"] = x.Mh2_final;
        o["save_Hcollsuccess"] = isHms ? x.save_coll : r["save_Hcoll_in"];
        o["save_Scollsuccess"] = isHms ? r["save_Scoll_in"] : x.save_coll;
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        o["nsig_last"] = rng.lastNsig();
        const optics::CollTrace& ct = arm.ct();
        o["n_absorb_calls"] = static_cast<double>(ct.n_absorb_calls);
        o["absorbed"] = ct.absorbed ? 1. : 0.;
        o["ab_p_last"] = ct.ab_p_last;
        o["ab_thick_last"] = ct.ab_thick_last;
        o["ab_trans_last"] = ct.ab_trans_last;
        o["ab_u_last"] = ct.ab_u_last;
        o["n_eloss_calls"] = static_cast<double>(ct.n_eloss_calls);
        o["eloss_sum"] = ct.eloss_sum;
        const optics::TransportTrace& tt = arm.tt();
        o["n_decay_checks"] = static_cast<double>(tt.n_decay_checks);
        o["n_steps"] = static_cast<double>(tt.n_steps);
        fillDecay(o, tt, true);
        std::size_t j = 0;
        for (int k = 1; k <= ncnt; ++k) {
            if (!isHms && k == shms::kStopIdSlot) continue;
            o[dnames[j++]] = static_cast<double>(x.dcnt[k]);
        }
        o["row"] = static_cast<double>(r.number());
        const long row = static_cast<long>(r.number());
        if (haveSteps) {
            auto it = stepsRef.find(row);
            checkSteps(row, arm.steps(), it == stepsRef.end() ? nullptr : &it->second, tc);
        }
        if (haveCalls) {
            auto it = callsRef.find(row);
            checkCalls(row, arm.calls(), it == callsRef.end() ? nullptr : &it->second, tc);
        }
    });
    if (rc == 77) return rc;
    bool ok = rc == 0;
    std::printf("  rows with Fortran write(6,*) messages (musc/enerloss warnings, expected for ranged-out rows): %ld\n",
                msgRows);
    if (haveSteps) {
        std::printf("  steps.txt.xz: %ld rows, %ld reference steps, %ld rows with a different step count, %ld steps "
                    "differ (bitwise, entry+exit state)\n",
                    tc.rows, tc.stepsRef, tc.rowsStepCountBad, tc.stepsBad);
        if (tc.rowsStepCountBad || tc.stepsBad) ok = false;
    } else {
        std::printf("  steps.txt.xz: not found (skipped)\n");
    }
    if (haveCalls) {
        std::printf("  coll_calls.txt: %ld reference calls, %ld rows with a different call count, %ld calls differ\n",
                    tc.callsRef, tc.rowsCallCountBad, tc.callsBad);
        if (tc.rowsCallCountBad || tc.callsBad) ok = false;
    }
    if (tc.firstBadRow) std::printf("  first row with a trace difference: %ld\n", tc.firstBadRow);
    std::printf("  => %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// ------------------------------------------------------------------------------------------------ tp tables
int runTp(const std::string& dir) {
    auto eng = std::make_shared<optics::CosyEngine>();
    std::ostringstream log;
    eng->setLog(&log);
    eng->setFortranCompat(true);
    eng->setCompat294f662(true);   // U.R: 294f662 reference (D-006 kaon decay mass Mk)
    const std::string mdir = hms::McHms::defaultMatrixDir();
    eng->initForward(1, mdir);
    eng->initForward(5, mdir);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    optics::TransportTrace tt;

    std::vector<ColumnSpec> c;
    for (const char* n : {"xs_out", "dxdzs_out", "ys_out", "dydzs_out", "dpps_out", "m2_out", "ph_out", "pathlen_out",
                          "decdist_out", "Mh2_final_out", "dflag_out", "n_unif", "rng_err"})
        c.push_back({n, 0.0, 0.0});
    for (const char* n : kDecCols) c.push_back({n, 0.0, 0.0});
    c.push_back({"adrift", 0.0, 0.0});
    c.push_back({"row", 0.0, 0.0});
    long stops = 0;
    int rc = simcl1::run("decay step (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        loadQueues(r, rng);
        optics::Track t;
        t.xs = r["xs"];
        t.dxdzs = r["dxdzs"];
        t.ys = r["ys"];
        t.dydzs = r["dydzs"];
        t.dpps = r["dpps"];
        optics::DecayState d;
        d.decay_flag = std::lround(r["decay_flag"]) != 0;
        d.dflag = std::lround(r["dflag_in"]) != 0;
        d.m2 = r["m2"];
        d.ph = r["ph"];
        d.ctau = r["ctau"];
        d.decdist = r["decdist_in"];
        d.Mh2_final = r["Mh2_final_in"];
        d.rng = &rng;
        d.fortran_compat = true;
        d.trace = &tt;
        tt.reset();
        double pathlen = r["pathlen_in"];
        const int routine = static_cast<int>(std::lround(r["routine"]));
        const int spectr = static_cast<int>(std::lround(r["spectr"]));
        const int cls = static_cast<int>(std::lround(r["class"]));
        double adrift = -1;
        try {
            if (routine == 1) {
                optics::project(t, r["zd"], pathlen, &d);
            } else {
                eng->transp(spectr, cls, t, r["zd"], pathlen, &d);
                adrift = eng->adrift(spectr, cls) ? 1. : 0.;
            }
        } catch (const FortranStop&) {
            ++stops;
        }
        o["xs_out"] = t.xs;
        o["dxdzs_out"] = t.dxdzs;
        o["ys_out"] = t.ys;
        o["dydzs_out"] = t.dydzs;
        o["dpps_out"] = t.dpps;
        o["m2_out"] = d.m2;
        o["ph_out"] = d.ph;
        o["pathlen_out"] = pathlen;
        o["decdist_out"] = d.decdist;
        o["Mh2_final_out"] = d.Mh2_final;
        o["dflag_out"] = d.dflag ? 1. : 0.;
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        fillDecay(o, tt, false);
        o["adrift"] = adrift;
        o["row"] = static_cast<double>(r.number());
    });
    if (stops) {
        std::printf("  unexpected FortranStop in %ld rows\n", stops);
        return 1;
    }
    return rc;
}

// ------------------------------------------------------------------------------------------------ P9.F1 decay tables
// decay_fixedrng (task P9.2): one project (mode 1) or transp (mode 2) call with decay, the P9.F1 column names.
int runPionDecayStep(const std::string& dir) {
    auto eng = std::make_shared<optics::CosyEngine>();
    std::ostringstream log;
    eng->setLog(&log);
    eng->setFortranCompat(true);
    eng->setCompat294f662(true);   // U.R: 294f662 reference (D-006 kaon decay mass Mk)
    const std::string mdir = hms::McHms::defaultMatrixDir();
    eng->initForward(1, mdir);
    eng->initForward(5, mdir);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<ColumnSpec> c;
    for (const char* n : {"xs_out", "dxdzs_out", "ys_out", "dydzs_out", "dpps_out", "dflag_out", "m2_out", "ph_out",
                          "pathlen_out", "decdist_out", "Mh2_final_out", "n_unif", "rng_err"})
        c.push_back({n, 0.0, 0.0});
    long stops = 0;
    int rc = simcl1::run("decay_fixedrng project/transp (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        loadQueues(r, rng);
        optics::Track t;
        t.xs = r["xs"];
        t.dxdzs = r["dxdzs"];
        t.ys = r["ys"];
        t.dydzs = r["dydzs"];
        t.dpps = r["dpps"];
        optics::DecayState d;
        d.decay_flag = std::lround(r["decay_flag"]) != 0;
        d.dflag = std::lround(r["dflag_in"]) != 0;
        d.m2 = r["m2_in"];
        d.ph = r["ph_in"];
        d.ctau = r["ctau"];
        d.decdist = r["decdist_in"];
        d.Mh2_final = r["Mh2_final_in"];
        d.rng = &rng;
        d.fortran_compat = true;
        double pathlen = r["pathlen_in"];
        try {
            if (std::lround(r["mode"]) == 1)
                optics::project(t, r["zd"], pathlen, &d);
            else
                eng->transp(static_cast<int>(std::lround(r["spectr"])), static_cast<int>(std::lround(r["class"])), t,
                            r["zd"], pathlen, &d);
        } catch (const FortranStop&) {
            ++stops;
        }
        o["xs_out"] = t.xs;
        o["dxdzs_out"] = t.dxdzs;
        o["ys_out"] = t.ys;
        o["dydzs_out"] = t.dydzs;
        o["dpps_out"] = t.dpps;
        o["dflag_out"] = d.dflag ? 1. : 0.;
        o["m2_out"] = d.m2;
        o["ph_out"] = d.ph;
        o["pathlen_out"] = pathlen;
        o["decdist_out"] = d.decdist;
        o["Mh2_final_out"] = d.Mh2_final;
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
    });
    if (stops) {
        std::printf("  unexpected FortranStop in %ld rows\n", stops);
        return 1;
    }
    return rc;
}

// mc_decay_fixedrng (task P9.2): mc_hms / mc_shms with decay for a pion (decdist = 0, Mh2_final = m2 per ray, as
// simc.f:174-177), the call counts of the harness wrappers (n_calls, n_project, n_transp, decay_call/kind/class/zd).
int runPionMcDecay(const std::string& dir) {
    Arm hmsArm(true, true), shmsArm(false, true);
    std::vector<optics::StepRecord> stepsH, stepsS;
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<ColumnSpec> c;
    for (const char* n : {"ok", "stop_code", "stop_id", "reached_hut", "x_fp", "dx_fp", "y_fp", "dy_fp", "dpp_out",
                          "dxdz_out", "dydz_out", "y_out", "pathlen", "m2_out", "resmult", "xs_end", "ys_end",
                          "dxdzs_end", "dydzs_end", "dpps_end", "n_gauss", "n_unif", "rng_err", "nsig_last", "decdist",
                          "Mh2_final", "n_calls", "n_project", "n_transp", "decay_call", "decay_kind", "decay_class",
                          "decay_zd"})
        c.push_back({n, 0.0, 0.0});
    return simcl1::run("mc_decay_fixedrng (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        loadQueues(r, rng);
        const int spectr = static_cast<int>(std::lround(r["spectr"]));
        Arm& arm = spectr == 1 ? hmsArm : shmsArm;
        CallIn in;
        in.spectr = spectr;
        in.p_spec = r["p_spec"];
        in.th_spec = r["th_spec"];
        in.dpp = r["dpp"];
        in.x = r["x"];
        in.y = r["y"];
        in.z = r["z"];
        in.dxdz = r["dxdz"];
        in.dydz = r["dydz"];
        in.m2 = r["m2"];
        in.fry = r["fry"];
        in.ms = std::lround(r["ms_flag"]) != 0;
        in.wcs = std::lround(r["wcs_flag"]) != 0;
        in.decay = std::lround(r["decay_flag"]) != 0;
        in.coll = std::lround(r["using_coll"]) != 0;
        in.ctau = r["ctau"];
        in.decdist = 0.;
        in.Mh2_final = in.m2;
        in.resmult = kUnset;
        in.pathlen = 0.;
        const CallOut x = arm.call(in, rng);
        arm.takeLog();
        o["ok"] = x.ok ? 1. : 0.;
        o["stop_code"] = x.stop_code;
        o["stop_id"] = x.stop_id;
        o["reached_hut"] = x.reached_hut ? 1. : 0.;
        o["x_fp"] = x.x_fp;
        o["dx_fp"] = x.dx_fp;
        o["y_fp"] = x.y_fp;
        o["dy_fp"] = x.dy_fp;
        o["dpp_out"] = x.dpp;
        o["dxdz_out"] = x.dxdz;
        o["dydz_out"] = x.dydz;
        o["y_out"] = x.y;
        o["pathlen"] = x.pathlen;
        o["m2_out"] = x.m2;
        o["resmult"] = x.resmult;
        o["xs_end"] = x.t.xs;
        o["ys_end"] = x.t.ys;
        o["dxdzs_end"] = x.t.dxdzs;
        o["dydzs_end"] = x.t.dydzs;
        o["dpps_end"] = x.t.dpps;
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        o["nsig_last"] = rng.lastNsig();
        o["decdist"] = x.decdist;
        o["Mh2_final"] = x.Mh2_final;
        const optics::TransportTrace& tt = arm.tt();
        long np = 0, nt = 0;
        for (const auto& st : arm.steps()) (st.kind == 1 ? np : nt)++;
        o["n_calls"] = static_cast<double>(tt.n_steps);
        o["n_project"] = static_cast<double>(np);
        o["n_transp"] = static_cast<double>(nt);
        o["decay_call"] = tt.decayed ? static_cast<double>(tt.dec_step) : 0.;
        o["decay_kind"] = tt.decayed ? (tt.dec.kind == 1 ? 1. : 2.) : 0.;
        o["decay_class"] = tt.decayed ? tt.dec.cls : 0.;
        o["decay_zd"] = tt.decayed ? tt.dec.step_len : 0.;
    });
}

// ------------------------------------------------------------------------------------------------ pion_coll_absorb
int runPca(const std::string& dir) {
    std::vector<ColumnSpec> c = {{"transmission", 0.0, 0.0}, {"row", 0.0, 0.0}};
    return simcl1::run("pion_coll_absorb (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        o["transmission"] = hms::pionCollAbsorb(r["ppi"], r["thick"]);
        o["row"] = static_cast<double>(r.number());
    });
}

// ------------------------------------------------------------------------------------------------ capture replay
struct CapCall {
    long n{0};
    std::string arm;
    std::vector<double> I, O;
    std::vector<long> K;
    std::vector<double> U, G;
};

bool readCapture(const std::string& path, std::vector<CapCall>& calls) {
    std::ifstream f(path);
    if (!f) return false;
    std::map<long, std::size_t> idx;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        std::string t;
        long n;
        if (!(ss >> t >> n)) continue;
        auto it = idx.find(n);
        if (it == idx.end()) {
            idx[n] = calls.size();
            calls.push_back(CapCall{});
            calls.back().n = n;
            it = idx.find(n);
        }
        CapCall& c = calls[it->second];
        if (t == "I") {
            ss >> c.arm;
            double x;
            while (ss >> x) c.I.push_back(x);
        } else if (t == "O") {
            double x;
            while (ss >> x) c.O.push_back(x);
        } else if (t == "K") {
            long nc, x;
            ss >> nc;
            while (ss >> x) c.K.push_back(x);
        } else if (t == "U" || t == "G") {
            long k;
            ss >> k;
            double x;
            std::vector<double>& v = (t == "U") ? c.U : c.G;
            while (ss >> x) v.push_back(x);
        }
    }
    return true;
}

// I: p_spec th_spec dpp x y z dxdz dydz m2 ms wcs decay fry coll spectr pathlen_in resmult_in decdist Mh2_final ctau
//    save_H save_S
// O: ok dpp dxdz dydz y x_fp dx_fp y_fp dy_fp m2 pathlen resmult decdist Mh2_final save_H save_S xs ys dxdzs dydzs dpps
int runReplay(const std::string& path, bool corrected) {
    std::vector<CapCall> calls;
    if (!readCapture(path, calls)) {
        std::printf("replay: SKIP (cannot read %s)\n", path.c_str());
        return 77;
    }
    Arm hmsArm(true, !corrected), shmsArm(false, !corrected);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    struct Stat {
        long calls{0}, bad{0}, okRef{0}, okCpp{0}, okChanged{0}, rngShort{0}, pathChanged{0}, decayed{0};
        long reconChanged{0}, rangedOut{0};
        double maxReconDiff[4]{0, 0, 0, 0};   // dpp (%), xptar, yptar, ytar (cm) of calls ok in both
        long slitRef[3]{0, 0, 0}, slitCpp[3]{0, 0, 0};
        std::map<int, long> stopRef, stopCpp;   // stop counter index -> count
        double maxPathDiff{0};
    };
    std::map<std::string, Stat> st;
    const char* oname[21] = {"ok",     "dpp",     "dxdz",    "dydz",      "y",      "x_fp",   "dx_fp",
                             "y_fp",   "dy_fp",   "m2",      "pathlen",   "resmult", "decdist", "Mh2_final",
                             "save_H", "save_S",  "xs",      "ys",        "dxdzs",  "dydzs",  "dpps"};
    std::map<std::string, long> badCol;
    long firstBad = 0;
    for (const CapCall& c : calls) {
        const bool isHms = c.arm == "hms";
        Arm& arm = isHms ? hmsArm : shmsArm;
        Stat& s = st[c.arm];
        ++s.calls;
        if (c.I.size() != 22 || c.O.size() != 21) {
            std::printf("replay: malformed call %ld\n", c.n);
            return 1;
        }
        rng.load(c.G, c.U);
        CallIn in;
        in.p_spec = c.I[0];
        in.th_spec = c.I[1];
        in.dpp = c.I[2];
        in.x = c.I[3];
        in.y = c.I[4];
        in.z = c.I[5];
        in.dxdz = c.I[6];
        in.dydz = c.I[7];
        in.m2 = c.I[8];
        in.ms = c.I[9] != 0;
        in.wcs = c.I[10] != 0;
        in.decay = c.I[11] != 0;
        in.fry = c.I[12];
        in.coll = c.I[13] != 0;
        in.spectr = static_cast<int>(c.I[14]);
        in.pathlen = c.I[15];
        in.resmult = c.I[16];
        in.decdist = c.I[17];
        in.Mh2_final = c.I[18];
        in.ctau = c.I[19];
        arm.setSave(isHms ? c.I[20] : c.I[21]);
        CallOut x;
        try {
            x = arm.call(in, rng);
        } catch (const std::exception& e) {
            std::printf("  call %ld: exception %s\n", c.n, e.what());
            ++s.bad;
            continue;
        }
        arm.takeLog();
        const double cv[21] = {x.ok ? 1. : 0., x.dpp,      x.dxdz, x.dydz,
                               x.y,            x.x_fp,     x.dx_fp, x.y_fp,
                               x.dy_fp,        x.m2,       x.pathlen, x.resmult,
                               x.decdist,      x.Mh2_final, isHms ? x.save_coll : c.I[20],
                               isHms ? c.I[21] : x.save_coll, x.t.xs, x.t.ys,
                               x.t.dxdzs,      x.t.dydzs,  x.t.dpps};
        // stop reason from the counter deltas
        auto stopOf = [&](const std::vector<long>& k, bool success) {
            if (success) return 0;
            const int n = static_cast<int>(k.size());
            for (int i = n; i >= 1; --i) {
                const int slot = i;
                if (isHms && (slot == hms::kTrials || slot == hms::kHut || slot == hms::kSlitHor ||
                              slot == hms::kSlitVert || slot == hms::kSlitOct))
                    continue;
                if (!isHms && (slot == shms::kTrials || slot == shms::kHut || slot == shms::kStopIdSlot ||
                               slot == shms::kSlitHor || slot == shms::kSlitVert || slot == shms::kSlitOct))
                    continue;
                if (k[static_cast<std::size_t>(i - 1)] > 0) return slot;
            }
            return -1;   // slit stop (fixed-slit path) or none
        };
        const bool okRef = c.O[0] != 0;
        s.okRef += okRef;
        s.okCpp += x.ok;
        if (okRef != x.ok) ++s.okChanged;
        if (rng.error()) ++s.rngShort;
        if (arm.tt().decayed) ++s.decayed;
        std::vector<long> kc(x.dcnt.begin() + 1, x.dcnt.end());
        ++s.stopRef[stopOf(c.K, okRef)];
        ++s.stopCpp[stopOf(kc, x.ok)];
        if (!sameBits(x.pathlen, c.O[10])) {
            ++s.pathChanged;
            s.maxPathDiff = std::max(s.maxPathDiff, std::fabs(x.pathlen - c.O[10]));
        }
        if (arm.ct().ranged_out) ++s.rangedOut;
        {
            const int sl[3] = {isHms ? int(hms::kSlitHor) : int(shms::kSlitHor),
                               isHms ? int(hms::kSlitVert) : int(shms::kSlitVert),
                               isHms ? int(hms::kSlitOct) : int(shms::kSlitOct)};
            for (int k = 0; k < 3; ++k) {
                s.slitRef[k] += c.K[static_cast<std::size_t>(sl[k] - 1)];
                s.slitCpp[k] += x.dcnt[static_cast<std::size_t>(sl[k])];
            }
        }
        if (okRef && x.ok) {
            const double dv[4] = {x.dpp - c.O[1], x.dxdz - c.O[2], x.dydz - c.O[3], x.y - c.O[4]};
            bool ch = false;
            for (int k = 0; k < 4; ++k) {
                if (dv[k] != 0) ch = true;
                s.maxReconDiff[k] = std::max(s.maxReconDiff[k], std::fabs(dv[k]));
            }
            if (ch) ++s.reconChanged;
        }
        if (corrected) continue;
        bool good = true;
        for (int k = 0; k < 21; ++k) {
            if (k >= 5 && k <= 8 && !okRef) continue;   // x_fp..dy_fp: simc keeps stale values unless ok
            if (!sameBits(cv[k], c.O[k])) {
                ++badCol[oname[k]];
                good = false;
            }
        }
        if (kc.size() != c.K.size()) {
            good = false;
            ++badCol["K size"];
        } else {
            for (std::size_t i = 0; i < kc.size(); ++i)
                if (kc[i] != c.K[i]) {
                    good = false;
                    ++badCol["counters"];
                    break;
                }
        }
        if (rng.nUniform() != c.U.size() || rng.nGauss() != c.G.size() || rng.error()) {
            good = false;
            ++badCol["rng"];
        }
        if (!good) {
            ++s.bad;
            if (!firstBad) firstBad = c.n;
        }
    }
    std::printf("replay %s (%s mode): %zu recorded mc calls\n", path.c_str(), corrected ? "CORRECTED" : "fortran_compat",
                calls.size());
    long bad = 0;
    for (const auto& kv : st) {
        const Stat& s = kv.second;
        bad += s.bad;
        std::printf("  %-5s calls %6ld  ok Fortran %6ld  ok C++ %6ld  ok changed %5ld  decayed(C++) %5ld  pathlen "
                    "differs %5ld (max %.4g cm)  RNG queue too short %ld  differing calls %ld\n",
                    kv.first.c_str(), s.calls, s.okRef, s.okCpp, s.okChanged, s.decayed, s.pathChanged, s.maxPathDiff,
                    s.rngShort, s.bad);
        std::printf("        slit_hor/vert/oct counters: Fortran %ld %ld %ld -> C++ %ld %ld %ld; ranged out in the collimator "
                    "(D-037 F4) %ld; ok in both with a different recon %ld (max |diff| delta %.3g %%, xptar %.3g, yptar "
                    "%.3g, ytar %.3g cm)\n",
                    s.slitRef[0], s.slitRef[1], s.slitRef[2], s.slitCpp[0], s.slitCpp[1], s.slitCpp[2], s.rangedOut,
                    s.reconChanged, s.maxReconDiff[0], s.maxReconDiff[1], s.maxReconDiff[2], s.maxReconDiff[3]);
        std::printf("        stop reasons (counter: Fortran -> C++):");
        std::map<int, int> keys;
        for (const auto& a : s.stopRef) keys[a.first] = 1;
        for (const auto& a : s.stopCpp) keys[a.first] = 1;
        for (const auto& k : keys) {
            const char* nm = k.first == 0 ? "success"
                             : k.first < 0 ? "slit/none"
                                           : (kv.first == "hms" ? kHmsNames[k.first] : kShmsNames[k.first]);
            auto get = [](const std::map<int, long>& m, int key) {
                auto it = m.find(key);
                return it == m.end() ? 0L : it->second;
            };
            std::printf(" %s %ld->%ld", nm, get(s.stopRef, k.first), get(s.stopCpp, k.first));
        }
        std::printf("\n");
    }
    if (corrected) return 0;
    for (const auto& kv : badCol) std::printf("  column %s differs in %ld calls\n", kv.first.c_str(), kv.second);
    if (firstBad) std::printf("  first differing call: %ld\n", firstBad);
    std::printf("  => %s (%ld differing calls)\n", bad ? "FAIL" : "PASS", bad);
    return bad ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s mc <dir> hms|shms | tp <dir> | pca <dir> | replay <cap_calls.txt> [corrected]\n",
                     argv[0]);
        return 2;
    }
    const std::string mode = argv[1], dir = argv[2];
    try {
        if (mode == "mc" && argc >= 4) return runMc(dir, std::string(argv[3]) == "hms");
        if (mode == "tp") return runTp(dir);
        if (mode == "pca") return runPca(dir);
        if (mode == "pdecay") return runPionDecayStep(dir);   // P9.2: decay_fixedrng
        if (mode == "mcdecay") return runPionMcDecay(dir);    // P9.2: mc_decay_fixedrng
        if (mode == "replay") return runReplay(dir, argc >= 4 && std::string(argv[3]) == "corrected");
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode\n");
    return 2;
}
