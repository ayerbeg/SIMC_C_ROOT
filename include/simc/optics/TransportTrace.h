// include/simc/optics/TransportTrace.h
//
// Optional diagnostic records of the transport kernels (task P6.4). They are filled only when a caller attaches them
// (optics::DecayState::trace, hms::McHms / shms::McShms ::setCollTrace) and never change the numerics or the random
// numbers. They mirror what the Fortran harness of P6.F2 records with its ld --wrap recorders around project, transp,
// loren, pion_coll_absorb and enerloss_new (reference/functions/mc_*_coll, mc_*_decay, transp_decay_*, project_decay:
// columns n_steps, n_decay_checks, dec_*, n_absorb_calls, ab_*, n_eloss_calls, eloss_sum; steps.txt.xz, coll_calls.txt).
#ifndef SIMC_OPTICS_TRANSPORTTRACE_H
#define SIMC_OPTICS_TRANSPORTTRACE_H

#include <vector>

namespace simc {
namespace optics {

/// One decay (project.f:63-116, transp.f:142-185 first half, transp.f:234-275 second half).
struct DecayRecord {
    int kind{0};             ///< 1 project, 2 transp 1st half (before the COSY map), 3 transp 2nd half (after it)
    int cls{0};              ///< COSY class (transp), 0 for project
    double x{0}, y{0};       ///< COMMON xs, ys when loren is called
    double dxdz_before{0}, dydz_before{0}, dpp_before{0}, p_before{0};   ///< parent slopes, dpps, ph
    double pathlen_at{0};    ///< pathlen when loren is called
    double decdist_at{0};    ///< decdist after the += z_decay
    double m2_before{0};
    double gamma{0}, bx{0}, by{0}, bz{0}, er{0}, pxr{0}, pyr{0}, pzr{0};   ///< loren inputs
    double ef{0}, pxf{0}, pyf{0}, pzf{0}, pf{0};                         ///< loren outputs
    double dxdz_after{0}, dydz_after{0};   ///< pxf/pzf, pyf/pzf
    double dpp_after{0}, m2_after{0}, p_after{0};   ///< dpps, m2, ph at the end of the decaying step
    double step_len{0};      ///< zd (transp) / z_drift (project)
    double z_decay{0};       ///< the sampled decay distance (cm along the path)
};

/// Entry or exit state of one project/transp call.
struct StepState {
    double xs{0}, ys{0}, dxdzs{0}, dydzs{0}, dpps{0}, m2{0}, ph{0}, pathlen{0}, decdist{0};
};

/// One project/transp call (steps.txt.xz of the P6.F2 tables).
struct StepRecord {
    int kind{0};   ///< 1 project, 2 transp
    int cls{0};    ///< transp class, 0 for project
    double zd{0};
    bool decay_check{false};   ///< entered with decay_flag .and. .not.dflag (one grnd)
    bool dflag_in{false}, dflag_out{false};
    StepState in, out;
};

/// Per-call counters of the transport (project/transp) kernels.
struct TransportTrace {
    long n_steps{0};          ///< every project/transp call
    long n_decay_checks{0};   ///< calls entered with decay_flag .and. .not.dflag
    bool decayed{false};
    long dec_check{0};        ///< 1-based index of the decaying check among the decay checks
    long dec_step{0};         ///< 1-based index of the decaying step among all steps
    DecayRecord dec;
    std::vector<StepRecord>* steps{nullptr};   ///< optional: every step
    void reset() {
        std::vector<StepRecord>* s = steps;
        *this = TransportTrace{};
        steps = s;
        if (steps) steps->clear();
    }
};

/// One call inside the collimator stepping (coll_calls.txt of the P6.F2 tables).
struct CollCall {
    int type{0};   ///< 1 pion_coll_absorb (a = p, thick, transmission, u, absorbed); 2 enerloss_new (len, epart, m, Eloss)
    double a[5]{0, 0, 0, 0, 0};
};

/// Per-call record of mc_hms_coll / mc_shms_coll.
struct CollTrace {
    long n_absorb_calls{0};
    bool absorbed{false};
    double ab_p_last{0}, ab_thick_last{0}, ab_trans_last{0}, ab_u_last{0};
    long n_eloss_calls{0};
    double eloss_sum{0};
    bool ranged_out{false};   ///< corrected mode only: stopped because enerloss_new clamped the loss (D-037 F4)
    std::vector<CollCall>* calls{nullptr};
    void reset() {
        std::vector<CollCall>* c = calls;
        *this = CollTrace{};
        calls = c;
        if (calls) calls->clear();
    }
};

}  // namespace optics
}  // namespace simc

#endif  // SIMC_OPTICS_TRANSPORTTRACE_H
