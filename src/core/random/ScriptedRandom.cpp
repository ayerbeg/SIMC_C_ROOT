// src/core/random/ScriptedRandom.cpp — see include/simc/core/random/ScriptedRandom.h.
#include "simc/core/random/ScriptedRandom.h"

#include <limits>
#include <utility>

#include "simc/core/random/Gauss1.h"

namespace simc {

// Fortran harness: rng_queue.f:hrng_load
void ScriptedRandom::load(std::vector<double> gaussians, std::vector<double> uniforms) {
    g_ = std::move(gaussians);
    u_ = std::move(uniforms);
    ig_ = iu_ = ncg_ = ncu_ = 0;
    nsig_last_ = 0;
    err_ = false;
}

// Fortran harness: rng_stub_grnd.f:grnd
double ScriptedRandom::uniform() {
    ++ncu_;
    if (iu_ >= u_.size()) {
        err_ = true;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return u_[iu_++];
}

// Fortran harness: rng_stub_gauss.f:gauss1 (Scripted) or gauss1.f on the uniform queue (Algorithm)
double ScriptedRandom::gauss1(double nsig_max) {
    ++ncg_;
    nsig_last_ = nsig_max;
    if (mode_ == GaussMode::Algorithm) return gauss1Polar(*this, nsig_max);
    if (ig_ >= g_.size()) {
        err_ = true;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return g_[ig_++];
}

}  // namespace simc
