// src/core/random/RanluxRandom.cpp — see include/simc/core/random/RanluxRandom.h.
// Bit-for-bit port of the CERNLIB RANLUX generator as SIMC calls it, so that a RANLUX run with the seed of the .inp
// reproduces the Fortran event by event (the "ranlux" IRandom; the shipped-input regression tests rely on it).
// The arrays keep the Fortran 1-based indexing (seeds_[1..24], next_[1..24]); element 0 is unused.
// Fortran: cern/ranlux.f (RANLUX, RLUXIN, RLUXUT, RLUXGO), call_ranlux.f (sgrnd, grnd), gauss1.f.
#include "simc/core/random/RanluxRandom.h"

#include <cstdlib>
#include <stdexcept>

#include "simc/core/Registry.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/Gauss1.h"

namespace simc {

namespace {
// ranlux.f PARAMETERs/DATA: luxury levels 0..4 and the number of numbers skipped per 24 at each level
constexpr int kMaxlev = 4;
constexpr int kNdskip[kMaxlev + 1] = {0, 24, 73, 199, 365};
constexpr int kIgiga = 1000000000, kJsdflt = 314159265, kItwo24 = 1 << 24, kIcons = 2147483563;
constexpr double kTwop12 = 4096.;
}  // namespace

RanluxRandom::RanluxRandom() = default;

// Configuration: an explicit 25-word RLUXUT state, else a seed (JSON, else the .inp random_seed, as simc.f:107).
void RanluxRandom::configure(const Json& config, const SimcParameters& params) {
    if (config.is_object() && config.contains("state")) {
        const auto& v = config.at("state");
        if (!v.is_array() || v.size() != 25) throw std::invalid_argument("rng ranlux: \"state\" needs 25 integers");
        std::array<int, 25> s{};
        for (std::size_t i = 0; i < 25; ++i) s[i] = v[i].get<int>();
        rluxin(s);
        latest_ = 0;
    } else if (config.is_object() && config.contains("seed")) {
        seed(config.at("seed").get<std::uint64_t>());
    } else {
        seed(static_cast<std::uint64_t>(params.random_seed));
    }
}

// call_ranlux.f:sgrnd: luxury level 3, no skipping; latest_ = 0 forces a fresh block at the next grnd
void RanluxRandom::seed(std::uint64_t s) {
    rluxgo(3, static_cast<int>(s), 0, 0);
    latest_ = 0;
}

// call_ranlux.f:grnd (l.58-72): rvec(1..999) of each RANLUX(rvec, 1000) block. Number 1000 of each block is never
// returned (the Fortran refills when latest reaches 1000): kept, it decides the stream.
double RanluxRandom::uniform() {
    if (latest_ <= 0 || latest_ >= 1000) {
        ranlux(rvec_ + 1, 1000);
        latest_ = 1;
    }
    const double g = rvec_[latest_];
    latest_ = latest_ + 1;
    return g;
}

double RanluxRandom::gauss1(double nsig_max) { return gauss1Polar(*this, nsig_max); }

// ranlux.f:RANLUX (incl. the default initialisation of the first call)
void RanluxRandom::ranlux(double* rvec, int lenv) {
    // --- first call without RLUXGO/RLUXIN: default initialisation (seed 314159265, level 3)
    if (notyet_) {
        notyet_ = false;
        int jseed = kJsdflt;
        inseed_ = jseed;
        luxlev_ = 3;
        nskip_ = kNdskip[luxlev_];
        in24_ = 0;
        kount_ = 0;
        mkount_ = 0;
        twom24_ = 1.;
        int iseeds[25];
        // the 24 seeds from a multiplicative congruential generator (L'Ecuyer, Schrage's factorisation avoids overflow)
        for (int i = 1; i <= 24; ++i) {
            twom24_ = twom24_ * 0.5;
            const int k = jseed / 53668;
            jseed = 40014 * (jseed - k * 53668) - k * 12211;
            if (jseed < 0) jseed = jseed + kIcons;
            iseeds[i] = jseed % kItwo24;
        }
        twom12_ = twom24_ * 4096.;
        for (int i = 1; i <= 24; ++i) {
            seeds_[i] = static_cast<double>(iseeds[i]) * twom24_;
            next_[i] = i - 1;   // circular list of the lags
        }
        next_[1] = 24;
        i24_ = 24;
        j24_ = 10;
        carry_ = 0.;
        if (seeds_[24] == 0.) carry_ = twom24_;
    }
    // --- subtract-with-borrow x(n) = x(n-10) - x(n-24) - carry, mod 1
    for (int ivec = 0; ivec < lenv; ++ivec) {
        double uni = seeds_[j24_] - seeds_[i24_] - carry_;
        if (uni < 0.) {
            uni = uni + 1.0;
            carry_ = twom24_;
        } else {
            carry_ = 0.;
        }
        seeds_[i24_] = uni;
        i24_ = next_[i24_];
        j24_ = next_[j24_];
        rvec[ivec] = uni;
        if (uni < twom12_) {   // small numbers are padded
            rvec[ivec] = rvec[ivec] + twom24_ * seeds_[j24_];
            if (rvec[ivec] == 0.) rvec[ivec] = twom24_ * twom24_;
        }
        in24_ = in24_ + 1;
        if (in24_ == 24) {   // skipping to luxury
            in24_ = 0;
            kount_ = kount_ + nskip_;
            for (int isk = 1; isk <= nskip_; ++isk) {
                uni = seeds_[j24_] - seeds_[i24_] - carry_;
                if (uni < 0.) {
                    uni = uni + 1.0;
                    carry_ = twom24_;
                } else {
                    carry_ = 0.;
                }
                seeds_[i24_] = uni;
                i24_ = next_[i24_];
                j24_ = next_[j24_];
            }
        }
    }
    kount_ = kount_ + lenv;
    if (kount_ >= kIgiga) {
        mkount_ = mkount_ + 1;
        kount_ = kount_ - kIgiga;
    }
}

// ranlux.f:RLUXIN: restore a state written by RLUXUT (24 seeds as 24-bit integers; word 25 packs i24, j24, in24,
// the luxury level and the carry as its sign)
void RanluxRandom::rluxin(const std::array<int, 25>& isdext) {
    notyet_ = false;
    twom24_ = 1.;
    for (int i = 1; i <= 24; ++i) {
        next_[i] = i - 1;
        twom24_ = twom24_ * 0.5;
    }
    next_[1] = 24;
    twom12_ = twom24_ * 4096.;
    for (int i = 1; i <= 24; ++i) seeds_[i] = static_cast<double>(isdext[static_cast<std::size_t>(i - 1)]) * twom24_;
    carry_ = 0.;
    if (isdext[24] < 0) carry_ = twom24_;
    int isd = std::abs(isdext[24]);
    i24_ = isd % 100;
    isd = isd / 100;
    j24_ = isd % 100;
    isd = isd / 100;
    in24_ = isd % 100;
    isd = isd / 100;
    luxlev_ = isd;
    if (luxlev_ <= kMaxlev) {
        nskip_ = kNdskip[luxlev_];
    } else if (luxlev_ >= 24) {
        nskip_ = luxlev_ - 24;
    } else {
        nskip_ = kNdskip[kMaxlev];
        luxlev_ = kMaxlev;
    }
    inseed_ = -1;
}

// ranlux.f:RLUXUT: the inverse of RLUXIN
std::array<int, 25> RanluxRandom::rluxut() const {
    std::array<int, 25> s{};
    for (int i = 1; i <= 24; ++i) s[static_cast<std::size_t>(i - 1)] = static_cast<int>(seeds_[i] * kTwop12 * kTwop12);
    s[24] = i24_ + 100 * j24_ + 10000 * in24_ + 1000000 * luxlev_;
    if (carry_ > 0.) s[24] = -s[24];
    return s;
}

// ranlux.f:RLUXGO: initialise with luxury level lux and seed ins, then skip k1 + k2*1e9 numbers
void RanluxRandom::rluxgo(int lux, int ins, int k1, int k2) {
    if (lux < 0) {
        luxlev_ = 3;
    } else if (lux <= kMaxlev) {
        luxlev_ = lux;
    } else if (lux < 24 || lux > 2000) {
        luxlev_ = kMaxlev;
    } else {
        luxlev_ = lux;
        for (int ilx = 0; ilx <= kMaxlev; ++ilx)
            if (lux == kNdskip[ilx] + 24) luxlev_ = ilx;
    }
    if (luxlev_ <= kMaxlev)
        nskip_ = kNdskip[luxlev_];
    else
        nskip_ = luxlev_ - 24;
    in24_ = 0;
    int jseed = ins > 0 ? ins : kJsdflt;
    inseed_ = jseed;
    notyet_ = false;
    twom24_ = 1.;
    int iseeds[25];
    for (int i = 1; i <= 24; ++i) {
        twom24_ = twom24_ * 0.5;
        const int k = jseed / 53668;
        jseed = 40014 * (jseed - k * 53668) - k * 12211;
        if (jseed < 0) jseed = jseed + kIcons;
        iseeds[i] = jseed % kItwo24;
    }
    twom12_ = twom24_ * 4096.;
    for (int i = 1; i <= 24; ++i) {
        seeds_[i] = static_cast<double>(iseeds[i]) * twom24_;
        next_[i] = i - 1;
    }
    next_[1] = 24;
    i24_ = 24;
    j24_ = 10;
    carry_ = 0.;
    if (seeds_[24] == 0.) carry_ = twom24_;
    kount_ = k1;
    mkount_ = k2;
    // skip k1 + k2*1e9 numbers (no luxury skipping during the advance; in24 is then recomputed)
    if (k1 + k2 != 0) {
        for (int iouter = 1; iouter <= k2 + 1; ++iouter) {
            const int inner = (iouter == k2 + 1) ? k1 : kIgiga;
            for (int isk = 1; isk <= inner; ++isk) {
                double uni = seeds_[j24_] - seeds_[i24_] - carry_;
                if (uni < 0.) {
                    uni = uni + 1.0;
                    carry_ = twom24_;
                } else {
                    carry_ = 0.;
                }
                seeds_[i24_] = uni;
                i24_ = next_[i24_];
                j24_ = next_[j24_];
            }
        }
        in24_ = kount_ % (nskip_ + 24);
        if (mkount_ > 0) {
            const int izip = kIgiga % (nskip_ + 24);
            const int izip2 = mkount_ * izip + in24_;
            in24_ = izip2 % (nskip_ + 24);
        }
        if (in24_ > 23) in24_ = 0;
    }
}

}  // namespace simc

SIMC_REGISTER(simc::IRandom, "ranlux", simc::RanluxRandom)
