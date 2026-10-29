// include/simc/core/random/RanluxRandom.h
//
// IRandom "ranlux" (task P13.2): the random-number generator of the Fortran SIMC, ported bit for bit, so that a SIMC++ run
// can reproduce the random stream of a Fortran run exactly (bitwise replays of whole runs; D-002 keeps mt64 the default).
//   cern/ranlux.f: RANLUX (Luescher's luxury subtract-with-borrow, F. James 1993), RLUXGO (initialisation from one integer),
//                  RLUXIN / RLUXUT (the 25-integer state); compiled with -fdefault-real-8 (REAL = double; every value is a
//                  multiple of 2^-48, so the arithmetic is exact in double precision).
//   call_ranlux.f: sgrnd(seed) = RLUXGO(3, seed, 0, 0); grnd() hands out RVEC(1..999) of each block of 1000 numbers and
//                  refills when latest reaches 1000 (the 1000th number of every block is never used).
//   gauss1.f:      gauss1 on this uniform stream (gauss1Polar).
// JSON config (optional): { "seed": <integer> } (default: the .inp random_seed, simc.f:107 sgrnd(random_seed)), or
// { "state": [25 integers] } (simc.f:109 start_file_random_state -> RLUXIN; the *_start_random_state.dat content).
#ifndef SIMC_CORE_RANDOM_RANLUXRANDOM_H
#define SIMC_CORE_RANDOM_RANLUXRANDOM_H

#include <array>
#include <cstdint>
#include <string>

#include "simc/stages/IRandom.h"

namespace simc {

class RanluxRandom : public IRandom {
public:
    RanluxRandom();

    std::string name() const override { return "ranlux"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext&) override {}

    /// call_ranlux.f:grnd
    double uniform() override;
    /// gauss1.f on grnd()
    double gauss1(double nsig_max) override;
    /// call_ranlux.f:sgrnd(seed) = RLUXGO(3, seed, 0, 0); also resets grnd's block pointer as a fresh process.
    void seed(std::uint64_t s) override;

    /// ranlux.f:RLUXGO(lux, ins, k1, k2)
    void rluxgo(int lux, int ins, int k1, int k2);
    /// ranlux.f:RLUXIN(isdext)
    void rluxin(const std::array<int, 25>& isdext);
    /// ranlux.f:RLUXUT(isdext)
    std::array<int, 25> rluxut() const;
    /// ranlux.f:RANLUX(rvec, lenv)
    void ranlux(double* rvec, int lenv);

private:
    double seeds_[25]{};   // SEEDS(1..24)
    int next_[25]{};       // NEXT(1..24)
    int i24_{24}, j24_{10}, in24_{0}, kount_{0}, mkount_{0}, luxlev_{3}, nskip_{199}, inseed_{0};
    double carry_{0}, twom24_{0}, twom12_{0};
    bool notyet_{true};
    // call_ranlux.f:grnd SAVEd locals (static zero under -fno-automatic)
    double rvec_[1001]{};
    int latest_{0};
};

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_RANLUXRANDOM_H
