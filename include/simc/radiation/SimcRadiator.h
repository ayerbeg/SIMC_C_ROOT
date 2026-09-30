// include/simc/radiation/SimcRadiator.h
//
// IRadiator "simc_rad": the SIMC radiative corrections (radc.f, brem.f, init.f radc_init_ev/basicrad_init_ev),
// task P5.1. One instance = the radiation COMMON blocks + SAVEd locals of one Fortran process (RadcState).
//
//   init(setup)          radc_init values from RunSetup (BuildRunSetup applies init.f:radc_init/limits_init)
//   initEvent(main, v)   init.f:radc_init_ev (l.655-728), called by complete_ev (event.f:1045) and calculate_central
//   radiate(...)         radc.f:generate_rad (l.120-519), with complete_ev supplied by the caller (CompleteEvFn)
//   peakedRadWeight(...) radc.f:peaked_rad_weight (l.523-646)
// The free functions (Radc.h, Brem.h) carry the numerics; this class only owns the state.
//
// JSON config (optional): { "verbose": bool }  (Fortran write(6,*) messages to std::cout; default true)
#ifndef SIMC_RADIATION_SIMCRADIATOR_H
#define SIMC_RADIATION_SIMCRADIATOR_H

#include <iosfwd>
#include <string>

#include "simc/radiation/Radc.h"
#include "simc/radiation/RadcState.h"
#include "simc/stages/IRadiator.h"

namespace simc {
namespace radiation {

class SimcRadiator : public IRadiator {
public:
    SimcRadiator();

    std::string name() const override { return "simc_rad"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(RunSetup& setup) override;

    /// Fortran: init.f:radc_init_ev(main, vertex).
    void initEvent(MainEvent& main, const Event& vertex) override;
    /// Fortran: radc.f:generate_rad(main, vertex, orig, success).
    bool radiate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng, const CompleteEvFn& complete_ev) override;
    /// Fortran: radc.f:peaked_rad_weight(vertex, Egamma, emin, emax, basicrad_val_reciprocal, basicrad_weight).
    double peakedRadWeight(const Event& vertex, double Egamma, double emin, double emax, double basicrad_val_reciprocal,
                           double basicrad_weight);

    /// Set the run constants directly (L1 tests, or a caller without a RunSetup).
    void setRun(const RadcRun& run) { run_ = run; }
    const RadcRun& run() const { return run_; }
    /// The radiation COMMON of the current event (ntail, Egamma_used, rad_proton_this_ev, hardcorfac, frac, ...).
    RadcState& state() { return st_; }
    const RadcState& state() const { return st_; }
    /// Where the Fortran write(6,*) messages go (nullptr = silent).
    void setLog(std::ostream* os) { log_ = os; }

private:
    RadcRun run_;
    RadcState st_;
    std::ostream* log_;
};

}  // namespace radiation
}  // namespace simc

#endif  // SIMC_RADIATION_SIMCRADIATOR_H
