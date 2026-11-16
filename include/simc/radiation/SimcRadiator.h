/**
 * @file
 * @ingroup radiation
 * @brief simc::radiation::SimcRadiator: IRadiator "simc_rad", the SIMC radiative corrections (radc.f, init.f:radc_init_ev).
 *
 * IRadiator "simc_rad": the SIMC radiative corrections (radc.f, brem.f, init.f radc_init_ev/basicrad_init_ev),
 * task P5.1. One instance = the radiation COMMON blocks + SAVEd locals of one Fortran process (RadcState).
 *
 *   init(setup)          radc_init values from RunSetup (BuildRunSetup applies init.f:radc_init/limits_init)
 *   initEvent(main, v)   init.f:radc_init_ev (l.655-728), called by complete_ev (event.f:1045) and calculate_central
 *   radiate(...)         radc.f:generate_rad (l.120-519), with complete_ev supplied by the caller (CompleteEvFn)
 *   peakedRadWeight(...) radc.f:peaked_rad_weight (l.523-646)
 * The free functions (Radc.h, Brem.h) carry the numerics; this class only owns the state.
 *
 * JSON config (optional): { "verbose": bool }  (Fortran write(6,*) messages to std::cout; default true)
 */
#ifndef SIMC_RADIATION_SIMCRADIATOR_H
#define SIMC_RADIATION_SIMCRADIATOR_H

#include <iosfwd>
#include <memory>
#include <string>

#include "simc/stages/IRadiator.h"

namespace simc {
namespace radiation {

/// @addtogroup radiation
/// @{

struct RadcRun;     // kernel (kernel/RadcState.h): run constants of radc.inc
struct RadcState;   // kernel (kernel/RadcState.h): /radccom/ + SAVEd locals

/// @brief IRadiator "simc_rad": internal and external bremsstrahlung tails of the three particles (radc.f).
///
/// Owns the run constants and the per-event radiation COMMON (RadcRun, RadcState). Not copyable, not thread-safe
/// (one radiator per worker).
class SimcRadiator : public IRadiator {
public:
    /// Creates an unconfigured radiator.
    SimcRadiator();
    /// Destroys the kernel state.
    ~SimcRadiator() override;
    /// Not copyable (owns Fortran COMMON state).
    SimcRadiator(const SimcRadiator&) = delete;
    /// Not copyable.
    /// @return (deleted)
    SimcRadiator& operator=(const SimcRadiator&) = delete;

    std::string name() const override { return "simc_rad"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;

    /// Fortran: init.f:radc_init_ev(main, vertex).
    /// @param main the main event record (teff read, radiation fields written).
    /// @param vertex the vertex.
    void initEvent(MainEvent& main, const Event& vertex) override;
    /// Fortran: radc.f:generate_rad(main, vertex, orig, success).
    /// @param main the main event record.
    /// @param vertex the generated vertex (radiated photon energies applied).
    /// @param orig the event after radiation (written on success).
    /// @param rng the random generator.
    /// @param complete_ev the reaction's complete_ev, called on the radiated kinematics.
    /// @return the Fortran success flag.
    bool radiate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng, const CompleteEvFn& complete_ev) override;
    /// Fortran: radc.f:peaked_rad_weight(vertex, Egamma, emin, emax, basicrad_val_reciprocal, basicrad_weight).
    /// @param vertex the vertex.
    /// @param Egamma photon energy [MeV].
    /// @param emin lower photon energy limit [MeV].
    /// @param emax upper photon energy limit [MeV].
    /// @param basicrad_val_reciprocal reciprocal of the basic radiation value.
    /// @param basicrad_weight basic radiation weight.
    /// @return the peaked radiation weight.
    double peakedRadWeight(const Event& vertex, double Egamma, double emin, double emax, double basicrad_val_reciprocal,
                           double basicrad_weight);

    RadiationValues currentValues() const override;
    std::array<double, 3> photonEnergies() const override;
    bool hadronRadiates() const override;

    // ---- kernel access (L1 tests, diagnostics; include kernel/RadcState.h to use them) ------------------------
    /// Set the run constants directly (L1 tests, or a caller without a RunSetup).
    /// @param run the run constants (copied).
    void setRun(const RadcRun& run);
    /// @return the run constants.
    const RadcRun& run() const { return *run_; }
    /// The radiation COMMON of the current event (ntail, Egamma_used, rad_proton_this_ev, hardcorfac, frac, ...).
    /// @return the radiation state (modifiable).
    RadcState& state() { return *st_; }
    /// @return the radiation state.
    const RadcState& state() const { return *st_; }
    /// Where the Fortran write(6,*) messages go (nullptr = silent).
    /// @param os the stream (not owned), or nullptr.
    void setLog(std::ostream* os) { log_ = os; }

private:
    std::unique_ptr<RadcRun> run_;
    std::unique_ptr<RadcState> st_;
    std::ostream* log_;
};

/// @}

}  // namespace radiation
}  // namespace simc

#endif  // SIMC_RADIATION_SIMCRADIATOR_H
