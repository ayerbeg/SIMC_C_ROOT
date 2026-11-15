/**
 * @file
 * @ingroup stages
 * @brief simc::IRadiator: the radiative-corrections stage, and the values it reports.
 * @note Origin: Fortran SIMC init.f:radc_init (l.576), radc_init_ev (l.655), radc.f:generate_rad (l.120),
 *       peaked_rad_weight (l.523).
 */
#ifndef SIMC_STAGES_IRADIATOR_H
#define SIMC_STAGES_IRADIATOR_H

#include <array>
#include <functional>

#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief The event-completion callback that the radiator calls after the incoming electron has radiated.
 *
 * Signature `bool(MainEvent& main, Event& vertex, IRandom& rng)`: the reaction's IReaction::completeEv. After the
 * incoming-electron photon has reduced vertex.Ein, the callee recomputes the vertex kinematics and the target quantities
 * teff[1..2] / Eloss[1..2] (with Gaussian draws from `rng`) and must call IRadiator::initEvent(main, vertex) on the
 * SAME radiator, exactly like the original program. It returns true on success.
 * @note Origin: Fortran SIMC event.f:complete_ev(main, vertex, success), as radc.f:generate_rad calls it
 *       (radc.f:322-328; event.f:1045 radc_init_ev).
 */
using CompleteEvFn = std::function<bool(MainEvent& main, Event& vertex, IRandom& rng)>;

/**
 * @brief The radiation quantities of the current event that the driver reports (the `.hist` "Central Values").
 *
 * The three tails (1 incoming electron, 2 scattered electron, 3 hadron) are indexed 0..2; c_int/c_ext/c/g keep the
 * original (0:3) indexing (index 0 = total).
 * @note Origin: Fortran SIMC radc.inc /radccom/, copied into central%rad by calculate_central (simc.f:1268-1281).
 */
struct RadiationValues {
    /// Hard-photon correction factor; eta*zeta of the target material; internal and external b parameters (g).
    double hardcorfac{0}, etatzai{0}, g_int{0}, g_ext{0};
    std::array<double, 3> frac{};     ///< frac(1:3): fraction of the radiation probability per tail
    std::array<double, 3> lambda{};   ///< lambda(1:3): radiation exponent of each tail
    std::array<double, 2> bt{};       ///< bt(1:2): b x t of the target material (incoming, outgoing)
    std::array<double, 4> c_int{}, c_ext{}, c{}, g{};   ///< (0:3): internal, external and total normalisations c, exponents g
};

/**
 * @brief Radiative corrections: generates the bremsstrahlung photons of a try (vertex -> orig) and weights the event.
 *
 * Per try: the reaction calls initEvent() from its event completion, then radiate(). The radiator keeps per-event
 * state between the two (and the reported values); one instance per run / worker, not thread-safe. Implementation:
 * simc::radiation::SimcRadiator ("simc_rad").
 */
class IRadiator : public IStage {
public:
    /**
     * @brief Per-event initialisation from the vertex kinematics (radiation probabilities of the three tails).
     * @param main The quantities of the try (target thicknesses are read).
     * @param vertex The completed vertex kinematics.
     * @note Origin: Fortran SIMC init.f:radc_init_ev.
     */
    virtual void initEvent(MainEvent& main, const Event& vertex) = 0;
    /**
     * @brief Radiates: generates the photon energies of the tails and builds the orig event from the vertex event.
     *
     * main.gen_weight is multiplied by the radiative weight. @p vertex is modified like in the original program
     * (vertex.Ein reduced by the incoming-electron photon and recompleted through @p complete_ev); the caller keeps
     * using it as the vertex event.
     * @param main The quantities of the try (gen_weight and the ntuple radiation fields are written).
     * @param vertex In/out: the vertex kinematics.
     * @param orig Out: the kinematics after radiation (written only on success, D-024).
     * @param rng The random source.
     * @param complete_ev The reaction's event completion (see simc::CompleteEvFn).
     * @return False if no valid radiation could be generated (the try fails).
     * @note Origin: Fortran SIMC radc.f:generate_rad.
     */
    virtual bool radiate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng,
                         const CompleteEvFn& complete_ev) = 0;

    // ---- what the driver reads back (D-022 revisited in P7.4; defaults suit a radiator without these notions) ----
    /**
     * @brief The radiation values after the last initEvent() (the driver takes them right after the central event's
     *        initEvent()).
     * @return The values; default: all zero.
     */
    virtual RadiationValues currentValues() const { return {}; }
    /**
     * @brief The photon energy radiated in each tail by the last successful radiate().
     * @return [MeV] per tail, index 0..2 = incoming electron, scattered electron, hadron; default: zeros.
     * @note Origin: Fortran SIMC Egamma_used(1:3), read by limits_update.
     */
    virtual std::array<double, 3> photonEnergies() const { return {}; }
    /**
     * @brief Tells whether the hadron radiates in the current event.
     *
     * Such successes are not counted in RunSummary::ncontribute_no_rad_proton (simc.f:296).
     * @return True if the hadron tail is on for this event; default: false.
     * @note Origin: Fortran SIMC rad_proton_this_ev, set by radc_init_ev.
     */
    virtual bool hadronRadiates() const { return false; }
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_IRADIATOR_H
