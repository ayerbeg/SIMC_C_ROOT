/**
 * @file
 * @ingroup stages
 * @brief simc::ISpectralFunction: the spectral function of a nuclear (e,e'p) reaction.
 *
 * The probability density of removing a proton with missing momentum Pm and missing energy Em. It is a component of
 * the A(e,e'p) reaction plug-ins ("deuterium_eep", "nucleus_eep"), selected in the reaction's JSON configuration:
 * @code{.json}
 * "reaction": {"type": "nucleus_eep", "spectral_function": {"type": "theory" | "benhar", "file": "<name>"}}
 * @endcode
 * Built-in plug-ins: "theory" (the theory files h2/c12/fe56/au197.theory), "benhar" (Benhar tables, use_benhar_sf).
 * The file is read in init() from RunSetup::data_dir. Task P8.1.
 * @note Origin: Fortran SIMC event.f:complete_main l.1397-1430 (main%SF_weight), init.f:theory_init, sf_lookup.f.
 */
#ifndef SIMC_STAGES_ISPECTRALFUNCTION_H
#define SIMC_STAGES_ISPECTRALFUNCTION_H

#include <string>

#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Spectral function S(Em, Pm) of a nucleus, evaluated for the vertex of a try (see the file description).
 *
 * Implementations: simc::nuclear::TheorySpectralFunction ("theory"), simc::nuclear::BenharSpectralFunction
 * ("benhar"). The tables are read once in init() and then only read; an implementation may keep state between calls
 * (not thread-safe).
 */
class ISpectralFunction : public IStage {
public:
    /**
     * @brief The spectral-function weight main.SF_weight for the vertex (Em, Pm) of a try.
     *
     * Not const: a faithful implementation may keep saved state between calls (sf_lookup in fortran_compat mode).
     * @param vertex The vertex kinematics (Em [MeV], Pm [MeV/c]).
     * @return The weight, in the units of the original program: (MeV/c)^-3 for deuterium, (MeV/c)^-3 MeV^-1 for
     *         A > 2, times the number of protons; 0 outside the table.
     */
    virtual double weight(const Event& vertex) = 0;
    /**
     * @brief The table file that was read.
     * @return The file name (for messages and reports).
     */
    virtual std::string file() const = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_ISPECTRALFUNCTION_H
