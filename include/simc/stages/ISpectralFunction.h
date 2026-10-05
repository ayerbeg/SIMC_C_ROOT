// include/simc/stages/ISpectralFunction.h
// Spectral function of a nuclear (e,e'p) reaction (task P8.1): the probability density of removing a proton with
// missing momentum Pm and missing energy Em. It is a component of the A(e,e'p) reaction plug-in ("deuterium_eep",
// "nucleus_eep"), selected in the reaction's JSON config:
//     "reaction": {"type": "nucleus_eep", "spectral_function": {"type": "theory" | "benhar", "file": "<name>"}}
// Fortran: event.f:complete_main l.1397-1430 (main%SF_weight), init.f:theory_init, sf_lookup.f.
// Built-in plug-ins: "theory" (the init.f theory files: h2/c12/fe56/au197.theory), "benhar" (Benhar tables,
// use_benhar_sf). The file is read in init() from RunSetup::data_dir.
#ifndef SIMC_STAGES_ISPECTRALFUNCTION_H
#define SIMC_STAGES_ISPECTRALFUNCTION_H

#include <string>

#include "simc/stages/IStage.h"

namespace simc {

class ISpectralFunction : public IStage {
public:
    /// main%SF_weight for the vertex (Em, Pm) of an event [the units of the Fortran: (MeV/c)^-3 for D,
    /// (MeV/c)^-3 MeV^-1 for A > 2, times the number of protons]. Not const: a Fortran-faithful implementation may
    /// keep SAVEd state between calls (sf_lookup, fortran_compat).
    virtual double weight(const Event& vertex) = 0;
    /// The file that was read (for messages and reports).
    virtual std::string file() const = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_ISPECTRALFUNCTION_H
