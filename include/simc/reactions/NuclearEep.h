// include/simc/reactions/NuclearEep.h
//
// IReaction "deuterium_eep" (D(e,e'p), doing_deuterium) and "nucleus_eep" (A(e,e'p), doing_heavy) — task P8.2.
//   generate()        event.f:generate l.214-428 after the beam part: electron angles, hadron angles, [A > 2: hadron
//                     energy], electron energy (gen_weight = generated fraction of each energy range), physics_angles,
//                     /pfermi/ reset (no Fermi momentum is generated for (e,e'p): Pm and Em follow from the e and p
//                     kinematics), complete_ev, main%Trec, generate_rad or (norad) the A > 2 VERTEXedge Em/Pm test and
//                     orig = vertex. Random numbers: D 4+1 grnd, A > 2 4+2 grnd (+ the beam part and complete_ev).
//   spectralWeight()  event.f:complete_main l.1397-1430 through the ISpectralFunction component ("theory": init.f theory
//                     files; "benhar": Benhar tables with use_benhar_sf; ISpectralFunction.h)
//   completeEv(), completeRecon(), completeMain(), genVolume(): EepReactionBase (D and heavy branches of event.f)
//
// JSON config (all optional):
//   {"spectral_function": {"type": "theory" | "benhar", "file": "<name>", "transparency": x}, "verbose": bool}
// The type defaults to the input (use_benhar_sf), the file to the Fortran choice for the target A (RunSetup::sf, with
// a warning when the Fortran falls back to C-12, D-033 F1). Both names create this class; the branch is chosen by the
// RunSetup flags, as in the Fortran (theory_init sets doing_deuterium from the file contents, init.f:889).
#ifndef SIMC_REACTIONS_NUCLEAREEP_H
#define SIMC_REACTIONS_NUCLEAREEP_H

#include <memory>
#include <string>

#include "simc/reactions/EepReactionBase.h"
#include "simc/core/SimcParameters.h"
#include "simc/stages/ISpectralFunction.h"

namespace simc {
namespace reactions {

class NuclearEep : public EepReactionBase {
public:
    explicit NuclearEep(std::string name = "nucleus_eep");

    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) override;
    double spectralWeight(const Event& vertex) override;

    ISpectralFunction& spectralFunction() { return *sf_; }

private:
    Json sfConfig_;
    SimcParameters params_;   ///< for the spectral function's configure()
    std::unique_ptr<ISpectralFunction> sf_;
};

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_NUCLEAREEP_H
