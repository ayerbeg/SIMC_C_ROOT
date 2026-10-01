// include/simc/reactions/StandardBeam.h
//
// IBeam "standard": the beam/vertex part of event.f:generate (l.150-208), task P4.1.
//   main%target%x/y   gaussian beam width (gauss1(3), +-3 sigma) + offsets            (l.156-157)
//   raster            fr_pattern 1 (bedpost), 2 (circle), 3 (flat), else none          (l.163-181)
//   main%target%z     ITarget::generateZ (= (0.5-grnd())*targ%length + targ%zoffset)  (l.185)
//   Eloss(1), teff(1) trip_thru_target(1, z-zoffset, Ebeam, 0, ..., Me, typeflag 1)  (l.193-195)
//   Coulomb           targ%Coulomb_constant if using_Coulomb                           (l.196-202)
//   vertex%Ein        Ebeam + (grnd()-0.5)*dEbeam + Coulomb - Eloss(1)                  (l.204-205)
//   Ein_shift/Ee_shift                                                                  (l.211-212)
// RNG order (D-024): 2 gauss1, 2 grnd (raster 1/2/3), 1 grnd (z), the target gaussians, 1 grnd (Ein, drawn even for
// dEbeam = 0). The target gaussians are drawn even with using_Eloss = 0 (Eloss(1) is zeroed afterwards).
// Collaborator: the run's ITarget, taken from the StageContext at init() (the driver owns it).
#ifndef SIMC_REACTIONS_STANDARDBEAM_H
#define SIMC_REACTIONS_STANDARDBEAM_H

#include <string>

#include "simc/core/RunSetup.h"
#include "simc/stages/IBeam.h"
#include "simc/stages/ITarget.h"

namespace simc {
namespace reactions {

class StandardBeam : public IBeam {
public:
    std::string name() const override { return "standard"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    /// Fortran: event.f:generate l.156-212.
    void generate(MainEvent& main, Event& vertex, IRandom& rng) override;

private:
    ITarget* target_{nullptr};   ///< not owned (StageContext::target)
    // run constants (RunSetup)
    double xwid_{0}, ywid_{0}, xoffset_{0}, yoffset_{0}, fr1_{0}, fr2_{0}, zoffset_{0};
    int fr_pattern_{0};
    double Ebeam_{0}, dEbeam_{0}, Ebeam_vertex_ave_{0};
    bool using_Eloss_{false}, using_Coulomb_{false};
    double Coulomb_constant_{0}, Coulomb_ave_{0};
};

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_STANDARDBEAM_H
