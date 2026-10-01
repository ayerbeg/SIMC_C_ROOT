// src/reactions/StandardBeam.cpp — see include/simc/reactions/StandardBeam.h.
#include "simc/reactions/StandardBeam.h"

#include <cmath>
#include <stdexcept>

#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace reactions {

void StandardBeam::configure(const Json&, const SimcParameters&) {}

void StandardBeam::init(RunSetup& s) {
    xwid_ = s.gen.xwid;
    ywid_ = s.gen.ywid;
    xoffset_ = s.targ.xoffset;
    yoffset_ = s.targ.yoffset;
    fr1_ = s.targ.fr1;
    fr2_ = s.targ.fr2;
    fr_pattern_ = s.targ.fr_pattern;
    zoffset_ = s.targ.zoffset;
    Ebeam_ = s.Ebeam;
    dEbeam_ = s.dEbeam;
    Ebeam_vertex_ave_ = s.Ebeam_vertex_ave;
    using_Eloss_ = s.sw.using_Eloss;
    using_Coulomb_ = s.sw.using_Coulomb;
    Coulomb_constant_ = s.targ.Coulomb_constant;
    Coulomb_ave_ = s.targ.Coulomb.ave;
}

// Fortran: event.f:generate (l.156-212)
void StandardBeam::generate(MainEvent& main, Event& vertex, IRandom& rng) {
    using constants::pi;
    if (!target_) throw std::logic_error("StandardBeam: no target (setTarget)");
    constexpr double nsig_max = 3.0e0;   // l.140-141

    main.target.x = rng.gauss1(nsig_max) * xwid_ + xoffset_;   // l.156-157
    main.target.y = rng.gauss1(nsig_max) * ywid_ + yoffset_;

    double t3, t4, t5, t6;   // l.163-181
    if (fr_pattern_ == 1) {   // old bedpost, square raster
        t3 = rng.uniform() * pi;
        t4 = rng.uniform() * pi;
        t5 = std::cos(t3) * fr1_;
        t6 = std::cos(t4) * fr2_;
    } else if (fr_pattern_ == 2) {   // circular raster
        t3 = rng.uniform() * 2. * pi;
        t4 = std::sqrt(rng.uniform()) * (fr2_ - fr1_) + fr1_;
        t5 = std::cos(t3) * t4;
        t6 = std::sin(t3) * t4;
    } else if (fr_pattern_ == 3) {   // new, flat square raster
        t3 = 2. * rng.uniform() - 1.0;
        t4 = 2. * rng.uniform() - 1.0;
        t5 = fr1_ * t3;
        t6 = fr2_ * t4;
    } else {
        t5 = 0.0;
        t6 = 0.0;
    }

    main.target.x = main.target.x + t5;   // l.183-187
    main.target.y = main.target.y + t6;
    main.target.z = target_->generateZ(rng);   // (0.5-grnd())*targ%length+targ%zoffset
    main.target.rastery = t6;
    main.target.rasterx = t5;

    // l.193-205
    target_->tripThruTarget(ITarget::kBeam, main.target.z - zoffset_, Ebeam_, 0.0e0, constants::Me, ITarget::kSample,
                            rng, main.target.Eloss[0], main.target.teff[0]);
    if (!using_Eloss_) main.target.Eloss[0] = 0.0;
    if (using_Coulomb_)
        main.target.Coulomb = Coulomb_constant_;
    else
        main.target.Coulomb = 0.0;
    // FORTRAN-QUIRK(D-024): the grnd is drawn even when dEbeam = 0.
    vertex.Ein = Ebeam_ + (rng.uniform() - 0.5) * dEbeam_ + main.target.Coulomb - main.target.Eloss[0];

    main.Ein_shift = vertex.Ein - Ebeam_vertex_ave_;   // l.211-212
    main.Ee_shift = main.target.Coulomb - Coulomb_ave_;
}

SIMC_REGISTER(IBeam, "standard", StandardBeam)

}  // namespace reactions
}  // namespace simc
