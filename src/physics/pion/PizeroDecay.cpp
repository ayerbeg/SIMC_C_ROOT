// src/physics/pion/PizeroDecay.cpp — Fortran: pizero_decay.f (l.1-94). See PizeroDecay.h.
// pi0 -> 2 gamma: an isotropic decay in the pi0 rest frame, boosted to the lab; the photons are then tracked to the
// calorimeter (driver/PizeroCalorimeter.cpp) or stored in the ntuple.
#include "simc/physics/pion/PizeroDecay.h"

#include <cmath>

#include "simc/core/SimcConstants.h"
#include "simc/setup/KinematicsUtils.h"

namespace simc {
namespace pion {

// Fortran: pizero_decay.f:pizero_decay (l.1-94)
PizeroDecayResult pizeroDecay(const Event& vertex, double Mh, IRandom& rng, std::ostream* log) {
    using constants::pi;
    const double Mgamma = 0.0;   // parameter(Mgamma=0.0)
    PizeroDecayResult r;

    // beta and gamma of the pi0 rest frame (l.37-41)
    const double ph = vertex.p.P;
    const double eh = std::sqrt(ph * ph + Mh * Mh);
    const double beta = ph / eh;
    const double gamma = 1. / std::sqrt(1. - beta * beta);

    // decay angles in the rest frame (l.44-46): 2 grnd, in this order
    const double rph = rng.uniform() * 2. * pi;
    const double rth1 = rng.uniform() * 2. - 1.;
    const double rth = std::acos(rth1);

    const double er = Mh / 2.0;
    if (er < Mgamma) {
        r.success = false;
        return r;
    }
    const double pr = std::sqrt(er * er - Mgamma * Mgamma);
    const double pxr1 = pr * std::sin(rth) * std::cos(rph);
    const double pyr1 = pr * std::sin(rth) * std::sin(rph);
    const double pzr1 = pr * std::cos(rth);
    const double pxr2 = -pxr1;
    const double pyr2 = -pyr1;
    const double pzr2 = -pzr1;

    // boost to the lab (l.66-74)
    const double bx = -beta * vertex.up.x;
    const double by = -beta * vertex.up.y;
    const double bz = -beta * vertex.up.z;
    LorenResult L = Loren(gamma, bx, by, bz, er, pxr1, pyr1, pzr1);
    r.gamma1 = {L.ef1, L.pxf, L.pyf, L.pzf};
    L = Loren(gamma, bx, by, bz, er, pxr2, pyr2, pzr2);
    r.gamma2 = {L.ef1, L.pxf, L.pyf, L.pzf};

    // energy and momentum conservation checks (l.79-93)
    if (log) {
        if (std::fabs(eh - r.gamma1[0] - r.gamma2[0]) > 0.001) *log << " pi0 decay: bad energy sum\n";
        if (std::fabs(ph * vertex.up.x - r.gamma1[1] - r.gamma2[1]) > 0.001) *log << " pi0 decay: bad px sum\n";
        if (std::fabs(ph * vertex.up.y - r.gamma1[2] - r.gamma2[2]) > 0.001) *log << " pi0 decay: bad py sum\n";
        if (std::fabs(ph * vertex.up.z - r.gamma1[3] - r.gamma2[3]) > 0.001) *log << " pi0 decay: bad pz sum\n";
    }
    r.success = true;
    return r;
}

}  // namespace pion
}  // namespace simc
