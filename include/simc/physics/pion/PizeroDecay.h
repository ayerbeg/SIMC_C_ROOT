// include/simc/physics/pion/PizeroDecay.h
//
// Fortran: pizero_decay.f:pizero_decay (l.1-94), task P9.1. The pi0 of the vertex decays into two photons, isotropic in
// the pi0 rest frame, boosted to the lab with loren.f (the single Loren port, setup/KinematicsUtils.h). Called by
// event.f:complete_ev (l.899-901) when doing_pizero; the photons go to ntup%gamma1/gamma2 (E, px, py, pz in the lab).
//
// Random numbers: exactly 2 uniforms per call, in this order: rph = 2*pi*u1 (grnd()*2.*pi), cos(theta*) = 2*u2 - 1.
#ifndef SIMC_PHYSICS_PION_PIZERODECAY_H
#define SIMC_PHYSICS_PION_PIZERODECAY_H

#include <array>
#include <ostream>

#include "simc/core/Event.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace pion {

struct PizeroDecayResult {
    bool success{false};
    std::array<double, 4> gamma1{};   ///< Fortran gamma1(1:4) = (E, px, py, pz) [MeV], lab
    std::array<double, 4> gamma2{};
};

/// Fortran: pizero_decay.f (l.1-94). Reads vertex%p%P and vertex%up; Mh = COMMON Mh (the pi0 mass for doing_pizero).
/// The conservation checks (0.001 MeV) print 'pi0 decay: bad ... sum' to `log` (if non-null), as the Fortran.
/// success = false only if Mh/2 < Mgamma (= 0): the Fortran then returns without writing gamma1/gamma2 (the caller
/// keeps its previous values; here they are zero and must not be copied).
PizeroDecayResult pizeroDecay(const Event& vertex, double Mh, IRandom& rng, std::ostream* log = nullptr);

}  // namespace pion
}  // namespace simc

#endif  // SIMC_PHYSICS_PION_PIZERODECAY_H
