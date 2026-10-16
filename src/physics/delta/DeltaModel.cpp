// src/physics/delta/DeltaModel.cpp — physics_delta.f:peedelta and the IPhysicsModel "delta" (task P11.1).
// See include/simc/physics/delta/DeltaModel.h.
//
// Numerics: as physics_delta.f compiled with -fdefault-real-8 (D-010): double literals except the DOUBLE PRECISION
// literals 1.d-8 (l.74) and 1.d6 (l.99), which are REAL(16) (binary128, __float128 here); x**2 -> x*x; pi**2 -> pi*pi.
#include "simc/physics/delta/DeltaModel.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/physics/pion/PionCrossSections.h"
#include "simc/reactions/EventChain.h"

namespace simc {
namespace delta {

namespace {
__extension__ typedef __float128 quad;
inline double sq(double x) { return x * x; }
using constants::pi;
}  // namespace

double binary128LowBytes(double x) {
    static_assert(sizeof(quad) == 16, "binary128 expected");
    const quad q = static_cast<quad>(x) / static_cast<quad>(1.e6);   // x/1.d6: __extenddftf2 + __divtf3
    double d;
    std::memcpy(&d, &q, sizeof d);   // the REAL*8 dummy reads the first 8 bytes (x86-64: low mantissa bits)
    return d;
}

// Fortran: physics_delta.f:peedelta (l.1-135)
double peedelta(const Event& vertex, MainEvent& main, const DeltaConfig& cfg, std::ostream* log, PeedeltaTrace* trace) {
    const FermiMomentum& f = main.pfermi;   // COMMON /pfermi/

    // l.50-62: photon-NUCLEON C.M. frame with the proton as the 'hadron' (jacobians.f)
    const reactions::CmTransform cm =
        reactions::transformToCm(vertex, main, f.pfer, f.pferx, f.pfery, f.pferz, f.efer, log);
    main.thetacm = cm.thetacm;
    main.phicm = cm.phicm;
    main.pcm = cm.phadcm;
    main.davejac = cm.jacobian;
    main.johnjac = cm.jac_old;   // approx. assuming collinear boost.

    // l.73-78. FORTRAN-QUIRK: `tfcos-1..lt.1.d-8` compares in binary128 (1.d-8 is a REAL(16) literal, __lttf2).
    double tfcos = f.pferx * vertex.uq.x + f.pfery * vertex.uq.y + f.pferz * vertex.uq.z;
    if (tfcos - 1. > 0. && static_cast<quad>(tfcos - 1.) < static_cast<quad>(1) / static_cast<quad>(100000000))
        tfcos = 1.0;
    const double tfsin = std::sqrt(1. - sq(tfcos));
    const double s = sq(vertex.nu + f.efer) - sq(vertex.q + f.pfer * tfcos) - sq(f.pfer * tfsin);
    main.wcm = std::sqrt(s);

    // l.91 (JRA 2001)
    const double k_eq = (sq(main.wcm) - sq(cfg.Mtar_struck)) / 2. / cfg.Mtar_struck;

    // l.99-103: the Blok fit (written to the ntuple only)
    double t_arg, q2_arg, s_arg;
    if (cfg.sigblok_binary128_quirk) {
        // fortran_compat: what the REAL*8 dummies of sig_blok read (D-048 F3)
        t_arg = binary128LowBytes(main.t);
        q2_arg = binary128LowBytes(vertex.Q2);
        s_arg = binary128LowBytes(s);
    } else {
        // FORTRAN-BUG(D-048 F3): original Fortran physics_delta.f:99-100:
        // `sig_blok(thetacm,phicm,main%t/1.d6,vertex%q2/1.d6,s/1.d6,...)`: 1.d6 is REAL(16) under -fdefault-real-8, so
        // binary128 temporaries are passed to REAL*8 dummies and sig_blok reads their low 8 bytes; corrected because
        // the intended arguments are t, Q2 and s in GeV^2: the doubles are passed (ntup%sigcm becomes the Blok fit).
        t_arg = main.t / 1.e6;
        q2_arg = vertex.Q2 / 1.e6;
        s_arg = s / 1.e6;
    }
    main.ntup.sigcm1 = pion::sigBlok(cm.thetacm, cm.phicm, t_arg, q2_arg, s_arg, main.epsilon, cfg.Mtar_struck / 1000.,
                                     cfg.which_pion);
    const double sigma_eepi = main.ntup.sigcm1;
    main.ntup.sigcm = sigma_eepi;   // sig_cm

    // l.127-128
    const double fac = 1. / (1. - f.pferz * f.pfer / f.efer) * cfg.Mtar_struck / f.efer;
    const double gtpr =
        constants::alpha / 2. / (pi * pi) * vertex.e.E / vertex.Ein * k_eq / vertex.Q2 / (1. - main.epsilon);

    if (trace) {
        trace->tfcos = tfcos;
        trace->s = s;
        trace->k_eq = k_eq;
        trace->fac = fac;
        trace->gtpr = gtpr;
        trace->t_arg = t_arg;
        trace->q2_arg = q2_arg;
        trace->s_arg = s_arg;
    }

    // l.130-132. FORTRAN-DIFF-CANDIDATE(D-048 F4): `peedelta = sigma_eepi*jacobian*(gtpr*fac)` is commented out; the
    // Fortran returns 1.0*jacobian*(gtpr*fac): the Delta weight contains NO cross section (phase space x photon flux x
    // the dt dphi_cm -> dOmega_lab jacobian). Question for the user/an expert: which Delta cross section should be used?
    // Kept in both modes.
    return 1.0 * cm.jacobian * (gtpr * fac);
}

// ------------------------------------------------------------------------------------------------- plug-in
void DeltaModel::configure(const Json& config, const SimcParameters&) {
    if (config.is_object() && config.contains("verbose")) verbose_ = config.at("verbose").get<bool>();
}

void DeltaModel::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    cfg_.which_pion = s.params.which_pion;
    cfg_.Mtar_struck = s.targ.Mtar_struck;
    cfg_.sigblok_binary128_quirk = s.fortran_compat;   // D-048 F3
    if (!(cfg_.Mtar_struck > 0))
        throw std::invalid_argument("delta: targ%Mtar_struck is not set (the \"delta\" model needs a doing_delta set-up)");
}

double DeltaModel::sigcc(const Event& vertex, const MainEvent& main) {
    MainEvent m = main;
    return peedelta(vertex, m, cfg_, verbose_ ? &std::cout : nullptr);
}

double DeltaModel::sigccMain(const Event& vertex, MainEvent& main) {
    return peedelta(vertex, main, cfg_, verbose_ ? &std::cout : nullptr);
}

}  // namespace delta
}  // namespace simc

SIMC_REGISTER(simc::IPhysicsModel, "delta", simc::delta::DeltaModel)
