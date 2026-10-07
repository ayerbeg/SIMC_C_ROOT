// src/driver/Central.cpp — see include/simc/driver/Central.h.
// Fortran: simc.f:calculate_central (l.1143-1309).
#include "simc/driver/Central.h"

#include <cmath>
#include <stdexcept>

#include "simc/core/Errors.h"

namespace simc {
namespace driver {

CentralResult calculateCentral(const RunSetup& s, IReaction& reaction, IRadiator& radiator, IPhysicsModel& model) {
    CentralResult res;
    MainEvent& main0 = res.main0;   // SAVEd local: zero except what is set below
    Event& vertex0 = res.vertex0;
    CentralEvent& central = res.central;

    // l.1166-1177
    main0.target.x = 0.0;
    main0.target.y = 0.0;
    main0.target.z = s.targ.zoffset;
    main0.target.rastery = 0.0;
    main0.target.Eloss[0] = s.targ.Eloss[0].ave;
    main0.target.Coulomb = s.targ.Coulomb.ave;
    main0.target.teff[0] = s.targ.teff[0].ave;
    vertex0.Ein = s.Ebeam_vertex_ave;
    main0.Ein_shift = 0.0;
    main0.Ee_shift = 0.0;
    main0.gen_weight = 1.0;

    // l.1181-1195
    vertex0.e.delta = 0.0;
    vertex0.e.yptar = 0.0;
    vertex0.e.xptar = 0.0;
    vertex0.e.theta = s.spec.e.theta;
    vertex0.e.phi = s.spec.e.phi;
    vertex0.e.P = s.spec.e.P * (1. + vertex0.e.delta / 100.);
    vertex0.e.E = vertex0.e.P;

    vertex0.p.delta = 0.0;
    vertex0.p.yptar = 0.0;
    vertex0.p.xptar = 0.0;
    vertex0.p.theta = s.spec.p.theta;
    vertex0.p.phi = s.spec.p.phi;
    vertex0.p.P = s.spec.p.P * (1. + vertex0.p.delta / 100.);
    vertex0.p.E = std::sqrt(vertex0.p.P * vertex0.p.P + s.Mh2);

    // l.1197-1204: /pfermi/ = 0, efer = M - sqrt(m_spec**2 + pfer**2) (= M - Mrec). COMMON /pfermi/ is
    // MainEvent::pfermi; the central peepi call (pion reaction) reads it.
    main0.pfermi.pfer = 0.0;
    main0.pfermi.pferx = 0.0;
    main0.pfermi.pfery = 0.0;
    main0.pfermi.pferz = 0.0;
    vertex0.Em = s.targ.Mtar_struck + s.targ.Mrec - s.targ.M;
    {
        const double m_spec = s.targ.M - s.targ.Mtar_struck + vertex0.Em;   // =targ.Mrec
        main0.pfermi.efer = s.targ.M - std::sqrt(m_spec * m_spec + main0.pfermi.pfer * main0.pfermi.pfer);
    }

    // l.1240-1242
    if (!reaction.completeRecon(vertex0)) throw FortranStop("COMPLETE_EV failed trying to complete a CENTRAL event! (simc.f:1240-1242)");
    central.Q2 = vertex0.Q2;
    central.q = vertex0.q;
    central.nu = vertex0.nu;
    central.Em = vertex0.Em;
    central.Pm = vertex0.Pm;
    central.W = vertex0.W;
    central.e.xptar = vertex0.e.xptar;
    central.e.yptar = vertex0.e.yptar;
    central.e.delta = vertex0.e.delta;
    central.p.xptar = vertex0.p.xptar;
    central.p.yptar = vertex0.p.yptar;
    central.p.delta = vertex0.p.delta;

    // l.1256-1260 (x**2 -> x*x)
    const double mm2 = central.Em * central.Em - central.Pm * central.Pm;
    if (mm2 < 0)
        central.MM = -std::sqrt(std::fabs(mm2));
    else
        central.MM = std::sqrt(mm2);

    // l.1262-1266
    main0.target.teff[1] = s.targ.teff[1].ave;
    main0.target.teff[2] = s.targ.teff[2].ave;
    radiator.initEvent(main0, vertex0);
    // l.1268-1281
    if (s.sw.using_rad) central.rad = radiator.currentValues();

    // l.1288-1290: recon0 = vertex0; complete_main(.true., main0, vertex0, vertex0, recon0, success)
    const Event recon0 = vertex0;
    reaction.prepareCentral(main0, vertex0);   // P9.2 hook (pion: D-038 F5 in corrected mode); no-op otherwise
    model.beginEvent(0);   // nevent = 0 before the loop
    reaction.completeMain(true, main0, vertex0, vertex0, recon0);
    central.sigcc = main0.sigcc;
    return res;
}

}  // namespace driver
}  // namespace simc
