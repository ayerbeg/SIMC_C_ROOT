// src/driver/TargetFieldAngles.cpp — see include/simc/driver/TargetFieldAngles.h.
// using_tgt_field: the angles between the polarized-target field (targ_Bangle, targ_Bphi) and each spectrometer arm,
// as simc.f computes them before trgInit, with the log lines the Fortran prints. Used by MonteCarlo::initTargetField.
#include "simc/driver/TargetFieldAngles.h"

#include <cmath>

#include "simc/core/SimcConstants.h"
#include "simc/material/FortranFormat.h"

namespace simc {
namespace driver {

// Fortran: simc.f:120-151
ArmFieldAngles TargetFieldArmAngles(double targ_Bangle, double targ_Bphi, double e_theta, double e_phi, double p_theta,
                         double p_phi, bool compat) {
    using constants::degrad;
    using fortran::listReal;
    ArmFieldAngles r;
    auto opposite = [&](double dphi) {
        // FORTRAN-BUG(D-067 T1): original Fortran simc.f:127/141 `degrad*abs(targ_Bphi-phi)-180.0 .lt. .01` has no
        // outer abs, so every |dphi| < 180.01 deg that is not "same side" counts as "opposite side"; corrected: the
        // outer abs. fortran_compat keeps the Fortran test.
        return compat ? (degrad * std::fabs(dphi) - 180.0 < .01) : (std::fabs(degrad * std::fabs(dphi) - 180.0) < .01);
    };
    if (degrad * std::fabs(targ_Bphi - e_phi) < .01) {
        if (targ_Bangle >= e_theta)
            r.earm = -1 * std::sin(e_phi) * (targ_Bangle - e_theta);
        else
            r.earm = +1 * std::sin(e_phi) * (e_theta - targ_Bangle);
    } else if (opposite(targ_Bphi - e_phi)) {
        r.earm = +1 * std::sin(e_phi) * (targ_Bangle + e_theta);
    } else {
        r.log.push_back("  Error determining angle between target and e-arm");
        r.log.push_back("  targ_Bangle =" + listReal(targ_Bangle * degrad) + " targ_Bphi =" + listReal(targ_Bphi * degrad));
        r.log.push_back("  Central e-arm angle =" + listReal(e_theta * degrad) + " e-arm phi =" + listReal(e_phi * degrad));
    }
    if (degrad * std::fabs(targ_Bphi - p_phi) < .01) {
        if (targ_Bangle >= p_theta)
            r.parm = -1 * std::sin(p_phi) * (targ_Bangle - p_theta);
        else if (compat)
            // FORTRAN-BUG(D-067 T2): original Fortran simc.f:139 `+1*sin(phi)*(targ_Bangle-theta)` = -sin(phi)*|B-theta|,
            // the sign of the B >= theta branch, unlike the e arm (simc.f:125); corrected: the e-arm formula.
            r.parm = +1 * std::sin(p_phi) * (targ_Bangle - p_theta);
        else
            r.parm = +1 * std::sin(p_phi) * (p_theta - targ_Bangle);
    } else if (opposite(targ_Bphi - p_phi)) {
        r.parm = +1 * std::sin(p_phi) * (targ_Bangle + p_theta);
    } else {
        r.log.push_back("  Error determining angle between target and p-arm");
        r.log.push_back("  targ_Bangle =" + listReal(targ_Bangle * degrad) + " targ_Bphi =" + listReal(targ_Bphi * degrad));
        // FORTRAN-QUIRK(D-067 T3): simc.f:146 prints spec%p%theta as the "p-arm phi"
        r.log.push_back("  Central p-arm angle =" + listReal(p_theta * degrad) + " p-arm phi =" +
                        listReal(p_theta * degrad));
    }
    r.log.push_back("  targ_Bangle =" + listReal(targ_Bangle * degrad) + " targ_Bphi =" + listReal(targ_Bphi * degrad));
    r.log.push_back("  Angle between target and e-arm" + listReal(r.earm * degrad));
    r.log.push_back("  Angle between target and p-arm" + listReal(r.parm * degrad));
    return r;
}

}  // namespace driver
}  // namespace simc
