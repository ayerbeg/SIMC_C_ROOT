// src/output/NtupleLayout.cpp
// Fortran: NtupleInit.f, results_write.f:results_ntu_write (l.1-271).
#include "simc/output/NtupleLayout.h"
#include "simc/core/Units.h"

#include <cmath>

namespace simc {

namespace {
constexpr double kMp = 938.27231;   // constants.inc Mp
}

const std::vector<NtupleColumn>& EepNtupleColumns() {
    // Columns 1-12 = right-side arm ("hs"), 13-24 = left-side arm ("ss").
    // "R" below = recon of the right-side particle (recon.e if electron is right, else recon.p).
    static const std::vector<NtupleColumn> cols = {
        {"hsdelta", "recon.R.delta", "%", "results_write.f:65/91"},
        {"hsyptar", "recon.R.yptar", "slope (rad)", "results_write.f:66/92"},
        {"hsxptar", "recon.R.xptar", "slope (rad)", "results_write.f:67/93"},
        {"hsytar", "recon.R.z", "cm", "results_write.f:68/94"},
        {"hsxfp", "main.FP.R.x", "cm", "results_write.f:69/95"},
        {"hsxpfp", "main.FP.R.dx", "slope (rad)", "results_write.f:70/96"},
        {"hsyfp", "main.FP.R.y", "cm", "results_write.f:71/97"},
        {"hsypfp", "main.FP.R.dy", "slope (rad)", "results_write.f:72/98"},
        {"hsdeltai", "orig.R.delta", "%", "results_write.f:73/102"},
        {"hsyptari", "orig.R.yptar", "slope (rad)", "results_write.f:74/103"},
        {"hsxptari", "orig.R.xptar", "slope (rad)", "results_write.f:75/104"},
        {"hsytari", "main.target.z * sin(theta_R)", "cm", "results_write.f:76/105"},
        {"ssdelta", "recon.L.delta", "%", "results_write.f:77/106"},
        {"ssyptar", "recon.L.yptar", "slope (rad)", "results_write.f:78/107"},
        {"ssxptar", "recon.L.xptar", "slope (rad)", "results_write.f:79/108"},
        {"ssytar", "recon.L.z", "cm", "results_write.f:80/109"},
        {"ssxfp", "main.FP.L.x", "cm", "results_write.f:81/110"},
        {"ssxpfp", "main.FP.L.dx", "slope (rad)", "results_write.f:82/111"},
        {"ssyfp", "main.FP.L.y", "cm", "results_write.f:83/112"},
        {"ssypfp", "main.FP.L.dy", "slope (rad)", "results_write.f:84/113"},
        {"ssdeltai", "orig.L.delta", "%", "results_write.f:85/114"},
        {"ssyptari", "orig.L.yptar", "slope (rad)", "results_write.f:86/115"},
        {"ssxptari", "orig.L.xptar", "slope (rad)", "results_write.f:87/116"},
        {"ssytari", "-main.target.z * sin(theta_L)", "cm", "results_write.f:88/117"},
        {"q", "recon.q/1000", "GeV/c", "results_write.f:121"},
        {"nu", "recon.nu/1000", "GeV", "results_write.f:122"},
        {"Q2", "recon.Q2/1e6", "(GeV/c)^2", "results_write.f:123"},
        {"W", "recon.W/1000", "GeV", "results_write.f:124"},
        {"epsilon", "recon.epsilon", "", "results_write.f:125"},
        {"Em", "recon.Em/1000", "GeV", "results_write.f:126"},
        {"Pm", "recon.Pm/1000", "GeV/c", "results_write.f:127"},
        {"thetapq", "recon.theta_pq", "rad", "results_write.f:128"},
        {"phipq", "recon.phi_pq", "rad", "results_write.f:129"},
        {"corrsing", "(recon.e.P - Mp*Ebeam/(2*Ebeam*sin^2(recon.e.theta/2)+Mp))/1000", "GeV",
         "results_write.f:55-56,247"},
        {"Pmx", "-recon.Pmx/1000", "GeV/c", "results_write.f:61,248"},
        {"Pmy", "(recon.Pmz*uq.y - recon.Pmy*uq.z)/sqrt(uq.y^2+uq.z^2)/1000", "GeV/c",
         "results_write.f:59-60,249"},
        {"Pmz", "-(recon.Pmy*uq.y + recon.Pmz*uq.z)/sqrt(uq.y^2+uq.z^2)/1000", "GeV/c",
         "results_write.f:57-58,250"},
        {"PmPar", "recon.PmPar/1000", "GeV/c", "results_write.f:251"},
        {"PmPer", "recon.PmPer/1000", "GeV/c", "results_write.f:252"},
        {"PmOop", "recon.PmOop/1000", "GeV/c", "results_write.f:253"},
        {"fry", "-main.target.rastery", "cm", "results_write.f:254"},
        {"radphot", "main.ntup.radphot/1000", "GeV", "results_write.f:255"},
        {"sigcc", "main.sigcc", "ub/sr (H) or ub*MeV^2/sr^2", "results_write.f:256"},
        {"Weight", "main.weight", "", "results_write.f:257"},
        {"theta_e", "recon.e.theta", "rad", "results_write.f:258"},
        {"theta_p", "recon.p.theta", "rad", "results_write.f:259"},
    };
    return cols;
}

ArmSide ElectronArmSide(int electron_arm) {
    if (electron_arm == 1 || electron_arm == 3 || electron_arm == 7) return ArmSide::ElectronRight;
    if (electron_arm == 2 || electron_arm == 4 || electron_arm == 5 || electron_arm == 6 || electron_arm == 8)
        return ArmSide::ElectronLeft;
    return ArmSide::Unknown;
}

bool UsesEepLayout(const ReactionFlags& f) {
    // NtupleInit.f: the pion/kaon/delta branch comes first, then semi/rho, then (e,e'p).
    if (f.doing_pion || f.doing_kaon || f.doing_delta) return false;
    if (f.doing_semi || f.doing_rho) return false;
    return f.doing_hyd_elast || f.doing_deuterium || f.doing_heavy;
}

bool FillEepNtuple(const RunSetup& s, const MainEvent& main, const Event& /*vertex*/, const Event& orig,
                   const Event& recon, std::array<double, kEepNtupleSize>& ntu) {
    const double Ebeam = s.params.Ebeam;
    // results_write.f:54-62
    const double poftheta = kMp * Ebeam / (2 * Ebeam * std::pow(std::sin(recon.e.theta / 2.), 2) + kMp);
    const double corrsing = recon.e.P - poftheta;
    const double uqyz = std::sqrt(recon.uq.y * recon.uq.y + recon.uq.z * recon.uq.z);
    const double Pm_Heepz = -(recon.Pmy * recon.uq.y + recon.Pmz * recon.uq.z) / uqyz;
    const double Pm_Heepy = (recon.Pmz * recon.uq.y - recon.Pmy * recon.uq.z) / uqyz;
    const double Pm_Heepx = -recon.Pmx;

    const ArmSide side = ElectronArmSide(s.params.electron_arm);
    auto fillArm = [&](int off, const ArmFull& r, const ArmFP& fp, const ArmFull& o, double ytari) {
        ntu[off + 0] = r.delta;
        ntu[off + 1] = r.yptar;
        ntu[off + 2] = r.xptar;
        ntu[off + 3] = r.z;
        ntu[off + 4] = fp.x;
        ntu[off + 5] = fp.dx;
        ntu[off + 6] = fp.y;
        ntu[off + 7] = fp.dy;
        ntu[off + 8] = o.delta;
        ntu[off + 9] = o.yptar;
        ntu[off + 10] = o.xptar;
        ntu[off + 11] = ytari;
    };
    if (side == ArmSide::ElectronRight) {          // results_write.f:64-88
        fillArm(0, recon.e, main.FP.e, orig.e, main.target.z * s.e_sin_th);
        fillArm(12, recon.p, main.FP.p, orig.p, -main.target.z * s.p_sin_th);
    } else if (side == ArmSide::ElectronLeft) {    // results_write.f:89-117
        fillArm(0, recon.p, main.FP.p, orig.p, main.target.z * s.p_sin_th);
        fillArm(12, recon.e, main.FP.e, orig.e, -main.target.z * s.e_sin_th);
    }
    ntu[24] = recon.q / units::kMeVPerGeV;
    ntu[25] = recon.nu / units::kMeVPerGeV;
    ntu[26] = recon.Q2 / 1.e6;
    ntu[27] = recon.W / units::kMeVPerGeV;
    ntu[28] = recon.epsilon;
    ntu[29] = recon.Em / units::kMeVPerGeV;
    ntu[30] = recon.Pm / units::kMeVPerGeV;
    ntu[31] = recon.theta_pq;
    ntu[32] = recon.phi_pq;
    ntu[33] = corrsing / units::kMeVPerGeV;
    ntu[34] = Pm_Heepx / units::kMeVPerGeV;
    ntu[35] = Pm_Heepy / units::kMeVPerGeV;
    ntu[36] = Pm_Heepz / units::kMeVPerGeV;
    ntu[37] = recon.PmPar / units::kMeVPerGeV;
    ntu[38] = recon.PmPer / units::kMeVPerGeV;
    ntu[39] = recon.PmOop / units::kMeVPerGeV;
    ntu[40] = -main.target.rastery;
    ntu[41] = main.ntup.radphot / units::kMeVPerGeV;
    ntu[42] = main.sigcc;
    ntu[43] = main.weight;
    ntu[44] = recon.e.theta;
    ntu[45] = recon.p.theta;
    return side != ArmSide::Unknown;
}

}  // namespace simc
