// src/output/NtupleLayout.cpp
// Fortran: NtupleInit.f, results_write.f:results_ntu_write (l.1-271).
#include "simc/output/NtupleLayout.h"

#include "simc/core/SpectrometerId.h"
#include "simc/core/Units.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

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

// results_write.f:64/89 through the spectrometer table (D-032): beam-right arms (HMS, HRS-R, calo HMS side) put the
// electron in the "hs" columns. results_write.f:90 also lists the retired code 6, which the input rejects.
ArmSide ElectronArmSide(int electron_arm) {
    const std::optional<SpectrometerId> id = SpectrometerFromCode(electron_arm);
    if (!id) return ArmSide::Unknown;
    return IsBeamRight(*id) ? ArmSide::ElectronRight : ArmSide::ElectronLeft;
}

bool UsesEepLayout(const ReactionFlags& f) {
    // NtupleInit.f: the pion/kaon/delta branch comes first, then semi/rho, then (e,e'p).
    if (f.doing_pion || f.doing_kaon || f.doing_delta) return false;
    if (f.doing_semi || f.doing_rho) return false;
    return f.doing_hyd_elast || f.doing_deuterium || f.doing_heavy;
}

namespace {
/// results_write.f:64-117: columns 1-24 (right-side "hs", left-side "ss" arm). Unknown arm: untouched.
ArmSide FillArms(const RunSetup& s, const MainEvent& main, const Event& orig, const Event& recon, double* ntu) {
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
    // results_write.f:121-129
    ntu[24] = recon.q / units::kMeVPerGeV;
    ntu[25] = recon.nu / units::kMeVPerGeV;
    ntu[26] = recon.Q2 / 1.e6;
    ntu[27] = recon.W / units::kMeVPerGeV;
    ntu[28] = recon.epsilon;
    ntu[29] = recon.Em / units::kMeVPerGeV;
    ntu[30] = recon.Pm / units::kMeVPerGeV;
    ntu[31] = recon.theta_pq;
    ntu[32] = recon.phi_pq;
    return side;
}
}  // namespace

const std::vector<NtupleColumn>& PionNtupleColumns(bool pizero) {
    static const std::vector<NtupleColumn> pion = [] {
        std::vector<NtupleColumn> c(EepNtupleColumns().begin(), EepNtupleColumns().begin() + 33);
        const std::vector<NtupleColumn> tail = {
            {"missmass", "main.ntup.mm/1000 (complete_recon_ev: sqrt|Em^2-Pm^2| with sign)", "GeV", "results_write.f:134"},
            {"mmnuc", "main.ntup.mmA/1000 (missing mass on the nucleus)", "GeV", "results_write.f:135"},
            {"phad", "recon.p.P/1000", "GeV/c", "results_write.f:136"},
            {"t", "main.ntup.t/1e6 (recon)", "GeV^2", "results_write.f:137"},
            {"pmpar", "recon.PmPar/1000", "GeV/c", "results_write.f:138"},
            {"pmper", "recon.PmPer/1000", "GeV/c", "results_write.f:139"},
            {"pmoop", "recon.PmOop/1000", "GeV/c", "results_write.f:140"},
            {"fry", "-main.target.rastery", "cm", "results_write.f:141"},
            {"radphot", "main.ntup.radphot/1000", "GeV", "results_write.f:142"},
            {"pfermi", "pfer/1000*sign(pfermi_hat.uq) (vertex uq)", "GeV/c", "results_write.f:143-145"},
            {"siglab", "main.sigcc", "ub/MeV/sr^2", "results_write.f:146"},
            {"sigcm", "main.ntup.sigcm", "ub/MeV^2/rad", "results_write.f:147"},
            {"Weight", "main.weight", "", "results_write.f:148"},
            {"decdist", "main.ntup.decdist", "cm", "results_write.f:149"},
            {"Mhadron", "sqrt(main.ntup.Mh2_final)", "MeV", "results_write.f:150"},
            {"pdotqhat", "pfer/1000*(pfermi_hat.uq)", "GeV/c", "results_write.f:151"},
            {"Q2i", "vertex.Q2/1e6", "GeV^2", "results_write.f:152"},
            {"Wi", "main.w/1000", "GeV", "results_write.f:153"},
            {"ti", "main.t/1e6", "GeV^2", "results_write.f:154"},
            {"phipqi", "main.phi_pq", "rad", "results_write.f:155"},
        };
        c.insert(c.end(), tail.begin(), tail.end());
        return c;
    }();
    static const std::vector<NtupleColumn> pi0 = [] {
        std::vector<NtupleColumn> c = pion;
        const std::vector<NtupleColumn> tail = {
            {"xcal_gamma1", "main.ntup.xcal_gamma1 (calorimeter only; 0 with HMS/SHMS)", "cm", "results_write.f:174"},
            {"ycal_gamma1", "main.ntup.ycal_gamma1", "cm", "results_write.f:175"},
            {"Egamma1", "main.ntup.gamma1[0]", "MeV", "results_write.f:176"},
            {"Pgamma1x", "main.ntup.gamma1[1]", "MeV/c", "results_write.f:177"},
            {"Pgamma1y", "main.ntup.gamma1[2]", "MeV/c", "results_write.f:178"},
            {"Pgamma1z", "main.ntup.gamma1[3]", "MeV/c", "results_write.f:179"},
            {"xcal_gamma2", "main.ntup.xcal_gamma2", "cm", "results_write.f:180"},
            {"ycal_gamma2", "main.ntup.ycal_gamma2", "cm", "results_write.f:181"},
            {"Egamma2", "main.ntup.gamma2[0]", "MeV", "results_write.f:182"},
            {"Pgamma2x", "main.ntup.gamma2[1]", "MeV/c", "results_write.f:183"},
            {"Pgamma2y", "main.ntup.gamma2[2]", "MeV/c", "results_write.f:184"},
            {"Pgamma2z", "main.ntup.gamma2[3]", "MeV/c", "results_write.f:185"},
        };
        c.insert(c.end(), tail.begin(), tail.end());
        return c;
    }();
    return pizero ? pi0 : pion;
}

// NtupleInit.f:160-166 (doing_kaon, no target field): the 53 pion/kaon/delta columns + saghai, factor (task P10.2)
const std::vector<NtupleColumn>& KaonNtupleColumns() {
    static const std::vector<NtupleColumn> kaon = [] {
        std::vector<NtupleColumn> c = PionNtupleColumns(false);
        c.push_back({"saghai", "main.ntup.sigcm1 (Saghai model; 0 as built, D-043 F1)", "ub/MeV^2/rad", "results_write.f:164"});
        c.push_back({"factor", "main.ntup.sigcm2 (factorized model, = sigcm)", "ub/MeV^2/rad", "results_write.f:165"});
        return c;
    }();
    return kaon;
}

bool UsesPionLayout(const ReactionFlags& f) { return f.doing_pion || f.doing_kaon || f.doing_delta; }

// Fortran: results_write.f:131-186 (doing_pion; no target field)
bool FillPionNtuple(const RunSetup& s, const MainEvent& main, const Event& vertex, const Event& orig,
                    const Event& recon, std::vector<double>& ntu) {
    const ArmSide side = FillArms(s, main, orig, recon, ntu.data());
    const FermiMomentum& F = main.pfermi;
    ntu[33] = main.ntup.mm / 1000.;
    ntu[34] = main.ntup.mmA / 1000.;
    ntu[35] = recon.p.P / 1000.;
    ntu[36] = main.ntup.t / 1.e6;
    ntu[37] = recon.PmPar / 1000.;
    ntu[38] = recon.PmPer / 1000.;
    ntu[39] = recon.PmOop / 1000.;
    ntu[40] = -main.target.rastery;
    ntu[41] = main.ntup.radphot / 1000.;
    double dummy = F.pferx * vertex.uq.x + F.pfery * vertex.uq.y + F.pferz * vertex.uq.z;
    if (dummy == 0) dummy = 1.e-20;
    ntu[42] = F.pfer / 1000. * std::fabs(dummy) / dummy;
    ntu[43] = main.sigcc;
    ntu[44] = main.ntup.sigcm;
    ntu[45] = main.weight;
    ntu[46] = main.ntup.decdist;
    ntu[47] = std::sqrt(main.ntup.Mh2_final);
    ntu[48] = F.pfer / 1000. * dummy;
    ntu[49] = vertex.Q2 / 1.e6;
    ntu[50] = main.w / 1.e3;
    ntu[51] = main.t / 1.e6;
    ntu[52] = main.phi_pq;
    if (s.flags.doing_kaon && ntu.size() >= static_cast<std::size_t>(kKaonNtupleSize)) {   // l.163-166 (no target field)
        ntu[53] = main.ntup.sigcm1;   // sigcm - saghai model
        ntu[54] = main.ntup.sigcm2;   // sigcm - factorized.
    }
    if (s.flags.doing_pizero && ntu.size() >= static_cast<std::size_t>(kPionPizeroNtupleSize)) {   // l.172-185
        ntu[53] = main.ntup.xcal_gamma1;
        ntu[54] = main.ntup.ycal_gamma1;
        for (int i = 0; i < 4; ++i) ntu[static_cast<std::size_t>(55 + i)] = main.ntup.gamma1[static_cast<std::size_t>(i)];
        ntu[59] = main.ntup.xcal_gamma2;
        ntu[60] = main.ntup.ycal_gamma2;
        for (int i = 0; i < 4; ++i) ntu[static_cast<std::size_t>(61 + i)] = main.ntup.gamma2[static_cast<std::size_t>(i)];
    }
    return side != ArmSide::Unknown;
}

const std::vector<NtupleColumn>& NtupleColumnsFor(const RunSetup& s) {
    // NtupleInit.f:101: doing_pion .or. doing_kaon .or. doing_delta share the 53 columns (doing_delta: task P11.1)
    if ((s.flags.doing_pion || s.flags.doing_delta) && !s.flags.using_tgt_field && !s.params.using_tgt_field)
        return PionNtupleColumns(s.flags.doing_pizero);
    if (s.flags.doing_kaon && !s.flags.using_tgt_field && !s.params.using_tgt_field) return KaonNtupleColumns();
    if (UsesEepLayout(s.flags)) return EepNtupleColumns();
    throw std::runtime_error("ntuple: only the (e,e'p)/phase-space (NtupleInit.f:296-321), the exclusive pion and the "
                             "kaon/delta (NtupleInit.f:101-190, no target field) layouts are implemented so far");
}

bool FillNtuple(const RunSetup& s, const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                std::vector<double>& ntu) {
    if (s.flags.doing_pion || s.flags.doing_kaon || s.flags.doing_delta)
        return FillPionNtuple(s, main, vertex, orig, recon, ntu);
    std::array<double, kEepNtupleSize> a{};
    std::copy(ntu.begin(), ntu.begin() + std::min<std::size_t>(ntu.size(), kEepNtupleSize), a.begin());
    const bool ok = FillEepNtuple(s, main, vertex, orig, recon, a);
    ntu.assign(a.begin(), a.end());
    return ok;
}

bool FillEepNtuple(const RunSetup& s, const MainEvent& main, const Event& /*vertex*/, const Event& orig,
                   const Event& recon, std::array<double, kEepNtupleSize>& ntu) {
    const double Ebeam = s.params.Ebeam;
    // results_write.f:54-62
    // sin(..)**2 is x*x in gfortran; std::pow(x,2) is NOT folded to x*x by gcc 16 under -fmath-errno (P10.4)
    const double sinh2 = std::sin(recon.e.theta / 2.);
    const double poftheta = kMp * Ebeam / (2 * Ebeam * (sinh2 * sinh2) + kMp);
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
