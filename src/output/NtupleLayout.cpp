// src/output/NtupleLayout.cpp — see include/simc/output/NtupleLayout.h.
// The ntuple layouts of the Fortran (column names, sources, units, origin) and the per-event fill: (e,e'p) 46 columns,
// pion/kaon/delta 53 (+2 kaon), SIDIS 56, rho 59, each optionally followed by the 8 polarized-target angles and the 12
// pi0-photon columns. Column order and values are those of NtupleInit.f / results_write.f, so the h10 trees of the
// two codes compare tag by tag. Used by RootTreeSink (the tree), the summary (column list) and the documentation.
// Fortran: NtupleInit.f, results_write.f:results_ntu_write (l.1-271).
#include "simc/output/NtupleLayout.h"

#include <mutex>

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

// results_write.f:64/89 through the spectrometer table: beam-right arms (HMS, HRS-R, calo HMS side) put the
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

// NtupleInit.f:160-166 (doing_kaon, no target field): the 53 pion/kaon/delta columns + saghai, factor
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

namespace {
// results_write.f:155-162 / 213-220: the polarized-target angles (using_tgt_field)
void FillTargetFieldColumns(const MainEvent& main, const Event& vertex, const Event& recon, std::vector<double>& ntu,
                            std::size_t at) {
    ntu[at + 0] = recon.theta_tarq;
    ntu[at + 1] = recon.phi_targ;
    ntu[at + 2] = recon.beta;
    ntu[at + 3] = recon.phi_s;
    ntu[at + 4] = recon.phi_c;
    ntu[at + 5] = main.beta;
    ntu[at + 6] = vertex.phi_s;
    ntu[at + 7] = vertex.phi_c;
}
// results_write.f:172-185 / 230-243: the pi0 photons (12 columns from `at`)
void FillPhotonColumns(const MainEvent& main, std::vector<double>& ntu, std::size_t at) {
    ntu[at + 0] = main.ntup.xcal_gamma1;
    ntu[at + 1] = main.ntup.ycal_gamma1;
    for (std::size_t i = 0; i < 4; ++i) ntu[at + 2 + i] = main.ntup.gamma1[i];
    ntu[at + 6] = main.ntup.xcal_gamma2;
    ntu[at + 7] = main.ntup.ycal_gamma2;
    for (std::size_t i = 0; i < 4; ++i) ntu[at + 8 + i] = main.ntup.gamma2[i];
}
}  // namespace

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
    // pfermi column = |pfermi| with the sign of pfermi.q_hat (|d|/d is that sign; d = 0 -> +)
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
    std::size_t next = 53;   // first column after phipqi
    if (s.flags.using_tgt_field) {   // l.154-162
        FillTargetFieldColumns(main, vertex, recon, ntu, 53);
        next = 61;
    }
    if (s.flags.doing_kaon && ntu.size() >= next + 2) {   // l.163-166 / 168-171
        ntu[next] = main.ntup.sigcm1;       // sigcm - saghai model
        ntu[next + 1] = main.ntup.sigcm2;   // sigcm - factorized.
        next += 2;
    }
    if (s.flags.doing_pizero && ntu.size() >= next + 12) {   // l.172-185
        // FORTRAN-BUG(D-067 T14): original Fortran results_write.f:172-185 writes the photons to ntu(54..65) also with
        // using_tgt_field, overwriting the 8 angle columns (tags 54-61), and leaves the photon tags 66-73 at 0;
        // corrected: the photons go to their own tags (62-73). fortran_compat keeps the Fortran columns.
        FillPhotonColumns(main, ntu, s.fortran_compat ? 53 : next);
    }
    return side != ArmSide::Unknown;
}

// NtupleInit.f:192-237 (doing_semi, no target field, no pi0)
const std::vector<NtupleColumn>& SemiNtupleColumns() {
    static const std::vector<NtupleColumn> semi = [] {
        std::vector<NtupleColumn> c(EepNtupleColumns().begin(), EepNtupleColumns().begin() + 33);
        const std::vector<NtupleColumn> tail = {
            {"missmass", "main.ntup.mm/1000 (complete_recon_ev, nucleon)", "GeV", "results_write.f:188"},
            {"ppi", "recon.p.P/1000", "GeV/c", "results_write.f:189"},
            {"t", "main.ntup.t/1e6 (recon)", "GeV^2", "results_write.f:190"},
            {"fry", "-main.target.rastery", "cm", "results_write.f:191"},
            {"radphot", "main.ntup.radphot/1000", "GeV", "results_write.f:192"},
            {"siglab", "main.sigcc (peepiX)", "ub/MeV^2/sr^2", "results_write.f:193"},
            {"sigcent", "main.sigcent (never assigned: 0, D-051 F5)", "", "results_write.f:194"},
            {"Weight", "main.weight", "", "results_write.f:195"},
            {"decdist", "main.ntup.decdist (decay-OFF: survival probability)", "cm or 1", "results_write.f:196"},
            {"Mhadron", "sqrt(main.ntup.Mh2_final)", "MeV", "results_write.f:197"},
            {"z", "recon.zhad", "", "results_write.f:198"},
            {"zi", "vertex.zhad", "", "results_write.f:199"},
            {"pt2", "recon.pt2/1e6", "GeV^2", "results_write.f:200"},
            {"pt2i", "vertex.pt2/1e6", "GeV^2", "results_write.f:201"},
            {"xbj", "recon.xbj", "", "results_write.f:202"},
            {"xbji", "vertex.xbj", "", "results_write.f:203"},
            {"thqi", "acos(vertex.uq.z)", "rad", "results_write.f:204"},
            {"sighad", "main.ntup.sigcm (peepiX sighad)", "GeV^-2", "results_write.f:205"},
            {"jacobian", "main.davejac (peepiX jacobian*1000)", "", "results_write.f:206"},
            {"centjac", "main.johnjac (never assigned: 0, D-051 F5)", "", "results_write.f:207"},
            {"pfermi", "pfer/1000*sign(pfermi_hat.uq) (vertex uq; 0/0 = NaN for pfer = 0 in compat, D-051 F3)", "GeV/c",
             "results_write.f:208-209"},
            {"xfermi", "main.ntup.xfermi (x_Bj with Fermi motion)", "", "results_write.f:210"},
            {"phipqi", "main.phi_pq", "rad", "results_write.f:211"},
        };
        c.insert(c.end(), tail.begin(), tail.end());
        return c;
    }();
    return semi;
}

// Fortran: results_write.f:187-215 (doing_semi; no target field, no pi0)
bool FillSemiNtuple(const RunSetup& s, const MainEvent& main, const Event& vertex, const Event& orig,
                    const Event& recon, std::vector<double>& ntu) {
    const ArmSide side = FillArms(s, main, orig, recon, ntu.data());
    const FermiMomentum& F = main.pfermi;
    ntu[33] = main.ntup.mm / 1000.;
    ntu[34] = recon.p.P / 1000.;
    ntu[35] = main.ntup.t / 1.e6;
    ntu[36] = -main.target.rastery;
    ntu[37] = main.ntup.radphot / 1000.;
    ntu[38] = main.sigcc;
    ntu[39] = main.sigcent;
    ntu[40] = main.weight;
    ntu[41] = main.ntup.decdist;
    ntu[42] = std::sqrt(main.ntup.Mh2_final);
    ntu[43] = recon.zhad;
    ntu[44] = vertex.zhad;
    ntu[45] = recon.pt2 / 1.e06;
    ntu[46] = vertex.pt2 / 1.e06;
    ntu[47] = recon.xbj;
    ntu[48] = vertex.xbj;
    ntu[49] = std::acos(vertex.uq.z);
    ntu[50] = main.ntup.sigcm;
    ntu[51] = main.davejac;
    ntu[52] = main.johnjac;
    double dummy = F.pferx * vertex.uq.x + F.pfery * vertex.uq.y + F.pferz * vertex.uq.z;
    // FORTRAN-BUG(D-051 F3): original Fortran results_write.f:208-209 has no `if (dummy.eq.0) dummy=1.e-20` (the pion
    // block, l.143, has it): pfer = 0 (every H and A > 2 event) gives 0/0 = NaN in the ntuple column pfermi; corrected
    // as the pion block (pfermi = 0).
    if (!s.fortran_compat && dummy == 0) dummy = 1.e-20;
    ntu[53] = F.pfer / 1000. * std::fabs(dummy) / dummy;
    ntu[54] = main.ntup.xfermi;
    ntu[55] = main.phi_pq;
    std::size_t next = 56;
    if (s.flags.using_tgt_field) {   // l.212-220
        FillTargetFieldColumns(main, vertex, recon, ntu, 56);
        next = 64;
    }
    if (s.flags.doing_rho) {   // l.221-229
        if (s.flags.using_tgt_field) {
            ntu[64] = main.ntup.rhomass;
            ntu[65] = main.ntup.rhotheta;
            // FORTRAN-BUG(D-067 T14): original Fortran results_write.f:221-224 never writes the tag mmnuc (67) with
            // using_tgt_field (always 0); corrected: ntup%mmA/1000 as without the field. fortran_compat keeps the 0.
            if (!s.fortran_compat) ntu[66] = main.ntup.mmA / 1000.;
            next = 67;
        } else {
            ntu[56] = main.ntup.rhomass;
            ntu[57] = main.ntup.rhotheta;
            ntu[58] = main.ntup.mmA / 1000.;   // missmass (nucleus)HSV
            next = 59;
        }
    }
    if (s.flags.doing_pizero && ntu.size() >= next + 12) {   // l.230-243
        // FORTRAN-BUG(D-067 T14): with using_tgt_field the Fortran writes the photons to ntu(57..68), over the angle
        // columns; corrected: their own tags. fortran_compat keeps the Fortran columns.
        FillPhotonColumns(main, ntu, s.fortran_compat ? 56 : next);
    }
    return side != ArmSide::Unknown;
}

// NtupleInit.f (doing_semi with doing_pizero, no target field): the 56 SIDIS columns + the 12 photon columns
const std::vector<NtupleColumn>& SemiPizeroNtupleColumns() {
    static const std::vector<NtupleColumn> c = [] {
        std::vector<NtupleColumn> v = SemiNtupleColumns();
        const std::vector<NtupleColumn>& pi0 = PionNtupleColumns(true);
        v.insert(v.end(), pi0.end() - 12, pi0.end());   // same tags and sources (the semi lines are results_write.f:230-243)
        return v;
    }();
    return c;
}

// NtupleInit.f:192-264 (doing_rho, no target field): the SIDIS columns (some unused for rho) + Mrho, Thrho, mmnuc
const std::vector<NtupleColumn>& RhoNtupleColumns() {
    static const std::vector<NtupleColumn> rho = [] {
        std::vector<NtupleColumn> c = SemiNtupleColumns();
        c[35].source = "main.ntup.t/1e6 (recon, the detected decay pion: Mh2 = Mpi2)";
        c[38].source = "main.sigcc (peerho of the rho vertex)";
        c[39].source = "main.sigcent (never assigned: 0)";
        c[44].source = "vertex.zhad (never assigned for rho: 0)";
        c[46].source = "vertex.pt2/1e6 (never assigned for rho: 0)";
        c[50].source = "main.ntup.sigcm (peerho dsigma/dOmega_cm)";
        c[51].source = "main.davejac (peerho, not applied)";
        c[52].source = "main.johnjac (peerho)";
        c[54].source = "main.ntup.xfermi (never assigned for rho: 0)";
        c[55].source = "main.phi_pq (of the rho)";
        const std::vector<NtupleColumn> tail = {
            {"Mrho", "main.ntup.rhomass (Breit-Wigner rho mass, generate_rho.f:79)", "MeV", "results_write.f:227"},
            {"Thrho", "main.ntup.rhotheta (decay theta in the rho rest frame w.r.t. the lab z axis, rho_decay.f:66)",
             "rad", "results_write.f:228"},
            {"mmnuc", "main.ntup.mmA/1000 (missing mass, nucleus)", "GeV", "results_write.f:229"},
        };
        c.insert(c.end(), tail.begin(), tail.end());
        return c;
    }();
    return rho;
}

// Fortran: results_write.f:187-229 (doing_rho; no target field)
bool FillRhoNtuple(const RunSetup& s, const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                   std::vector<double>& ntu) {
    return FillSemiNtuple(s, main, vertex, orig, recon, ntu);   // (the rho columns are written there, l.221-229)
}

// NtupleInit.f:142-159 / 239-256 (using_tgt_field): the 8 polarized-target angle columns
const std::vector<NtupleColumn>& TargetFieldNtupleColumns() {
    static const std::vector<NtupleColumn> c = {
        {"th_tarq", "recon.theta_tarq (angle q - polarization; NaN for targ_pol = 0, D-067 T16)", "rad",
         "results_write.f:155/213"},
        {"phitarq", "recon.phi_targ", "rad", "results_write.f:156/214"},
        {"beta", "recon.beta", "rad", "results_write.f:157/215"},
        {"phis", "recon.phi_s (Sivers angle)", "rad", "results_write.f:158/216"},
        {"phic", "recon.phi_c (Collins angle)", "rad", "results_write.f:159/217"},
        {"betai", "main.beta", "rad", "results_write.f:160/218"},
        {"phisi", "vertex.phi_s", "rad", "results_write.f:161/219"},
        {"phici", "vertex.phi_c", "rad", "results_write.f:162/220"},
    };
    return c;
}

namespace {
// NtupleInit.f:101-191 with using_tgt_field: pion 53 + 8 angles [+ saghai, factor (kaon)] [+ 12 photons (pi0)]
const std::vector<NtupleColumn>& PionFieldNtupleColumns(bool kaon, bool pizero) {
    static std::vector<NtupleColumn> cache[2][2];
    static std::mutex mutex;   // the sinks of the parallel loop's worker threads may ask at the same time
    std::lock_guard<std::mutex> lock(mutex);
    std::vector<NtupleColumn>& c = cache[kaon][pizero];
    if (c.empty()) {
        c = PionNtupleColumns(false);
        const auto& f = TargetFieldNtupleColumns();
        c.insert(c.end(), f.begin(), f.end());
        if (kaon) {
            const auto& k = KaonNtupleColumns();
            c.insert(c.end(), k.end() - 2, k.end());
            c[c.size() - 2].fortran = "results_write.f:163";
            c[c.size() - 1].fortran = "results_write.f:164";
        }
        if (pizero) {
            const auto& g = PionNtupleColumns(true);
            c.insert(c.end(), g.end() - 12, g.end());
        }
    }
    return c;
}
// NtupleInit.f:192-264 with using_tgt_field: semi 56 + 8 angles [+ Mrho, Thrho, mmnuc (rho)] [+ 12 photons (pi0)]
const std::vector<NtupleColumn>& SemiFieldNtupleColumns(bool rho, bool pizero) {
    static std::vector<NtupleColumn> cache[2][2];
    static std::mutex mutex;   // the sinks of the parallel loop's worker threads may ask at the same time
    std::lock_guard<std::mutex> lock(mutex);
    std::vector<NtupleColumn>& c = cache[rho][pizero];
    if (c.empty()) {
        c = rho ? std::vector<NtupleColumn>(RhoNtupleColumns().begin(), RhoNtupleColumns().begin() + kSemiNtupleSize)
                : SemiNtupleColumns();
        const auto& f = TargetFieldNtupleColumns();
        c.insert(c.end(), f.begin(), f.end());
        if (rho) {
            const auto& r = RhoNtupleColumns();
            c.insert(c.end(), r.end() - 3, r.end());
        }
        if (pizero) {
            const auto& g = PionNtupleColumns(true);
            c.insert(c.end(), g.end() - 12, g.end());
        }
    }
    return c;
}
}  // namespace

const std::vector<NtupleColumn>& NtupleColumnsFor(const RunSetup& s) {
    const bool field = s.flags.using_tgt_field || s.params.using_tgt_field;
    // NtupleInit.f:101: doing_pion .or. doing_kaon .or. doing_delta share the 53 columns.
    // (using_tgt_field is checked on the flags and on the raw parameter: either one selects the field layouts.)
    if (field && (s.flags.doing_pion || s.flags.doing_kaon || s.flags.doing_delta))
        return PionFieldNtupleColumns(s.flags.doing_kaon, s.flags.doing_pizero);
    if (field && (s.flags.doing_semi || (s.flags.doing_rho && !s.flags.doing_pizero)))
        return SemiFieldNtupleColumns(!s.flags.doing_semi, s.flags.doing_pizero);
    if ((s.flags.doing_pion || s.flags.doing_delta) && !field) return PionNtupleColumns(s.flags.doing_pizero);
    if (s.flags.doing_kaon && !field) return KaonNtupleColumns();
    if (s.flags.doing_semi && !field) return s.flags.doing_pizero ? SemiPizeroNtupleColumns() : SemiNtupleColumns();
    if (s.flags.doing_rho && !s.flags.doing_pizero && !field) return RhoNtupleColumns();
    if (UsesEepLayout(s.flags)) return EepNtupleColumns();
    throw std::runtime_error("ntuple: no layout for this reaction (rho with pi0 is not set up)");
}

bool FillNtuple(const RunSetup& s, const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                std::vector<double>& ntu) {
    if (s.flags.doing_pion || s.flags.doing_kaon || s.flags.doing_delta)
        return FillPionNtuple(s, main, vertex, orig, recon, ntu);
    if (s.flags.doing_semi) return FillSemiNtuple(s, main, vertex, orig, recon, ntu);
    if (s.flags.doing_rho) return FillRhoNtuple(s, main, vertex, orig, recon, ntu);
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
    // sin(..)**2 is x*x in gfortran; std::pow(x,2) is NOT folded to x*x by gcc 16 under -fmath-errno (not bitwise)
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
