// include/simc/reactions/PionProduction.h
//
// IReaction "pion_production" (task P9.2): exclusive pion electroproduction, the doing_pion branches of event.f.
// One class for every final state the Fortran supports with the HMS/SHMS (which_pion, doing_pizero, the target):
//   which_pion 0 / 1     pi+ (pi0 from p) / pi- (pi0 from n), quasifree: H (doing_hydpi), D (doing_deutpi, Fermi
//                        momentum from deut.dat), A > 2 (doing_hepi, Fermi momentum + Benhar missing energy)
//   which_pion 2 / 3     pi Delta final states (relativistic Breit-Wigner Delta mass per complete_ev call)
//   which_pion 10 / 11   coherent production on A > 2 (treated as hydrogen with Mtar_struck = M_A)
// The branch is chosen by the RunSetup flags exactly as in the Fortran.
//
//   generate()        event.f:generate l.214-428 after the beam part: electron and hadron angles, electron energy in
//                     the sumEgen window (gen_weight), physics_angles, COMMON /pfermi/ (MainEvent::pfermi) and Em
//                     (deut.dat / he3.dat momentum distributions, generate_em), complete_ev, generate_rad or orig = vertex.
//   completeEv()      event.f:complete_ev l.432-1054, pion branch: Delta mass (l.478-496, 611-631), hadron energy from
//                     the quadratic incl. Fermi motion (l.633-698), main%W/epsilon/theta_pq/t/tmin/q2/phi_pq
//                     (l.720-786), pizero_decay (l.899-901), Pm components, Trec, krel (l.914-977), both 1/cos**3
//                     jacobians (l.1013-1024), trip_thru_target 2/3, radc_init_ev.
//   completeRecon()   EepReactionBase (complete_recon_ev is reaction-generic; it also returns ntup%mm/mmA/t).
//   completeMain()    event.f:complete_main l.1363-1569, pion branch: SF_weight = 1, sigcc = peepi(vertex, main) (the
//                     run's IPhysicsModel, "pion", through sigccMain: it also writes main%thetacm ... and ntup%sigcm),
//                     the Delta factors (l.1464-1494), sigcc_recon = 1, tgtweight = targ%N (which_pion 1/11) or targ%Z,
//                     the Coulomb focusing factor and the weight. Also copies ntup%mm/mmA/t of complete_recon_ev into
//                     main.ntup (the Fortran COMMON /ntuple_variables/).
//   prepareCentral()  calculate_central (simc.f:1236-1290), see D-038 F5 below.
//   genVolume()       EepReactionBase (simc.f:383-396).
// Decay in flight of the pion (doing_decay, ctau) happens in the hadron spectrometer (ISpectrometer, DecayState of the
// optics engine; D-037 corrections there); the reaction has no part in it.
//
// Random numbers (as the Fortran, D-024 order): generate 4 grnd (angles) + 1 (e energy) [+ 4 D/A: ranprob, pfer,
// cos theta, phi; + 2 generate_em for A > 2]; complete_ev: 1 grnd for the Delta mass (which_pion 2/3), 2 for the pi0
// decay, 8 gauss1 (trip_thru_target 2/3).
//
// JSON config (all optional): {"data_dir": "<dir with delta_relativistic_bw.inp>" (default data/pion), "verbose": bool}.
// The momentum distributions and the Benhar table come from RunSetup::data_dir (data/spectral).
//
// fortran_compat (D-031, docs/FORTRAN_DIFFERENCES.md):
//   FORTRAN-BUG(D-038 F5)  simc.f:1236-1290: the central cross section (.hist CENTRAL.sigcc) is computed with
//                          main0%epsilon = t = phi_pq = 0 (only complete_recon_ev(vertex0) is called). Corrected: main0
//                          gets the central W, epsilon, theta_pq, t, tmin, q2, phi_pq (complete_ev formulas on vertex0).
//   FORTRAN-DIFF-CANDIDATE(D-038 F9) complete_ev overwrites targ%Mrec_struck with a new Breit-Wigner Delta mass on every
//                          call, and radiation tail 1 calls complete_ev again (a second Delta mass and pi0 decay): kept.
//   FORTRAN-DIFF-CANDIDATE(D-038 F4) coherent production with the nucleon model (in the "pion" model): kept.
//   FORTRAN-QUIRK(D-038 F11) ntup%krel = the SAVEd local krel, never set by the H/coherent branch (0).
#ifndef SIMC_REACTIONS_PIONPRODUCTION_H
#define SIMC_REACTIONS_PIONPRODUCTION_H

#include <memory>
#include <string>
#include <vector>

#include "simc/reactions/EepReactionBase.h"

namespace simc {
namespace nuclear {
class SfTable;
}
namespace reactions {

/// Fortran: event.f:complete_ev l.478-496, the Delta(1232) mass table (SAVEd locals mass_save, probsum of the first
/// call). Loaded from delta_relativistic_bw.inp.
struct DeltaMassTable {
    std::vector<double> mass;      ///< mass_save(i) [GeV]
    std::vector<double> probsum;   ///< normalised cumulative probability
    bool load(const std::string& path, std::vector<std::string>* msgs, std::string* err);
    /// event.f:621-631: the mass [MeV] whose probsum is closest to rn (first minimum of |rn - probsum|).
    double sample(double rn) const;
};

class PionProduction : public EepReactionBase {
public:
    explicit PionProduction(std::string name = "pion_production");
    ~PionProduction() override;

    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) override;
    bool completeEv(MainEvent& main, Event& vertex, IRandom& rng) override;
    bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                      const Event& recon) override;
    void prepareCentral(MainEvent& main0, const Event& vertex0) override;
    double spectralWeight(const Event&) override { return 1.0; }   // event.f:1397-1398
    /// simc.f:844 prints targ%Mrec_struck at the end of the run: the last Delta mass for which_pion 2/3.
    void finalTargetState(TargetInfo& targ) const override;

    /// COMMON targ%Mrec_struck as complete_ev leaves it (Breit-Wigner Delta mass for which_pion 2/3).
    double mrecStruck() const { return mrecStruck_; }
    /// Test hook (L1 tables): set the COMMON targ%Mrec_struck (e.g. the value of an earlier call).
    void setMrecStruck(double m) { mrecStruck_ = m; }
    const DeltaMassTable& deltaTable() const { return delta_; }
    nuclear::SfTable* sfTable() { return sf_.get(); }

private:
    /// event.f:720-786: main%W, epsilon, theta_pq, t, tmin, q2, phi_pq from the completed vertex.
    void pionMainKinematics(MainEvent& main, const Event& vertex) const;

    std::string dataDir_;
    bool verbose_{true};
    double mrecStruck_{0};             ///< COMMON targ%Mrec_struck (overwritten per call for Delta final states)
    DeltaMassTable delta_;
    std::unique_ptr<nuclear::SfTable> sf_;   ///< doing_hepi: sf_lookup_init table for generate_em
};

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_PIONPRODUCTION_H
