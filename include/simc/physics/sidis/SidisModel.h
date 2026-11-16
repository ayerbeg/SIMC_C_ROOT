/**
 * @file
 * @ingroup physics
 * @brief simc::sidis::SidisModel: IPhysicsModel "sidis", the semi-inclusive (e,e'h)X cross section (semi_physics.f:peepiX).
 *
 * IPhysicsModel plug-in "sidis" (task P12.1): the semi-inclusive (e,e'h)X cross section of the Fortran SIMC,
 * semi_physics.f:peepiX, called by event.f:complete_main (l.1520-1522) for doing_semi:
 *     main%sigcc       = peepiX(vertex, vertex0, main, survivalprob, .FALSE.)   (sigccMain / sigcc)
 *     main%sigcc_recon = 1.0                                                   (sigccRecon)
 * peepiX also writes main%davejac, ntup%sigcm (= sighad), ntup%xfermi (do_fermi) and, when the decay is not simulated
 * (doing_decay = 0), the survival probability of the hadron into the COMMON decdist (main.ntup.decdist), from where
 * the reaction applies it to the weight for doing_semika (event.f:1560-1561; P12.2).
 * Data: the CTEQ5M table cteq5m.tbl (SetCtq5(1)) and, for kaons, the fDSS NLO grid KANLO.GRID, read at init() from the
 * JSON key "data_dir" (default data/sidis of the source tree, SIMC_SIDIS_DATA_DIR; byte copies of the Fortran files).
 *
 * JSON config (all optional):
 *     { "data_dir": "<dir with cteq5m.tbl and KANLO.GRID>", "verbose": true }
 * Run state read from the set-up: targ%A, targ%Z, targ%Mtar_struck, doing_semipi/semika/pizero/hplus, do_fermi,
 * doing_decay, hadron_arm, ctau, spec%p%P, Mh2, sigc_* (dead doing_cent branch only).
 * fortran_compat (D-031) restores the Fortran as built; the default corrects D-051 F1 (0/0 when all PDFs are 0),
 * F2 (negative Fermi x_Bj: weight 0 instead of the STOP), F4 (pi0 dbar FF and mass) and F6 (survival probability:
 * hadron momentum, zaero for the SHMS, sqrt in the path). Unsupported
 * combinations are refused in init() (sidis::checkSupported). docs/translation/markdown/FORTRAN_DIFFERENCES.md.
 */
#ifndef SIMC_PHYSICS_SIDIS_SIDISMODEL_H
#define SIMC_PHYSICS_SIDIS_SIDISMODEL_H

#include <string>
#include <vector>

#include "simc/physics/sidis/Cteq5Pdf.h"
#include "simc/physics/sidis/Fdss.h"
#include "simc/physics/sidis/SidisCrossSection.h"
#include "simc/stages/IPhysicsModel.h"

namespace simc {

struct RunSetup;

namespace sidis {

/// @addtogroup physics
/// @{

/// The SIDIS input combinations SIMC++ refuses (each entry is an error message); empty = supported.
/// Checked: doing_semi set, exactly one of doing_semipi/doing_semika (dbase.f:125-132), no doing_pizero (pi0 SIDIS:
/// pi0 SIDIS is supported since P14.3, with the calorimeter arm), no doing_phsp (sigcc = 1 there, event.f:1446).
/// @param s the run set-up.
/// @return the error messages.
std::vector<std::string> checkSupported(const RunSetup& s);
/// The warnings for SIDIS inputs that the Fortran reads and ignores (pt_b_param, sigc_*).
/// @param s the run set-up.
/// @return the warnings.
std::vector<std::string> ignoredInputNotes(const RunSetup& s);
/// PeepiXConfig from the set-up (fortran_compat selects the D-031 switches).
/// @param s the run set-up.
/// @return the peepiX constants.
PeepiXConfig makeConfig(const RunSetup& s);

/// @brief IPhysicsModel "sidis": peepiX with CTEQ5 and fDSS. Owns the tables; not thread-safe.
class SidisModel : public IPhysicsModel {
public:
    std::string name() const override { return "sidis"; }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Checks the inputs (checkSupported: std::invalid_argument / simc::NotPorted), reads cteq5m.tbl and, for kaons,
    /// KANLO.GRID (std::invalid_argument if they cannot be read).
    /// @param ctx the stage context.
    void init(const StageContext& ctx) override;

    /// peepiX(vertex, vertex, main) on copies (main, the SAVEd state and the SAVEd PDF warnings are not changed).
    /// @param vertex the vertex.
    /// @param main the main event record (copied).
    /// @return the cross section.
    double sigcc(const Event& vertex, const MainEvent& main) override;
    /// main%sigcc = peepiX(vertex, vertex0 (unused for doing_cent = .FALSE.), main, survivalprob, .FALSE.): also writes
    /// main%davejac, ntup%sigcm, ntup%xfermi and (doing_decay = 0) ntup%decdist; updates the SAVEd locals.
    /// @param vertex the vertex.
    /// @param main the main event record (written).
    /// @return the cross section.
    double sigccMain(const Event& vertex, MainEvent& main) override;
    /// event.f:1522 main%sigcc_recon = 1.0
    /// @return 1.
    double sigccRecon(const Event&, const MainEvent&) override { return 1.0; }

    /// @return the peepiX constants.
    const PeepiXConfig& config() const { return cfg_; }
    /// @return the peepiX constants (modifiable).
    PeepiXConfig& config() { return cfg_; }   ///< test hook
    /// @return the SAVEd state.
    PeepiXState& state() { return st_; }
    /// @return the CTEQ5 set.
    Cteq5Pdf& pdf() { return pdf_; }
    /// @return the fDSS grid.
    const FdssGrid& fdss() const { return fdss_; }
    /// The SAVEd local survivalprob of complete_main (event.f:1372), updated by peepiX when doing_decay = 0.
    /// @return the survival probability.
    double survivalProbability() const { return survivalprob_; }
    /// @param p the survival probability.
    void setSurvivalProbability(double p) { survivalprob_ = p; }

private:
    std::string data_dir_;
    bool verbose_{true};
    PeepiXConfig cfg_;
    PeepiXState st_;
    Cteq5Pdf pdf_;
    FdssGrid fdss_;
    double survivalprob_{0};
};

/// @}

}  // namespace sidis
}  // namespace simc

#endif  // SIMC_PHYSICS_SIDIS_SIDISMODEL_H
