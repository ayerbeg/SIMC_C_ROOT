/**
 * @file
 * @ingroup output
 * @brief Column layouts and fillers of the Fortran diagnostic ntuple (tree "h10"); ROOT-free.
 *
 * Column layout and filling of the Fortran diagnostic ntuple (tree "h10" after
 * util/root_tree/make_root_tree). ROOT-free so it can be unit-tested and reused.
 *
 * Fortran: NtupleInit.f (tags, l.36-321), results_write.f:results_ntu_write (l.1-271),
 *          make_root_tree.f + froot.c (one Float_t branch per tag, tree "h10").
 *
 * Implemented layouts: the (e,e'p) family, i.e. doing_hyd_elast / doing_deuterium /
 * doing_heavy (this includes doing_phsp, which keeps doing_hyd_elast & co.): 46 columns; exclusive pion production
 * (doing_pion, task P9.2): 53 columns (65 with doing_pizero); kaon (doing_kaon, P10.2): 55; Delta (P11.1) as pion; SIDIS
 * (doing_semi, P12.2): 56. rho comes with that reaction.
 */
#ifndef SIMC_OUTPUT_NTUPLELAYOUT_H
#define SIMC_OUTPUT_NTUPLELAYOUT_H

#include <array>
#include <string>
#include <vector>

#include "simc/core/Event.h"
#include "simc/core/RunSetup.h"

namespace simc {

/// @addtogroup output
/// @{

/// Description of one ntuple column (a static table entry; strings are literals).
struct NtupleColumn {
    const char* name;     ///< tag = ROOT branch name
    const char* source;   ///< C++ expression (for the right-side/left-side assignment see arm_note)
    const char* units;  ///< physical units of the column.
    const char* fortran;  ///< results_write.f line(s)
};

/// Number of columns of the (e,e'p) layout.
constexpr int kEepNtupleSize = 46;

/// The 46 (e,e'p) columns in Fortran order.
/// @return the static column table.
const std::vector<NtupleColumn>& EepNtupleColumns();

/// True when the electron arm is on the beam-right side (HMS, HRSR, calo on HMS side):
/// results_write.f:64 (SpectrometerId::IsBeamRight). Then the "hs*" columns carry the electron.
/// For the beam-left arms (SOS, HRSL, SHMS, calo on SHMS side; results_write.f:89) the "hs*" columns carry the hadron.
/// "hs" = right-side spectrometer, "ss" = left-side spectrometer.
/// Side of the electron arm (see ElectronArmSide()).
enum class ArmSide { ElectronRight /**< electron on the beam-right side */, ElectronLeft /**< electron on the beam-left side */, Unknown /**< unknown arm code */ };
/// @brief The side of the electron arm, which decides whether the "hs*" columns hold the electron or the hadron.
/// @param electron_arm arm code of the electron arm.
/// @return ElectronRight, ElectronLeft, or Unknown for an unknown code.
ArmSide ElectronArmSide(int electron_arm);

/// results_ntu_write for the (e,e'p) family. `ntu` keeps its previous content for
/// columns 1-24 when the electron arm is unknown, exactly like the Fortran static array
/// (results_write.f:118-119, -fno-automatic). Returns false in that case.
/// @param setup the run set-up.
/// @param main the main event record.
/// @param vertex the generated vertex.
/// @param orig the event after radiation.
/// @param recon the reconstructed event.
/// @param[in,out] ntu the ntuple row (static between calls, as in the Fortran).
/// @return false for an unknown electron arm.
bool FillEepNtuple(const RunSetup& setup, const MainEvent& main, const Event& vertex, const Event& orig,
                   const Event& recon, std::array<double, kEepNtupleSize>& ntu);

/// True if the reaction flags select the (e,e'p) ntuple layout (NtupleInit.f:296).
/// @param f the reaction flags.
/// @return true for the (e,e'p) layout.
bool UsesEepLayout(const ReactionFlags& f);

// ------------------------------------------------------------------------------------------- pion layout (P9.2)
constexpr int kPionNtupleSize = 53;        ///< NtupleInit.f:101-141 (doing_pion, no target field, no pi0)
constexpr int kPionPizeroNtupleSize = 65;  ///< + 12 pi0 photon columns (NtupleInit.f:165-190)
constexpr int kKaonNtupleSize = 55;        ///< doing_kaon: + saghai, factor (NtupleInit.f:160-166, task P10.2)
/// The kaon columns (55): the 53 pion/kaon/delta columns + saghai (ntup%sigcm1), factor (ntup%sigcm2).
/// @return the static column table.
const std::vector<NtupleColumn>& KaonNtupleColumns();

/// The pion/kaon/delta columns in Fortran order (53, or 65 with doing_pizero). Columns 1-33 are those of the (e,e'p)
/// layout.
/// @param pizero true with doing_pizero.
/// @return the static column table.
const std::vector<NtupleColumn>& PionNtupleColumns(bool pizero);

/// True if the reaction flags select the pion/kaon layout (NtupleInit.f:101: doing_pion, doing_kaon; delta is not ported).
/// @param f the reaction flags.
/// @return true for the pion/kaon layout.
bool UsesPionLayout(const ReactionFlags& f);

/// results_ntu_write for doing_pion (results_write.f:64-186). `ntu` must have PionNtupleColumns(pizero).size()
/// elements; columns 1-24 keep their previous content for an unknown electron arm (as FillEepNtuple).
/// @param setup the run set-up.
/// @param main the main event record.
/// @param vertex the generated vertex.
/// @param orig the event after radiation.
/// @param recon the reconstructed event.
/// @param[in,out] ntu the ntuple row.
/// @return false for an unknown electron arm.
bool FillPionNtuple(const RunSetup& setup, const MainEvent& main, const Event& vertex, const Event& orig,
                    const Event& recon, std::vector<double>& ntu);

// ------------------------------------------------------------------------------------------- semi layout (P12.2)
constexpr int kSemiNtupleSize = 56;   ///< NtupleInit.f:192-237 (doing_semi, no target field, no pi0)
/// The SIDIS columns (56) in Fortran order; columns 1-33 are those of the (e,e'p) layout.
/// @return the static column table.
const std::vector<NtupleColumn>& SemiNtupleColumns();
constexpr int kSemiPizeroNtupleSize = 68;   ///< doing_semi with doing_pizero: + xcal/ycal/E/p of both photons (P14.3)  ///< doing_semi with doing_pizero: + xcal/ycal/E/p of both photons (P14.3).
/// The SIDIS columns with doing_pizero (68).
/// @return the static column table.
const std::vector<NtupleColumn>& SemiPizeroNtupleColumns();
/// results_ntu_write for doing_semi (results_write.f:187-215). `ntu` must have 56 elements. setup.fortran_compat
/// selects the Fortran pfermi = 0/0 for pfer = 0 (D-051 F3); the default guards the division as the pion block.
/// @param setup the run set-up.
/// @param main the main event record.
/// @param vertex the generated vertex.
/// @param orig the event after radiation.
/// @param recon the reconstructed event.
/// @param[in,out] ntu the ntuple row.
/// @return false for an unknown electron arm.
bool FillSemiNtuple(const RunSetup& setup, const MainEvent& main, const Event& vertex, const Event& orig,
                    const Event& recon, std::vector<double>& ntu);

// ------------------------------------------------------------------------------------------- rho layout (P13.2)
constexpr int kRhoNtupleSize = 59;   ///< NtupleInit.f:192-264 (doing_rho, no target field): the 56 semi + Mrho, Thrho, mmnuc
/// The rho columns (59) in Fortran order.
/// @return the static column table.
const std::vector<NtupleColumn>& RhoNtupleColumns();
/// results_ntu_write for doing_rho (results_write.f:187-228): FillSemiNtuple + ntup%rhomass, ntup%rhotheta, ntup%mmA/1000.
/// @param setup the run set-up.
/// @param main the main event record.
/// @param vertex the generated vertex.
/// @param orig the event after radiation.
/// @param recon the reconstructed event.
/// @param[in,out] ntu the ntuple row.
/// @return false for an unknown electron arm.
bool FillRhoNtuple(const RunSetup& setup, const MainEvent& main, const Event& vertex, const Event& orig,
                   const Event& recon, std::vector<double>& ntu);

// ------------------------------------------------------------------------- polarized-target field (P15.2b)
/// The 8 columns th_tarq, phitarq, beta, phis, phic, betai, phisi, phici that follow phipqi (pion/kaon/delta: 54-61,
/// semi/rho: 57-64) when using_tgt_field (NtupleInit.f:142-159, 239-256).
/// @return the static column table.
const std::vector<NtupleColumn>& TargetFieldNtupleColumns();

/// The layout of a run: the columns and the filler (FillEepNtuple or FillPionNtuple). Throws std::runtime_error for
/// the layouts that are not implemented (semi/rho, kaon, delta, target field).
/// @param setup the run set-up.
/// @return the static column table of the run.
/// @throws std::runtime_error for a layout that is not implemented.
const std::vector<NtupleColumn>& NtupleColumnsFor(const RunSetup& setup);
/// @brief Fill the ntuple row with the filler of the run's layout.
/// @param setup the run set-up.
/// @param main the main event record.
/// @param vertex the generated vertex.
/// @param orig the event after radiation.
/// @param recon the reconstructed event.
/// @param[in,out] ntu the ntuple row (size of NtupleColumnsFor()).
/// @return false for an unknown electron arm.
bool FillNtuple(const RunSetup& setup, const MainEvent& main, const Event& vertex, const Event& orig,
                const Event& recon, std::vector<double>& ntu);

/// @}

}  // namespace simc

#endif  // SIMC_OUTPUT_NTUPLELAYOUT_H
