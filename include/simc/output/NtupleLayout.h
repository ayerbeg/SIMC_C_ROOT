// include/simc/output/NtupleLayout.h
//
// Column layout and filling of the Fortran diagnostic ntuple (tree "h10" after
// util/root_tree/make_root_tree). ROOT-free so it can be unit-tested and reused.
//
// Fortran: NtupleInit.f (tags, l.36-321), results_write.f:results_ntu_write (l.1-271),
//          make_root_tree.f + froot.c (one Float_t branch per tag, tree "h10").
//
// Implemented layouts: the (e,e'p) family, i.e. doing_hyd_elast / doing_deuterium /
// doing_heavy (this includes doing_phsp, which keeps doing_hyd_elast & co.): 46 columns; exclusive pion production
// (doing_pion, task P9.2): 53 columns (65 with doing_pizero). Kaon/delta and semi/rho come with those reactions.
#ifndef SIMC_OUTPUT_NTUPLELAYOUT_H
#define SIMC_OUTPUT_NTUPLELAYOUT_H

#include <array>
#include <string>
#include <vector>

#include "simc/core/Event.h"
#include "simc/core/RunSetup.h"

namespace simc {

struct NtupleColumn {
    const char* name;     ///< tag = ROOT branch name
    const char* source;   ///< C++ expression (for the right-side/left-side assignment see arm_note)
    const char* units;
    const char* fortran;  ///< results_write.f line(s)
};

constexpr int kEepNtupleSize = 46;

/// The 46 (e,e'p) columns in Fortran order.
const std::vector<NtupleColumn>& EepNtupleColumns();

/// True when the electron arm is on the beam-right side (HMS, HRSR, calo on HMS side):
/// results_write.f:64 (SpectrometerId::IsBeamRight). Then the "hs*" columns carry the electron.
/// For the beam-left arms (SOS, HRSL, SHMS, calo on SHMS side; results_write.f:89) the "hs*" columns carry the hadron.
/// "hs" = right-side spectrometer, "ss" = left-side spectrometer.
enum class ArmSide { ElectronRight, ElectronLeft, Unknown };
ArmSide ElectronArmSide(int electron_arm);

/// results_ntu_write for the (e,e'p) family. `ntu` keeps its previous content for
/// columns 1-24 when the electron arm is unknown, exactly like the Fortran static array
/// (results_write.f:118-119, -fno-automatic). Returns false in that case.
bool FillEepNtuple(const RunSetup& setup, const MainEvent& main, const Event& vertex, const Event& orig,
                   const Event& recon, std::array<double, kEepNtupleSize>& ntu);

/// True if the reaction flags select the (e,e'p) ntuple layout (NtupleInit.f:296).
bool UsesEepLayout(const ReactionFlags& f);

// ------------------------------------------------------------------------------------------- pion layout (P9.2)
constexpr int kPionNtupleSize = 53;        ///< NtupleInit.f:101-141 (doing_pion, no target field, no pi0)
constexpr int kPionPizeroNtupleSize = 65;  ///< + 12 pi0 photon columns (NtupleInit.f:165-190)

/// The pion/kaon/delta columns in Fortran order (53, or 65 with doing_pizero). Columns 1-33 are those of the (e,e'p)
/// layout.
const std::vector<NtupleColumn>& PionNtupleColumns(bool pizero);

/// True if the reaction flags select the pion layout (NtupleInit.f:101: doing_pion; kaon/delta are not ported).
bool UsesPionLayout(const ReactionFlags& f);

/// results_ntu_write for doing_pion (results_write.f:64-186). `ntu` must have PionNtupleColumns(pizero).size()
/// elements; columns 1-24 keep their previous content for an unknown electron arm (as FillEepNtuple).
bool FillPionNtuple(const RunSetup& setup, const MainEvent& main, const Event& vertex, const Event& orig,
                    const Event& recon, std::vector<double>& ntu);

/// The layout of a run: the columns and the filler (FillEepNtuple or FillPionNtuple). Throws std::runtime_error for
/// the layouts that are not implemented (semi/rho, kaon, delta, target field).
const std::vector<NtupleColumn>& NtupleColumnsFor(const RunSetup& setup);
bool FillNtuple(const RunSetup& setup, const MainEvent& main, const Event& vertex, const Event& orig,
                const Event& recon, std::vector<double>& ntu);

}  // namespace simc

#endif  // SIMC_OUTPUT_NTUPLELAYOUT_H
