// include/simc/output/NtupleLayout.h
//
// Column layout and filling of the Fortran diagnostic ntuple (tree "h10" after
// util/root_tree/make_root_tree). ROOT-free so it can be unit-tested and reused.
//
// Fortran: NtupleInit.f (tags, l.36-321), results_write.f:results_ntu_write (l.1-271),
//          make_root_tree.f + froot.c (one Float_t branch per tag, tree "h10").
//
// Implemented layout: the (e,e'p) family, i.e. doing_hyd_elast / doing_deuterium /
// doing_heavy (this includes doing_phsp, which keeps doing_hyd_elast & co.): 46 columns.
// The pion/kaon/delta, semi/rho layouts will be added with those reactions.
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
/// results_write.f:64 (electron_arm = 1, 3, 7). Then the "hs*" columns carry the electron.
/// For electron_arm = 2, 4, 5, 6, 8 (results_write.f:89) the "hs*" columns carry the hadron.
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

}  // namespace simc

#endif  // SIMC_OUTPUT_NTUPLELAYOUT_H
