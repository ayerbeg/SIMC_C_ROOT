// include/simc/optics/ForwardChains.h
//
// Aperture-free forward chains target -> focal plane, built on CosyEngine (task P6.1).
//   hmsChain()  : the project()/transp()/rotate_haxis() sequence of hms/mc_hms.f:201-413
//                 (fixed-slit branch, i.e. using_HMScoll = .false.), every aperture test removed, ending with
//                 project(driftdist(1,12) - z_dip3) to the focal plane (mc_hms.f:413; SIMC hands this drift to
//                 mc_hms_hut instead).
//   shmsChain() : shms/mc_shms.f:322-1047 (use_coll = .true., skip_hb = .false., using_SHMScoll = .false.), no
//                 apertures, ending with project(zd_fp = 307.95) (SIMC projects 307.95 + hcer_1_zentrance and steps
//                 through the hut).
//   cosyChain() : transp(spectr, k) for k = 1..n_classes (pure COSY maps, no project).
// Each step reports the track, the accumulated path and the (xt, yt) the corresponding mc_* aperture test uses
// (rotated/offset coordinates where the Fortran rotates, otherwise xs, ys). A step callback returning false
// stops the chain (hook for the apertures of P6.2/P6.3).
// The mc_* constants are copied verbatim (they also define the reference tables transp_{hms,shms}_chain).
// FORTRAN-DIFF-CANDIDATE (D-015 F5): the SHMS drifts use hard-coded lengths (39.47, 60.68, 307.95) that differ from the
// COSY drift classes (39.4715, 60.68058, 307.9514).
#ifndef SIMC_OPTICS_FORWARDCHAINS_H
#define SIMC_OPTICS_FORWARDCHAINS_H

#include <functional>
#include <string>
#include <vector>

#include "simc/optics/CosyEngine.h"

namespace simc {
namespace optics {

struct ChainStep {
    std::string name;   ///< step name as in the reference tables (coll_ent, q1_in, ..., fp)
    Track track;        ///< COMMON /track/ after the step
    double path{0};     ///< accumulated pathlen (cm)
    double xt{0}, yt{0};   ///< coordinates used by the aperture test of this step
};

/// Called after each step; return false to stop the chain.
using ChainCallback = std::function<bool(const ChainStep&)>;

/// HMS (spectr = 1) aperture-free chain. Returns false if the callback stopped it.
bool hmsChain(CosyEngine& eng, Track& t, double& pathlen, const ChainCallback& cb, DecayState* decay = nullptr);
/// SHMS (spectr = 5) aperture-free chain.
bool shmsChain(CosyEngine& eng, Track& t, double& pathlen, const ChainCallback& cb, DecayState* decay = nullptr);

/// zd that mc_hms / mc_shms pass to transp for a magnet class (0 if the class is not applied with transp).
double mcMagnetZd(int spectr, int cls);
/// Pure COSY chain: transp(spectr, k) for k = 1..n_classes, zd = driftdist for drift classes, else mcMagnetZd.
void cosyChain(CosyEngine& eng, int spectr, Track& t, double& pathlen);

/// Convenience: run a chain and collect every step.
std::vector<ChainStep> collectHmsChain(CosyEngine& eng, Track t);
std::vector<ChainStep> collectShmsChain(CosyEngine& eng, Track t);

}  // namespace optics
}  // namespace simc

#endif  // SIMC_OPTICS_FORWARDCHAINS_H
