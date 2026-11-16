/**
 * @file
 * @ingroup physics
 * @brief simc::rho::RhoModel: IPhysicsModel "rho", the exclusive rho0 cross section (rho_physics.f:peerho).
 *
 * IPhysicsModel plug-in "rho" (task P13.1): the exclusive rho0 cross section of the Fortran SIMC, rho_physics.f:peerho,
 * called by event.f:complete_main (l.1515-1518) for doing_rho:
 *     main%sigcc       = peerho(vertex, main)     (sigccMain / sigcc)
 *     main%sigcc_recon = 1.0                     (sigccRecon)
 *     tgtweight        = targ%Z + targ%N         (applied by the reaction, P13.2)
 * peerho also writes main%w [GeV], main%t [GeV^2], main%pcm, thetacm, phicm, davejac, johnjac and ntup%sigcm, and in
 * the Fortran mode the COMMON efer (main.pfermi.efer, D-057 F6).
 *
 * JSON config (all optional): { "verbose": true }. No data files.
 * Run state read from the set-up: targ%Mtar_struck, targ%M, the reaction flags, fortran_compat (D-057 F6, F20). Unsupported input
 * combinations are refused in init() (rho::checkSupported). docs/FORTRAN_DIFFERENCES.md.
 */
#ifndef SIMC_PHYSICS_RHO_RHOMODEL_H
#define SIMC_PHYSICS_RHO_RHOMODEL_H

#include <sstream>
#include <string>
#include <vector>

#include "simc/physics/rho/RhoPhysics.h"
#include "simc/stages/IPhysicsModel.h"

namespace simc {

struct RunSetup;

namespace rho {

/// @addtogroup physics
/// @{

/// The rho input combinations SIMC++ refuses (each entry is an error message); empty = supported.
/// Checked: doing_rho set (and not switched off by doing_phsp); the hadron arm is 1-5 (HMS, SOS, HRS-R, HRS-L, SHMS):
/// the calorimeter arms (7, 8) make rho_decay STOP ('Unknown spectrometer setup dude!', rho_decay.f:96-98).
/// @param s the run set-up.
/// @return the error messages.
std::vector<std::string> checkSupported(const RunSetup& s);
/// true if a checkSupported problem is a NotPorted feature (HRS hadron arm) rather than an input error.
/// @param problem one message of checkSupported().
/// @return true for a not-ported feature.
bool isNotPortedProblem(const std::string& problem);
/// PeerhoConfig from the set-up (fortran_compat selects the D-057 F6, F20 switches).
/// @param s the run set-up.
/// @return the peerho constants.
PeerhoConfig makePeerhoConfig(const RunSetup& s);
/// GenerateRhoConfig / RhoDecayConfig from the set-up, for the reaction (P13.2) and the tests.
/// @param s the run set-up.
/// @return the generate_rho constants.
GenerateRhoConfig makeGenerateRhoConfig(const RunSetup& s);
/// RhoDecayConfig from the set-up.
/// @param s the run set-up.
/// @return the rho_decay constants.
RhoDecayConfig makeRhoDecayConfig(const RunSetup& s);

/// @brief IPhysicsModel "rho": peerho. Not thread-safe.
class RhoModel : public IPhysicsModel {
public:
    std::string name() const override { return "rho"; }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Checks the inputs (checkSupported: simc::NotPorted for HRS, std::invalid_argument otherwise).
    /// @param ctx the stage context.
    /// @throws simc::NotPorted or std::invalid_argument for unsupported inputs.
    void init(const StageContext& ctx) override;

    /// peerho(vertex, main) on a copy of main (nothing written).
    /// @param vertex the vertex.
    /// @param main the main event record (copied).
    /// @return the cross section [ub/MeV/sr^2].
    double sigcc(const Event& vertex, const MainEvent& main) override;
    /// main%sigcc = peerho(vertex, main): also writes main%w, t, pcm, thetacm, phicm, davejac, johnjac, ntup%sigcm (and
    /// the COMMON efer in the Fortran mode).
    /// @param vertex the vertex.
    /// @param main the main event record (written).
    /// @return the cross section [ub/MeV/sr^2].
    double sigccMain(const Event& vertex, MainEvent& main) override;
    /// event.f:1517 main%sigcc_recon = 1.0
    /// @return 1.
    double sigccRecon(const Event&, const MainEvent&) override { return 1.0; }

    /// @return the peerho constants.
    const PeerhoConfig& config() const { return cfg_; }
    /// @return the peerho constants (modifiable).
    PeerhoConfig& config() { return cfg_; }   ///< test hook
    /// The fort.222 diagnostics (rho_physics.f:252, 262, 293, 390) written so far (not a file in SIMC++).
    /// @return the diagnostics text.
    std::string diagnostics() const { return diag_.str(); }

private:
    bool verbose_{true};
    PeerhoConfig cfg_;
    std::ostringstream diag_;
};

/// @}

}  // namespace rho
}  // namespace simc

#endif  // SIMC_PHYSICS_RHO_RHOMODEL_H
