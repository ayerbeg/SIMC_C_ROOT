// include/simc/driver/Pipeline.h
//
// The plug-in choice of a SIMC++ run (task P7.2, D-026): one registry name (+ its JSON config) per stage.
//
// JSON ("pipeline" section of the run JSON, or the --pipeline file of the CLI). Every key is optional; a stage is
// either a string (the registry name) or an object {"type": name, <config keys>...}. Config keys are merged over
// the default config of that stage.
//   {
//     "rng":      {"type": "mt64", "seed": 12345},
//     "target":   "hallc_cryo",
//     "beam":     "standard",
//     "reaction": "hydrogen_elastic",          // or "phase_space", "deuterium_eep", "nucleus_eep"
//                                              // D/A: {"type": "nucleus_eep", "spectral_function":
//                                              //       {"type": "theory"|"benhar", "file": "fe56.theory"}}
//     "physics":  "sigep",                     // or "deforest", "deforest_cc1", ...
//     "radiator": "simc_rad",
//     "spectrometers": {"hadron": "hms", "electron": {"type": "shms"}},
//     "sinks":    ["root_tree", "run_summary", "gen_hist", "fortran_hist"]   // replaces the whole list
//   }
// Defaults() reproduces the Fortran choices from the run flags (dbase.f / event.f / simc.f dispatch):
//   reaction    doing_phsp -> "phase_space", doing_hyd_elast -> "hydrogen_elastic" (D-009: H + phsp is phase_space,
//               which runs the elastic kinematics with sigcc = 1), doing_deuterium -> "deuterium_eep", doing_heavy ->
//               "nucleus_eep" (P8.2); other reactions are not ported yet (NotPorted).
//   physics     "sigep" for H(e,e'p) (event.f:1451), "deforest" for D/A(e,e'p) (event.f:1455)
//   radiator    "simc_rad" (radc_init_ev runs for every event even with using_rad = 0, event.f:1045)
//   spectrometers  from electron_arm / hadron_arm (MonteCarlo::DefaultArm)
//   sinks       root_tree (only if Nntu > 0, simc.f:69), run_summary, gen_hist, fortran_hist (.hist/.geni text)
//   rng         "mt64" seeded with random_seed (D-002: the RANLUX sequence is not reproduced)
#ifndef SIMC_DRIVER_PIPELINE_H
#define SIMC_DRIVER_PIPELINE_H

#include <string>
#include <vector>

#include "simc/core/Json.h"
#include "simc/core/RunSetup.h"

namespace simc {
namespace driver {

struct StageSpec {
    std::string type;
    Json config = Json::object();
    /// Parse "name" or {"type": name, ...} and merge it over *this (an object without "type" keeps the type).
    void apply(const Json& j, const std::string& what);
    Json toJson() const;
};

struct PipelineConfig {
    StageSpec rng, target, beam, reaction, physics, radiator;
    StageSpec hadron_arm, electron_arm;   ///< type "" = from the arm code
    std::vector<StageSpec> sinks;

    static PipelineConfig Defaults(const RunSetup& setup);
    /// The target type a "pipeline" JSON section selects ("hallc_cryo" if it does not say). Needed BEFORE the set-up
    /// is built: the set-up takes its start-up losses from that target model (ITargetLossProvider).
    static std::string TargetType(const Json& section);
    /// The spectral function a "pipeline" JSON section names for D/A(e,e'p) (reaction.spectral_function: a type name or
    /// {"type", "file"}); empty strings if it does not. Needed BEFORE the set-up (task P8.1): the set-up reads the
    /// theory file and uses it for the generation limits.
    struct SpectralFunctionSpec {
        std::string type, file;
    };
    static SpectralFunctionSpec SpectralFunction(const Json& section);
    /// Apply a "pipeline" JSON section (see the header comment). Unknown keys throw std::invalid_argument.
    void apply(const Json& section);
    Json toJson() const;
};

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_PIPELINE_H
