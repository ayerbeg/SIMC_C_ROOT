# SIMC++ architecture (translation branch)

Status: design v1, 2026-09-29 (task P1.1). This document replaces the older status and architecture notes
in `docs/`. Those notes are historical and unverified (D-004).

## 1. Principles

1. **The Fortran is the specification of the physics.** Every component reproduces the algorithm of
   the Fortran routine(s) it replaces and names them in a comment:
   `// Fortran: event.f:complete_ev (l.432-1055)`. Where the Fortran looks wrong, we reproduce it, flag it with
   `// FORTRAN-QUIRK(D-xxx): …` and record it in project_control/DECISIONS.md. Known cases: simc.f:237, `lambda`
   overwritten, and the `doing_phsp` caveats.
2. **Everything is pluggable.** The event chain is a sequence of stages, and each stage is an
   abstract interface with a name-keyed registry. The JSON selects an implementation by name and passes
   it its own config subtree. Stages are not a closed list. A new kind of stage (for example a detector-response
   model or an event filter) is added by defining an interface and inserting it into the driver's chain.
3. **Validation is statistical (L2).** Deterministic pieces also get function-level tests (L1)
   against Fortran tables. The random-number sequence need not match; the distributions must.
4. **Explicit data flow, no hidden global state.** Components receive the parameters, the event and the RNG
   explicitly. They may cache derived setup data (for example COSY matrices) created at `init()`.
5. **Units follow the Fortran internally**: MeV, MeV/c, cm, rad. Fortran conventions (delta in %, slopes
   dx/dz, mrad in some outputs) are kept on the event structs and documented per field. Converting to the
   output units is the output module's job.

## 2. Data model (mirrors structures.inc)

| C++ (`include/simc/core/`) | Fortran | notes |
|---|---|---|
| `SimcParameters` | all 106 variables registered in dbase.f `regallvars` + `extra_dbase_file` | Grouped like the Fortran (`gen`, `targ`, `spec.e/p`, `cuts`, `SPedge`, `slop`, flags…). The defaults are the **true** Fortran defaults (mostly 0/false; see reference/functions/_inputs/inp_variables.csv). |
| `RunSetup` | quantities derived in init.f / dbase.f (edges, slops, Mh, targ derived, radc init, central event) | Built once by `SetupBuilder` from `SimcParameters`. |
| `ArmFull`, `Event` | `arm_full`, `event` | `vertex`, `orig`, `recon` are `Event`s. |
| `EventTarget`, `MainEvent` | `event_target`, `event_main` | weight, jacobian, sigcc, SP/RECON/FP per arm… |
| `CentralEvent` | `event_central(_rad)` | Output of `calculate_central`. |
| `Limits`, `Range`, `Cuts` | `cutstype`, `rangetype`, `arm_*` | |

## 3. Stage interfaces (`include/simc/stages/`)

Each interface has `name()`, `configure(const Json&, const SimcParameters&)` and `init(RunSetup&)`, plus
the stage-specific methods below.

| Interface | Responsibility | First implementations (Fortran origin) |
|---|---|---|
| `IBeam` | beam energy, spread, raster, offsets for the vertex | `StandardBeam` (event.f:generate beam part, l.150-210) |
| `ITarget` | geometry, materials, energy loss and multiple scattering along the path (`trip_thru_target`, `target_musc`, extreme trips for limits) | `HallCCryoTarget` (target.f, enerloss_new.f); later `PolarizedTarget` (trg_track.f) |
| `IReaction` | generates the physics vertex: variables, limits, phase-space volume (genvol), Jacobian, `complete_ev`/`complete_recon_ev` kinematics | `HydrogenElastic`, `PhaseSpace` (doing_phsp); later `AeepSpectralFunction`, `PionElectroproduction`, `Kaon`, `Delta`, `Semi`, `Rho`… |
| `IPhysicsModel` | cross section `sigcc(vertex)` and `sigcc_recon` for a reaction | `DeForestCC1/CC2`, `SigEP` (physics_proton.f), `Unity` (phsp) |
| `IRadiator` | radiative corrections: `radc_init`, `radc_init_ev`, `generate_rad`, `peaked_rad_weight` | `MoTsaiPeaking` (radc.f, brem.f), `NoRadiation` |
| `ISpectrometer` | forward transport through the magnets and apertures, hut/detectors, reconstruction; returns ok + **stop reason**, FP and recon | `HMS`, `SHMS` on a shared `CosyEngine`; later `SOS`, `HRSL`, `HRSR`, `Calo` |
| `ICollimator` / `IAperture` (inside spectrometers) | reusable aperture elements (octagon, slit, dipole boundary…) with their own stop counters | `HmsCollimator`, `ShmsOctagon`, `HitDipole`… |
| `IMaterialModel` | energy loss (most probable / Landau) and multiple scattering (`musc`, `musc_ext`) | `SimcEnerloss`, `SimcMusc` |
| `IOutputSink` | per event and end of run | `RootTreeSink` (make_root_tree branch names), `HistReportSink` (.hist, simc.f:report), `GenHistSink` (.gen/.geni) |
| `IRandom` | uniform/gauss1 | `Mt64Random` (default), with `gauss1(nsig_max)` reproducing gauss1.f truncation |

Registration:
```cpp
SIMC_REGISTER(ISpectrometer, "hms", HMS);   // static self-registration
auto spec = Registry<ISpectrometer>::create(json["spectrometers"]["hadron"]["type"]);
```
A new plug-in is a new `.cpp` with a `SIMC_REGISTER` line. The driver doesn't change.

## 4. Driver (`SimcRun`, follows simc.f)

Implemented in P7.2; the full description (pipeline JSON, EventState, outputs, validation) is in [DRIVER.md](DRIVER.md).

```
configure: JSON (or .inp → JSON converter) → SimcParameters → RunSetup (init.f equivalents)
calculate_central()                                   // simc.f:1143
loop until |ngen| successes (ngen>0) or tries (ngen<0) // simc.f:169, 341-350
   reaction.generate(beam, target, radiator, rng)     // event.f:generate + complete_ev
   montecarlo(): per arm target exit + spectrometer.transport() // simc.f:1310
   reaction.completeRecon()                           // event.f:complete_recon_ev
   completeMain(): sigcc, radiative weight, jacobians, weight   // event.f:complete_main
   stop-reason and limits accounting; sinks.fill()
normalisation: luminosity, genvol, normfac, ntried    // simc.f:360-400
sinks.finish()  → .root, .hist, .gen, .geni
```
**Acceptance mode** (D-008): `doing_phsp` selects the `PhaseSpace` reaction and the `Unity` model, forces
radiation off and normfac = 1, exactly as dbase.f:227-236, event.f:709-713/1132-1135/1446-1448 and simc.f:397 do.
With `ngen<0` it counts accepted events out of |ngen| tries. The per-aperture stop counters are reported for both arms.

## 5. Input

- The JSON schema holds every `SimcParameters` field (D-003) plus C++-only settings (plug-in choice, output paths, RNG).
- `tools/inp2json`: Fortran `.inp` (+ `extra_dbase_file`) → JSON. It accepts `%` and obsolete `.` names,
  is case-insensitive like CTP, and **warns on unknown names**. A completeness test fails if any `regallvars`
  variable has no JSON field.

## 6. Output parity

- The ROOT tree `h10` has the 46 `make_root_tree` branches (hsdelta … theta_p), with the same names and units, as
  Float_t like the Fortran tool. Extra C++ branches may be added under a prefix.
- The `.hist` equivalent carries every quantity `report()` prints that the comparison tool needs
  (Ngen, Ntried, Ncontribute, wtcontr, sigcc central/average, luminosity, genvol, normfac, stop counters…).

## 7. Layout

```
include/simc/core/        parameters, event structs, setup, registry, random, constants
include/simc/stages/      the interfaces
src/beam/ src/target/ src/reactions/ src/physics/ src/radiation/ src/spectrometers/{common,hms,shms,…} src/output/
src/driver/               SimcRun, main
tools/                    inp2json, run_simc.sh
validation/               comparison tool (pulls, χ²/KS, PDFs)
tests/unit/ tests/l1/     CTest; L1 runners read reference/functions/<routine>/ tables
data/                     COSY files, tables, JSON configs
```

## 8. Migration from the current code (P0.3 audit)

- **Kept** (moved into the new layout): SimcConstants, CrossSection proton part (sigep, deForest, fofa, sigmott),
  ExternalRadiation, the CosyMatrix parser (seed of `CosyEngine`), the RandomGenerator core.
- **Refactored**: event structs, OutputManager → sinks, EnergyLoss, RadiativeCorrections/Bremsstrahlung (bugs fixed against the Fortran), HMS/SHMS (on `CosyEngine`).
- **Rewritten**: main/driver, EventGenerator → Beam/Target/Reaction, MonteCarloTransport, ConfigManager, tests.
- **Removed** once superseded: SpectrometerOptics, Kinematics, pion/kaon stubs, the duplicate RadiativeCorrections,
  `*.backup`/`*.v1`, stale docs and tests. They remain in git history.
- SOS/HRS classes: left in place but out of the build until Phase 14.
