# SIMC++ architecture (v2)

Status: v2, task P7.4 (2026-10-04). It replaces v1 (P1.1). The driver itself is described in
[DRIVER.md](DRIVER.md), the output files in [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md), the input mapping in
[INP_JSON_MAPPING.md](INP_JSON_MAPPING.md), and how to add plug-ins in [EXTENDING.md](EXTENDING.md).

## 1. Principles

1. **The Fortran is the specification of the physics** (D-001). Every ported routine reproduces the algorithm, the
   operation order (bitwise results) and the random-number calls of the Fortran routine it replaces, and cites it:
   `// Fortran: event.f:complete_ev (l.432-1055)`. Fortran behaviour that looks wrong is reproduced, marked
   `// FORTRAN-QUIRK(D-xxx)` and recorded in project_control/DECISIONS.md (D-026: the default is always the Fortran
   behaviour; a corrected behaviour may later become an opt-in).
2. **A modern, object-oriented framework around faithful kernels** (D-026). The Fortran-faithful code ("kernels")
   is the numerical core. The rest — stage interfaces, plug-ins, driver, I/O — is ordinary C++17: explicit inputs
   and outputs, RAII and `unique_ptr` ownership, value types for events, const-correct interfaces, no global state.
3. **Everything is a plug-in** chosen by name from the JSON pipeline (§4). The driver depends only on interfaces.
4. **Validation is layered** (D-002, §7): L1 tables per routine, bitwise replays of captured Fortran runs,
   before/after output identity for refactors, and the statistical L2 comparison with the Fortran.
5. **One unit convention** (`core/Units.h`): MeV, MeV/c, MeV^2, cm, rad, g/cm^2, delta in %, slopes. Input units
   (.inp/JSON: degrees, mrad, mg/cm^2) are converted once by the set-up; output units (ntuple GeV, mrad) by the
   output code.

## 2. Layers

```
 tools        simc++ (CLI), inp2json, hist2json, simc_root_identity        validation/  simc_compare (L2)
 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 driver       SimcRun (simc.f main), MonteCarlo (simc.f:montecarlo), Pipeline, Central, RunReport, ReportSink
 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 io           input/  InpReader (CTP .inp), EffectiveValues      output/  sinks, NtupleLayout, SummaryJson
 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 plug-ins     material/HallCCryoTarget, SimcMaterial   radiation/SimcRadiator   reactions/StandardBeam,
 (stages)     HydrogenElastic   physics/proton models   spectrometers/{hms,shms}/*Spectrometer   output sinks
              examples/ideal_box
 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 kernels      radiation/kernel (radc.f, brem.f, init.f radc_*)   spectrometers/*/kernel (mc_hms.f, mc_shms.f ...)
              optics (COSY engine: transp.f, project.f, recon)   setup (dbase.f/init.f -> RunSetup)
              reactions/EventChain (complete_main, limits_update, transform_to_cm)   physics/proton (sigep, deForest)
 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 stages       include/simc/stages: IStage, StageContext and the stage interfaces (§3)
 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 core         SimcParameters, RunSetup, Event/MainEvent, RunSummary, Registry, Errors, Units, constants, random
```

A layer uses only the layers below it. In particular the driver uses only `core`, `stages`, `setup`,
`reactions/EventChain` and `io`. It includes no plug-in header and casts no stage to a concrete class.

**Kernel visibility rule.** Fortran-mirror state lives in headers under a `kernel/` directory
(`include/simc/radiation/kernel/`, `include/simc/spectrometers/{hms,shms}/kernel/`): COMMON-block structs
(`RadcState`, `HmsCounters`, `ShmsCounters`, `/track/`), argument-list mirrors (`McHmsArgs`, `McShmsArgs`),
1-based arrays and SAVEd locals. Only the module's own sources, its plug-in `.cpp` and its L1 tests include them,
and they are not installed. A plug-in header exposes its kernel only as an opaque `kernel()` accessor (pimpl),
for the L1 tests. The optics (COSY) engine and the set-up builder are kernel libraries shared by several plug-ins
or by the driver. Their headers are public, but they are Fortran-structured (spectrometer ids, 1-based classes).

## 3. Stage interfaces (`include/simc/stages/`)

Life cycle: `Registry<I>::create(name)` → `configure(json, params)` → `init(StageContext)` → per-event calls.

| Interface | Per-event role (Fortran origin) | Built-in plug-ins |
|---|---|---|
| `IRandom` | `uniform()` = grnd, `gauss1(nsig)` = gauss1.f | `mt64` (`ScriptedRandom` for replays) |
| `ITarget` (+ `ITargetLossProvider`) | vertex z, `tripThruTarget` → {eloss, radlen}, `multipleScattering` → {dyptar, dxptar}; start-up `targetInit`, `extremeLosses` (target.f, init.f:target_init) | `hallc_cryo` |
| `IBeam` | beam spread, raster, vertex position, incoming energy (event.f:generate l.150-212) | `standard` |
| `IReaction` | `generate`, `completeEv`, `completeRecon`, `genVolume` (event.f) | `hydrogen_elastic`, `phase_space` |
| `IRadiator` | `initEvent` (radc_init_ev), `radiate` (generate_rad); read-back `currentValues`, `photonEnergies`, `hadronRadiates` | `simc_rad` |
| `IPhysicsModel` | `sigcc`, `sigccRecon` (physics_*.f, called by complete_main) | `sigep`, `deforest`, `deforest_cc1/cc2/cc1_onshell` |
| `ISpectrometer` | `transport(SpecInput) → SpecResult`, `stopCounters()` (mc_*.f) | `hms`, `shms`, `ideal_box` (example) |
| `IMaterialModel` | energy loss, musc, musc_ext (enerloss_new.f, musc*.f) | `simc_material` (inside `hallc_cryo`) |
| `IOutputSink` | `beginRun`, `fill` (every try), `endRun` (+ optional `driver::IReportSink`) | `root_tree`, `run_summary`, `gen_hist`, `fortran_hist` |

`StageContext` = {`const RunSetup& setup`, `ITarget* target`, `IRadiator* radiator`}. It is implicitly
constructible from a `RunSetup`, so `stage.init(setup)` works for stages without collaborators.

Errors (`core/Errors.h`): `simc::FortranStop` (a Fortran STOP reached), `simc::NotPorted` (a Fortran feature not
ported yet), both derived from `simc::Error`. Misuse is `std::logic_error`, bad configuration `std::invalid_argument`.

## 4. Plug-in wiring

```cpp
SIMC_REGISTER(simc::ISpectrometer, "hms", simc::hms::HmsSpectrometer);   // in the plug-in's .cpp
```

Registration happens during static initialisation, so plug-ins are compiled into CMake OBJECT libraries that
the executables link directly (`SIMC_DRIVER_LIBS`, `SIMC_EXAMPLE_PLUGINS`). The JSON `pipeline` section (or
`simc++ --pipeline FILE`) names the plug-in of each stage. A stage is `"name"` or `{"type": "name", <config>}`, and
the config is merged over the defaults:

```json
"pipeline": {"rng": {"type": "mt64"}, "target": "hallc_cryo", "beam": "standard", "reaction": "hydrogen_elastic",
             "physics": "sigep", "radiator": "simc_rad",
             "spectrometers": {"hadron": "hms", "electron": {"type": "shms"}},
             "sinks": ["root_tree", "run_summary", "gen_hist", "fortran_hist"]}
```

`PipelineConfig::Defaults(setup)` reproduces the Fortran choices from the input flags (reaction from `doing_*`, model
from the target, arms from `electron_arm`/`hadron_arm`). `simc++ --print-pipeline` prints the resolved pipeline.

`SimcRun::initialise` (simc.f l.53-160):

1. creates a set-up instance of the pipeline target and builds `RunSetup` with `BuildRunSetup(params, opts,
   target)` (dbase.f / init.f; the target supplies target_init and extreme_trip_thru_target);
2. resolves the pipeline, creates `rng`, `target`, `radiator`, then `beam`, `reaction`, `physics` with
   `StageContext{setup, target, radiator}`, the `MonteCarlo` with its arm spectrometers, and the sinks;
3. computes the central event (`calculateCentral`, which also takes central%rad from `IRadiator::currentValues`).

## 5. Data flow per event (one try, `SimcRun::step` = simc.f l.169-352)

```
EventState ev {main, vertex, orig, recon}       ONE object for the whole run (Fortran stale-value semantics, D-024)
 resetForTry                                    simc.f:174-183 (ntup fields)
 IBeam::generate(main, vertex)                  beam x/y, raster, z (ITarget::generateZ), Eloss(1) (tripThruTarget), Ein
 IReaction::generate(main, vertex, orig)        angles, completeEv (tripThruTarget e/p, IRadiator::initEvent),
                                                IRadiator::radiate(..., completeEv) or orig = vertex
 MonteCarlo::run(orig, main, recon)             per arm: SP quantities, ITarget::multipleScattering, TRANSPORT frame,
                                                ISpectrometer::transport → RECON/FP, recon kinematics,
                                                most-probable energy-loss correction (tripThruTarget)
 IReaction::completeRecon(recon)                complete_recon_ev
 reactions::completeMain(...)                   weight: IPhysicsModel::sigcc/sigccRecon, jacobian, gen_weight
 pass_cuts, hard_cuts                           simc.f:226-243
 IOutputSink::fill(main, vertex, orig, recon, success)   for every try
 totals, limitsUpdate (IRadiator::photonEnergies, hadronRadiates)
finish: Normalise → RunReport → RunSummary → sinks endRun / writeReport (.root, .summary.json/.txt, .gen.root,
        .hist/.geni/.gen)
```

The order of these calls and of every random-number draw is the Fortran order. That is what makes the bitwise
replays (§7) possible.

## 6. Fortran correspondence

* Each ported function or block cites its origin, `// Fortran: file.f:routine (l.a-b)`, and each file header lists
  its routines. project_control/TRANSLATION_MAP.md maps the Fortran routines to the C++ ones.
* Kernels keep the Fortran structure: the same branches, the same operation order (`x**2` → `x*x`, `x**n` →
  `optics::powi`, D-024), the Fortran literals (for example brem.f works in GeV). COMMON blocks and SAVEd locals
  become members of the object that represents "one Fortran process" (one `McHms`, one `SimcRadiator` ...), and
  the zero-initialised COMMON semantics are reproduced by value-initialisation.
* The stage interfaces carry the physics content of the Fortran argument lists, without the Fortran-isms: return
  values instead of output arguments, `std::optional` for "written or not" (`SpecResult::resolution_factor`),
  0-based tails (`RadiationValues`, `photonEnergies`).
* `FORTRAN-QUIRK` marks reproduced odd behaviour. Persistent quirks (stale values) live in an object, never in a
  global: `EventState` (stale event fields), `MonteCarlo::tmpfact_`, `HydrogenElastic::krel_`,
  `HallCCryoTarget::thCan_`, the kernels' SAVEd members.

## 7. Validation and tests (CTest labels)

| label | what | reference |
|---|---|---|
| `unit` | parameters, .inp reader, sinks, summary/hist2json, stage headers and registry | reference/runs .inp/.hist |
| `l1` | routine tables, bitwise or with stated tolerance | reference/functions/<routine>/*.tsv |
| `replay` | the complete MonteCarlo and SimcRun driven by a captured RANLUX stream: bitwise tree/.bin/.hist | reference/functions/event_capture, runs/p4f1_*_cap |
| `cli` | `simc++` with .inp vs JSON inputs, `inp2json`, plug-in selection by `--pipeline` | — |
| `regression` | before/after output identity (`tests/regression`): all h10 branches and entries bitwise, text outputs equal apart from time stamps | SIMC++_Project/baselines/f9cea51_{rad,norad,phsp}_s4242 |
| `examples` | example plug-ins run through the full driver | — |
| `slow` | long distribution tests | — |

L2 (statistical, D-002): `validation/simc_compare` compares a C++ run with the Fortran A/B runs (.hist pulls,
weighted χ²/KS per tree variable, PDF). Results are in reference/comparisons/.

## 8. Source layout

```
include/simc/core/            parameters, RunSetup, events, Registry, Errors, Units, constants, random sources
include/simc/stages/          the interfaces, StageContext
include/simc/<module>/        plug-in classes and kernel libraries; <module>/kernel/ = private Fortran mirrors
src/<module>/                 one CMake target per module (simcfw_<module>; OBJECT if it registers plug-ins)
src/driver/                   simcfw_driver (+ simcfw_driver_sinks): SimcRun, MonteCarlo, Pipeline, report
tools/                        simc_main.cpp (simc++), inp2json, hist2json, check_all.sh
examples/                     example plug-ins (ideal_box) and their tests
tests/{unit,l1,regression}/   CTest suites (labels above)
validation/                   simc_compare (own CMake project)
data/                         COSY maps (byte copies of the Fortran files), tables, configs (data/config/*.json)
```

## 9. Known limitations of the design (to be lifted with the ports that need them)

* `reactions::completeMain` (event.f:complete_main) assembles the weight in the driver and has only the
  H(e,e'p)/phase-space branch. When the other reactions are ported (P8+), the reaction-specific part should move
  behind `IReaction`.
* The arm codes `electron_arm`/`hadron_arm` still drive Fortran logic outside the spectrometer: the window
  materials in `hallc_cryo`, the `.geni` blocks and the default plug-in choice.
* The set-up (`BuildRunSetup`) is one large port of dbase.f/init.f with reaction-specific branches. It is not
  split per stage.
* Out-of-tree plug-ins: the install rules install the programs, the public headers and the data, but no CMake
  package config yet. Plug-ins are built in-tree (EXTENDING.md).
