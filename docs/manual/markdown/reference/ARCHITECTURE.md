# SIMC++ architecture (v2)

Status: v2, task P7.4 (2026-10-04); updated in P8.2 (reaction-owned weight, D/A(e,e'p), spectral functions), P15.3/P15.4 (layout, package, threads) and P15.D1 (labels, limitations). It replaces v1 (P1.1). The developer manual ([manual/developer/index.md](../developer/index.md)) is the narrative entry point; this page is its short reference. The driver itself is described in
[DRIVER.md](DRIVER.md), the output files in [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md), the input mapping in
[INP_JSON_MAPPING.md](INP_JSON_MAPPING.md), and how to add plug-ins in [EXTENDING.md](EXTENDING.md).

## 1. Principles

1. **The Fortran is the specification of the physics** (D-001). Every ported routine reproduces the algorithm, the
   operation order (bitwise results) and the random-number calls of the Fortran routine it replaces, and cites it:
   `// Fortran: event.f:complete_ev (l.432-1055)`. Since D-031 (P7.5) clear Fortran bugs are CORRECTED by default and
   marked `// FORTRAN-BUG(D-xxx): original <file:line>: <code>; corrected because ...`; the run option `fortran_compat`
   (`--fortran-compat`, JSON `"run": {"fortran_compat": true}`, `RunSetup::fortran_compat`) restores the exact Fortran
   behaviour for the bitwise L1/replay/regression tests. Ambiguous cases keep the Fortran behaviour and are marked
   `FORTRAN-DIFF-CANDIDATE(D-xxx)` with the question; the open acceptance-mode items `FORTRAN-OPEN(D-031 flag 1)`;
   documented faithful oddities `FORTRAN-QUIRK(D-xxx)`. Every item is listed in [FORTRAN_DIFFERENCES.md](../../../translation/markdown/FORTRAN_DIFFERENCES.md).
2. **A modern, object-oriented framework around faithful kernels** (D-026). The Fortran-faithful code ("kernels")
   is the numerical core. The rest — stage interfaces, plug-ins, driver, I/O — is ordinary C++20: explicit inputs
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
 (stages)     HydrogenElastic, NuclearEep (on EepReactionBase), Pion/Kaon/Delta/Sidis/RhoProduction
              physics/{proton,nuclear,pion,kaon,delta,sidis,rho} models and spectral functions
              spectrometers/{hms,shms,sos,hrs,calo}/*Spectrometer   output sinks   examples/ideal_box
 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 kernels      radiation/kernel (radc.f, brem.f, init.f radc_*)   spectrometers/*/kernel (mc_hms.f, mc_shms.f ...)
              optics (COSY engine: transp.f, project.f, recon; TargetField: trg_track.f)   setup (dbase.f/init.f -> RunSetup)
              reactions/EventChain (assembleWeight, limits_update, transform_to_cm)   physics/proton (sigep, deForest)
              physics/nuclear (sf_lookup.f SfTable, theory-file SF weight)
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
| `IRandom` | `uniform()` = grnd, `gauss1(nsig)` = gauss1.f | `mt64`, `ranlux` (the Fortran RANLUX/grnd stream, bit for bit, P13.2) (`ScriptedRandom` for replays) |
| `ITarget` (+ `ITargetLossProvider`) | vertex z, `tripThruTarget` → {eloss, radlen}, `multipleScattering` → {dyptar, dxptar}; start-up `targetInit`, `extremeLosses` (target.f, init.f:target_init) | `hallc_cryo` |
| `IBeam` | beam spread, raster, vertex position, incoming energy (event.f:generate l.150-212) | `standard` |
| `IReaction` | `generate`, `completeEv`, `completeRecon`, `completeMain` (the event weight), `genVolume` (event.f, simc.f) | `hydrogen_elastic`, `phase_space`, `deuterium_eep`, `nucleus_eep`, `pion_production` (P9.2), `kaon_production` (P10.2), `delta_production` (P11.1), `sidis_production` (P12.2), `rho_production` (P13.2) |
| `ISpectralFunction` | `weight(vertex)` = main%SF_weight (event.f:complete_main l.1397-1430); a component of the D/A(e,e'p) reactions | `theory` (init.f theory files), `benhar` (sf_lookup.f) |
| `IRadiator` | `initEvent` (radc_init_ev), `radiate` (generate_rad); read-back `currentValues`, `photonEnergies`, `hadronRadiates` | `simc_rad` |
| `IPhysicsModel` | `sigcc`, `sigccRecon` (physics_*.f, called by complete_main) | `sigep`, `deforest`, `deforest_cc1/cc2/cc1_onshell`, `pion` (physics_pion.f:peepi; `sigccMain` also writes main%thetacm ... and ntup%sigcm*), `kaon` (physics_kaon.f:peeK, Saghai + factorized models, CERNLIB fint; P10.1), `delta` (physics_delta.f:peedelta: phase space × flux, NO cross section, D-048 F4; P11.1), `sidis` (semi_physics.f:peepiX: CTEQ5M PDFs × Bosted π / fDSS K fragmentation functions × F1F2IN21 (Christy 2021); `sigccMain` writes main%davejac, ntup%sigcm, xfermi, decdist; P12.1), `rho` (rho_physics.f:peerho; `sigccMain` writes main%w, t [GeV units], pcm, thetacm, phicm, davejac, johnjac, ntup%sigcm; P13.1; kernels `rho::generateRho`, `rho::rhoDecay` called by `rho_production`) |
| `ISpectrometer` | `transport(SpecInput) → SpecResult`, `stopCounters()` (mc_*.f) | `hms`, `shms`, `sos` (P14.1), `hrsr`, `hrsl` (P14.2), `calo` (P14.3; the π0 two-photon logic is `driver::PizeroCalorimeter`), `ideal_box` (example) |
| `IMaterialModel` | energy loss, musc, musc_ext (enerloss_new.f, musc*.f) | `simc_material` (inside `hallc_cryo`) |
| `IOutputSink` | `beginRun`, `fill` (every try), `endRun` (+ optional `driver::IReportSink`) | `root_tree`, `run_summary`, `gen_hist`, `fortran_hist` |

`StageContext` = {`const RunSetup& setup`, `ITarget* target`, `IRadiator* radiator`, `IPhysicsModel* physics`}. It is
implicitly constructible from a `RunSetup`, so `stage.init(setup)` works for stages without collaborators. Since P8.2
(D-030) the reaction owns the event weight: `IReaction::completeMain` sets main%SF_weight from its own spectral function
and assembles the weight with the run's physics model (`StageContext::physics`) through the kernel
`reactions::assembleWeight` (event.f:complete_main l.1432-1563).

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
   The physics model is created before the reaction and handed to it (`StageContext::physics`); the D/A(e,e'p)
   reaction creates its `ISpectralFunction` (JSON `reaction.spectral_function`) in its `init()`;
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
 IReaction::completeMain(...)                   weight: ISpectralFunction (D/A), IPhysicsModel::sigcc/sigccRecon,
                                                Coulomb factor, jacobian, gen_weight
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
* `FORTRAN-QUIRK` marks reproduced odd behaviour (FORTRAN_DIFFERENCES.md Category D); `FORTRAN-BUG`, `FORTRAN-DIFF-CANDIDATE`
  and `FORTRAN-OPEN` Categories A, B, C. The `fortran_compat` flag is read once into `RunSetup` and reaches every stage
  through `StageContext::setup`; kernels get it explicitly (function argument, `ExtremeTripArgs`, `DeForestState`,
  `CosyEngine::setFortranCompat`). Persistent quirks (stale values) live in an object, never in a
  global: `EventState` (stale event fields), `MonteCarlo::tmpfact_`, `EepReactionBase::krel_`, `nuclear::SfTable` (the SAVEd sf_lookup locals),
  `HallCCryoTarget::thCan_`, the kernels' SAVEd members. Optional diagnostics of the transport kernels (`optics/TransportTrace.h`:
  every project/transp step, the decay record, the collimator calls) are attached by the L1 tests and never change the numerics.

## 7. Validation and tests (CTest labels)

| label | what | reference |
|---|---|---|
| `unit` | parameters, .inp reader, sinks, summary/hist2json, stage headers and registry | reference/runs .inp/.hist |
| `l1` | routine tables, bitwise or with stated tolerance | reference/functions/<routine>/*.tsv (P8: theory_init, sf_lookup*, generate_em, eep_state, *_eepA; P6.4: mc_{hms,shms}_{coll,decay}, transp_decay_*, project_decay, pion_coll_absorb, bitwise incl. steps.txt.xz / coll_calls.txt; P9.1: pion cross sections, transform_to_cm_pion, pizero_decay_fixedrng; P9.2: pion_state, generate_pion_fixedrng, complete_ev_pion, complete_recon_ev_pion, complete_main_pion, decay_fixedrng, mc_decay_fixedrng; P10.1/P10.2: saghai_tables, sig_factorized, eekeek/eekeeks/phspwght/peek (+ `_fintabi` for the corrected fint path), transform_to_cm_kaon, kaon_state, generate_kaon_fixedrng, complete_ev_kaon, complete_recon_ev_kaon, complete_main_kaon, decay_kaon_fixedrng, mc_decay_kaon_fixedrng; P11.1: delta_state, generate_delta_fixedrng, complete_ev_delta, complete_recon_ev_delta, complete_main_delta, peedelta, transform_to_cm_delta; P12.1: sidis/ctq5pdf (+ ctq5_state), fdss, f1f2in21, slacemcfit_rhad_global, peepix, `l1_sidis_data_files`; unit `sidis_tables`, `sidis_sanity`; P12.2: sidis/generate_semi_fixedrng, complete_ev_semi, complete_recon_ev_semi, complete_main_semi; P13.1/P13.2: rho/peerho, generate_rho_fixedrng, rho_decay_fixedrng, rho_state, generate_rho_branch_fixedrng, complete_ev_rho, complete_recon_ev_rho, complete_main_rho, `l1_ranlux_stream_rho_h`) |
| `replay` | the complete MonteCarlo and SimcRun driven by a captured RANLUX stream: bitwise tree/.bin/.hist; every mc_hms/mc_shms call of the P6.F2 pion captures (9464 calls, `l1_replay_p6f2_pi_*_cap`); the two P9.F1 pion captures (π+ H, π− D, decay ON): full run (`l1_simcrun_replay_pion_*`), montecarlo (`l1_montecarlo_replay_pion_*`) and per routine (`l1_replay_event_capture_pion_*`); the two P10.F1 kaon captures (K+Λ H, K+Σ− D, decay ON): `l1_simcrun_replay_kaon_*`, `l1_montecarlo_replay_kaon_*`, `l1_replay_event_capture_kaon_*`; the two P11.F1 Delta captures (H(e,e'p)π rad, D1/D2): `l1_simcrun_replay_delta_*`, `l1_montecarlo_replay_delta_*`, `l1_replay_event_capture_delta_*`; the two P12.F1 SIDIS captures (H/D(e,e'π+)X rad): `l1_simcrun_replay_semi_{h,d}`, `l1_montecarlo_replay_semi_{h,d}`, `l1_replay_event_capture_semi_{h,d}`; the two P13.F1 rho captures (H/D(e,e'π+)[ρ] rad, two-pass format): per routine and montecarlo over the full-logged tries, and the full SimcRun with the whole stream regenerated by the RANLUX port (`l1_simcrun_replay_rho_{h,d}`: every one of the 2.9e6 tries, tree, .hist/.geni/.gen) | reference/functions/event_capture{,_eepA,_pion}, runs/p4f1_*_cap, runs/p8f1_{d,fe56}_rad_cap, runs/p6f2_pi_{hms,shms}_cap (call files in fortran_harness/spectrometers_coll_decay/work) |
| `cli` | `simc++` with .inp vs JSON inputs, `inp2json`, plug-in selection by `--pipeline` | — |
| `regression` | before/after output identity (`tests/regression`): all h10 branches and entries bitwise, text outputs equal apart from time stamps | CURRENT (P15.1, D-054): every family regenerated once as SIMC++_Project/baselines/55cf7a0_glibc243_<suffix> (`SIMC_*_BASELINE_COMMIT`), plus the calorimeter family; history: SIMC++_Project/baselines/f9cea51_{rad,norad,phsp}_s4242 (`--fortran-compat`), 08812e6_corrected_{rad,norad,phsp}_s4242 (default); D/A(e,e'p): 508633a_[corrected_]eep_{d,fe56}_rad_n20k_s4242 (P8.2); pion: 3b27d5c_[corrected_]pion_{h_k1_rad_dec1,d_k1_pim_rad_dec1}_n20k_s4242 (P9.2); kaon: 1275770_[corrected_]kaon_{h_k1_rad_dec1,d_k1_sigm_rad_dec1}_n20k_s4242 (P10.2). Since P10.4 (D-044) the tests use the same set regenerated on glibc 2.43 / gcc 16.2.1: 4408493_glibc243_<same suffix> (validation/scripts/regen_baselines.sh; README.txt records the toolchain, tests/toolchain_guard.sh prints it on a failure); Delta (P11.1): 074cced_glibc243_[corrected_]delta_h_d1_rad_n20k_s4242; SIDIS (P12.2): 31b4232_glibc243_[corrected_]sidis_h_k1_pip_rad_n20k_s4242; rho (P13.2): e6feaf0_glibc243_[corrected_]rho_h_k1_rad_n1k_s4242 |
| `fortran_diff` | one test per corrected Fortran bug, both modes (`tests/fortran_diff`), and the D-032 arm-code checks | reference/functions, reference/runs |
| `examples` | example plug-ins run through the full driver (`ideal_box`, the developer-manual plug-ins `examples/dev_manual`); `package`: install + out-of-tree plug-in build (P15.3) | — |
| `regression_shipped` | full-run identity (RANLUX, `--fortran-compat`) of `simc++` against every shipped-input / P15.F1 Fortran reference run (`validation/scripts/p15_shipped_regression.sh`; also `slow`) | reference/runs/shipped_*, p15f1_* |
| `threads` | `--threads N`: reproducibility, `--threads 1` = sequential, refusals (P15.4) | — |
| `docs` | user and developer manuals: parameter reference, header coverage, verbatim snippets, Markdown links (P15.D1/D2) | — |
| `slow` | long tests (distributions, the shipped-input suite) | — |

L2 (statistical, D-002): `validation/simc_compare` compares a C++ run with the Fortran A/B runs (.hist pulls,
weighted χ²/KS per tree variable, PDF). Results are in reference/comparisons/.

## 8. Source layout

```
include/simc/core/            parameters, RunSetup, events, Registry, Errors, Units, constants, random sources
include/simc/stages/          the interfaces, StageContext
include/simc/<module>/        plug-in classes and kernel libraries; <module>/kernel/ = private Fortran mirrors
src/<module>/                 one CMake target per module (simcfw_<module>; OBJECT if it registers plug-ins)
src/setup/                    simcfw_setup: SetupBuilder.cpp (BuildRunSetup sequence, dbase.f STOP checks, radc_init,
                              luminosity), SetupReaction.cpp (dbase.f:124-547), SetupTheory.cpp (target_init, theory /
                              momentum distribution / spectral-function choice and the file readers), SetupLimits.cpp
                              (limits_init); private SetupInternal.h (P15.3)
src/driver/                   simcfw_driver (+ simcfw_driver_sinks): SimcRun (+ SimcRunParallel.cpp, --threads N, P15.4),
                              MonteCarlo, Pipeline, CommandLine (the simc++
                              command line, CommandLineMain), TargetFieldAngles, report: RunReport.cpp (normalisation,
                              summary), RunReportHist.cpp (.hist), RunReportGeni.cpp (.geni), RunReportRecords.cpp
                              (reaction lines, FORMAT records; private RunReportInternal.h)
tools/                        simc_main.cpp (simc++ = CommandLineMain), inp2json, hist2json, check_all.sh
cmake/                        SIMCConfig.cmake.in (the installed CMake package)
examples/                     example plug-ins: ideal_box (in-tree), out_of_tree_plugin (find_package(SIMC)); tests
tests/{unit,l1,regression,fortran_diff}/   CTest suites (labels above)
validation/                   simc_compare (own CMake project): simc_compare.cpp (CLI, loading, exit status),
                              CompareCommon.{h,cpp} (config, tree, records, session), CompareRunLevel.cpp,
                              CompareTreeLevel.cpp, CompareReport.cpp (verdict, summary.txt/json, PDF); scripts/
data/                         byte copies of the Fortran data files with MD5SUMS: matrices/simc_fortran (COSY maps),
                              spectral, pion, kaon, sidis, field; configs (data/config/*.json)
```

**Installation and package (P15.3).** `cmake --install <build> --prefix <p>` installs `simc++`, `inp2json`,
`hist2json`, `simc_root_identity` (bin/), the public headers (include/simc, without the private kernel/ headers), the
libraries (lib/simc: the static libraries and, for the OBJECT libraries that carry SIMC_REGISTER registrations, their
object files in lib/simc/objects), the data (share/simc/data), docs and examples, and the CMake package
lib/cmake/SIMC (SIMCConfig.cmake, SIMCConfigVersion.cmake, SIMCTargets.cmake). `find_package(SIMC 2.0)` gives the
target `SIMC::simc` (everything a complete run links; dependencies ROOT 6.40 and nlohmann_json are found by the
config). The default data directory (P15.4, D-069; `include/simc/core/DataPath.h`, `simc++ --print-data-dir`) is: the
environment variable `SIMC_DATA_DIR`, else for a program built with find_package(SIMC) the package's share/simc/data
(registered by the interface source share/simc/cmake/SimcInstalledDataRoot.cpp that SIMCConfig.cmake adds to
`SIMC::simc`; relocatable), else for an installed program `<exe>/../share/simc/data`, else the source tree (in-tree
builds and tests). `--data-dir` and the stages' `matrix_dir`/`data_dir` JSON keys still override it. CTest
`package_out_of_tree_plugin` installs into a temporary prefix, builds examples/out_of_tree_plugin with
find_package(SIMC), runs it, and checks that the installed `simc++` and `my_simc` use the installed data without
`--data-dir` (and that `SIMC_DATA_DIR` overrides it).

**Performance and the parallel event loop (P15.4).** docs/manual/markdown/reference/PERFORMANCE.md: the profile, the bitwise-neutral
single-thread optimisations (tgamma memo, COSY power tables) and `--threads N` (`SimcRun::runThreaded`,
src/driver/SimcRunParallel.cpp): N worker `SimcRun`s with their own stages and mt64 streams, merged in the fixed worker
order (`MergeTotals`, `IParallelSink` sinks), refused with `--fortran-compat` / RANLUX.

## 9. Known limitations of the design (to be lifted with the ports that need them)

* (lifted in P8.2, D-030) the event weight is assembled by the reaction (`IReaction::completeMain`); `reactions::assembleWeight`
  serves the (e,e'p) family only (every other reaction assembles its own weight).
  The pion reaction (P9.2) assembles its weight itself (`PionProduction::completeMain`: peepi, Delta factors,
  tgtweight Z/N, Coulomb factor), the kaon reaction (P10.2) too (`KaonProduction::completeMain`: peeK, tgtweight,
  Coulomb factor, the survival probability for decay-OFF runs, which the "kaon" model leaves in the COMMON decdist,
  `main.ntup.decdist`, as physics_kaon.f:164), the Delta reaction (P11.1) too (`DeltaProduction::completeMain`: peedelta,
  tgtweight 1, Coulomb factor).
* COMMON /gnrl/ (Mh, Mh2, Mh2_final) is a set-up constant except for doing_rho, where generate_rho (the Breit-Wigner rho
  mass) and rho_decay (the decay pion) change it per try: `MainEvent::gnrl` with `gnrl_per_event` (set by `rho_production`)
  is then used by the radiator (generate_rad), the montecarlo (spectrometer mass, recon energy, energy-loss correction) and
  complete_recon_ev (P13.2).
* Two optional `IReaction` hooks exist for Fortran COMMON state that the driver otherwise would not see (P9.2):
  `prepareCentral(main0, vertex0)` (calculate_central: the pion reaction fills main0 epsilon/t/phi_pq in corrected
  mode, D-038 F5) and `finalTargetState(targ)` (the report prints targ%Mrec_struck as the run leaves it, the last
  Breit-Wigner Delta mass, simc.f:844). COMMON /pfermi/ is `MainEvent::pfermi`, set by the reaction's generate (and
  by calculate_central for the central event).
* The spectral-function file of D/A(e,e'p) is chosen before the set-up (the set-up reads the theory file, which also
  fixes the generation limits): `SimcRun` reads `pipeline.reaction.spectral_function` and passes it to
  `BuildRunSetup` (`SetupOptions::sf_file`, `use_benhar_sf`).
* The arm codes `electron_arm`/`hadron_arm` (one table, `core/SpectrometerId.h`, D-032) still drive Fortran logic outside the spectrometer: the window
  materials in `hallc_cryo`, the `.geni` blocks and the default plug-in choice.
* The set-up (`BuildRunSetup`) is a port of dbase.f/init.f with reaction-specific branches. Since P15.3 it is split
  into phases by responsibility (§8), not per stage plug-in.
* Plug-ins are linked, not loaded at run time (SIMC_REGISTER at static initialisation): an out-of-tree plug-in is built
  into its own executable against the installed package (EXTENDING.md §1, examples/out_of_tree_plugin).
