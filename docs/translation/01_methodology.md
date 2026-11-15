# 1. Methodology

This chapter explains what was translated, how the translation was organised, and how it was proved. The later
chapters give the details: the validation framework (chapter 2), the routine map (chapter 3), the numerical rules
(chapter 4), the compatibility switch (chapter 5) and the catalogue of every Fortran bug and quirk (chapter 6).

## 1.1 The Fortran SIMC in one page

SIMC is the Hall C Monte Carlo for coincidence (e,e'X) experiments. The version translated here is
`/home/ayerbe/NewDisk/simc_gfortran` (gfortran build, flags `-O -w -ffixed-line-length-132 -ff2c -fno-automatic
-fdefault-real-8`). Its structure:

* **Main program** `simc.f` (l.1-623):
  1. start-up (l.53-160): `dbase_read` (dbase.f: the CTP `.inp` reader and every derived quantity), `init.f`
     (`min_max_init`, `target_init`, `limits_init`, `radc_init`, `theory_init`), the target factor and luminosity,
     `calculate_central` (the central event and `CENTRAL.sigcc`), the `using_tgt_field` arm angles (l.120-151);
  2. the event loop (l.169-352): `do while (nevent < |ngen|)`. One pass is one **try**:
     `generate` (event.f) → `montecarlo` (simc.f l.1310-1853: both spectrometer arms) → `complete_recon_ev` →
     `complete_main` (the weight) → `pass_cuts` (l.226-238) and `hard_cuts` → histograms (`inc`, l.250-287) →
     `results_ntu_write` (the ntuple) → the accumulators and `limits_update`;
  3. the end (l.354-620): normalisation (`normfac`, `genvol`, resolutions), the `.geni` and `.gen` files, and
     `report` (the `.hist` file).
* **Event generation** `event.f`: `generate` (beam, vertex, reaction variables, `complete_ev`, `generate_rad`),
  `complete_ev`, `complete_recon_ev`, `complete_main`, `limits_update`, `physics_angles`, `spectrometer_angles`.
  Every reaction is a set of `if (doing_pion) ...` branches inside these routines.
* **Physics**: `physics_proton.f` (sigep, deForest), `physics_pion.f`, `physics_kaon.f`, `physics_delta.f`,
  `semi_physics.f` (+ CTEQ5, fDSS, F1F2IN21), `rho_physics.f`, `generate_rho.f`, `rho_decay.f`, `pizero_decay.f`,
  `sf_lookup.f`, `jacobians.f`.
* **Radiation** `radc.f`, `brem.f`, `init.f:radc_init*`. **Target** `target.f`, `enerloss_new.f`, `shared/musc*.f`.
* **Spectrometers** `hms/`, `shms/`, `sos/`, `hrsl/`, `hrsr/`, `calo/` (one `mc_<arm>` single-arm Monte Carlo each,
  with hut and reconstruction), sharing the COSY transport `shared/transp.f`, `project.f`, `rotate_*axis.f`.
* **Polarized-target field** `trg_track.f`. **Random numbers** `call_ranlux.f` + `cern/ranlux.f` (RANLUX luxury 3),
  `gauss1.f`.
* **State** lives in COMMON blocks and in static locals (`-fno-automatic`: every local is SAVEd):
  `simulate.inc` (`/gnrl/` masses and beam, `/theory/`, `/decd/` decay, `/ntuple_variables/` ntup and decdist,
  `/sigkaon/`, `/pfermi/` Fermi motion), `structures.inc` (the derived types `event`, `event_main`, cuts, limits),
  `target.inc` (`/target_info/ targ`), `radc.inc` (`/radccom/`), `brem.inc`, `sf_lookup.inc` (`/sftable/`),
  `spectrometers.inc` (`/track/`, `/driftdistances/`), `hms/struct_hms.inc` and `shms/struct_shms.inc` (the
  aperture stop counters), `histograms.inc` (the `.gen` histograms), `constants.inc`.
* **Outputs**: `<run>.hist` (report), `<run>.geni` (stop counters), `<run>.gen` (histograms), `<run>.bin` (HBOOK-like
  ntuple, converted to ROOT by `util/root_tree/make_root_tree`).

The machine-readable routine list is [fortran_routine_map.md](fortran_routine_map.md) (172 routines of the files
the Makefile links into `simc`).

## 1.2 How SIMC++ maps it

| Fortran | SIMC++ |
|---|---|
| simc.f main program | `driver::SimcRun` (`initialise` / `step` / `finish`, include/simc/driver/SimcRun.h) |
| simc.f:montecarlo | `driver::MonteCarlo::run` (+ `driver::PizeroCalorimeter` for π0 → 2γ) |
| simc.f:calculate_central | `driver::calculateCentral` (Central.h) |
| simc.f:report, .geni, .gen, normalisation | `driver::RunReport` (`Normalise`, `FormatHist`, `FormatGeni`), sink `fortran_hist` |
| dbase.f + init.f start-up | `BuildRunSetup` → `RunSetup` (include/simc/setup/SetupBuilder.h, core/RunSetup.h) |
| CTP `.inp` reader | `InpReader` → `SimcParameters` (input/InpReader.h, core/SimcParameters.h) |
| event.f:generate (beam part) | `IBeam` plug-in `standard` |
| event.f:generate/complete_ev/complete_recon_ev/complete_main per reaction | one `IReaction` plug-in per reaction (`hydrogen_elastic`, `phase_space`, `deuterium_eep`, `nucleus_eep`, `pion_production`, `kaon_production`, `delta_production`, `sidis_production`, `rho_production`) |
| physics_*.f | `IPhysicsModel` plug-ins (`sigep`, `deforest*`, `pion`, `kaon`, `delta`, `sidis`, `rho`); `ISpectralFunction` (`theory`, `benhar`) |
| radc.f, brem.f | `IRadiator` plug-in `simc_rad` around the `radiation/kernel` |
| target.f, enerloss_new.f, musc*.f | `ITarget` plug-in `hallc_cryo` with `IMaterialModel` `simc_material` |
| mc_hms.f, mc_shms.f, mc_sos.f, mc_hrs[lr].f, mc_calo.f | `ISpectrometer` plug-ins `hms`, `shms`, `sos`, `hrsr`, `hrsl`, `calo` around their kernels |
| transp.f, project.f, recon | `optics::CosyEngine` |
| trg_track.f | `optics::TargetField` (owned by `MonteCarlo`) |
| call_ranlux.f, ranlux.f, gauss1.f | `IRandom`: `mt64` (default), `ranlux` (bit-for-bit Fortran stream), `ScriptedRandom` (tests); `Gauss1` |
| results_write.f + NtupleInit.f + make_root_tree | sink `root_tree` (`RootTreeSink`, `NtupleLayout`) |
| (none) | sinks `run_summary` (`.summary.json/.txt`), `gen_hist` (`.gen.root`); `simc_compare` |

COMMON blocks become members of the object that plays the role of "one Fortran process": `RunSetup` (start-up
constants: `/gnrl/` except the per-try masses, `/target_info/`, cuts, limits), `EventState` (the main-loop locals
`main`, `vertex`, `orig`, `recon`), `MainEvent::ntup` (`/ntuple_variables/`), `MainEvent::pfermi` (`/pfermi/`),
`MainEvent::gnrl` (the per-try `/gnrl/` masses of doing_rho), `radiation::RadcState` (`/radccom/`), the spectrometer
kernels (`/track/`, stop counters), `nuclear::SfTable` (`/sftable/` and the SAVEd locals of sf_lookup),
`optics::TargetField` (the trg_track COMMONs). Nothing is a global variable. The developer manual ([chapter 3](../manual/developer/03_architecture.md))
describes these classes; chapter 3 here gives the routine-level map.

## 1.3 Design principles

* **D-001 / D-007 / D-026: an object-oriented, plug-in framework around faithful kernels.** Spectrometers,
  reactions, targets, physics models and outputs are plug-ins (abstract interface + registration by name), and
  "the design must let any stage of the chain be replaced" (D-007). The Fortran is the specification of the physics:
  each ported routine reproduces the algorithm, the operation order (bitwise results) and the random-number calls
  of the routine it replaces, and cites it (`// Fortran: file.f:routine (l.a-b)`). Fortran-isms (COMMON mirrors,
  1-based arrays, SAVEd locals, stale values) are confined to the **kernel layer** behind the stage interfaces
  (D-026), so they can be refactored without touching users.
* **D-031: clear Fortran bugs are corrected by default; `--fortran-compat` restores the Fortran.** Each correction is
  marked `// FORTRAN-BUG(D-xxx): original Fortran <file:line>: <code>; corrected because ...`. The compat switch exists
  so that the bitwise L1, replay and regression tests can still prove the translation. Ambiguous cases keep the
  Fortran behaviour in both modes (`FORTRAN-DIFF-CANDIDATE`, the "expert list"), the acceptance-mode question is
  `FORTRAN-OPEN(D-031 flag 1)`, faithful oddities are `FORTRAN-QUIRK`. Chapter 5.
* **D-002: validation is statistical.** The acceptance test of a run is L2: `.hist` normalisation pulls plus
  weighted χ²/KS per tree variable against the Fortran (`simc_compare`). L1 (routine tables, bitwise) is the
  debugging and unit-test tool. In practice SIMC++ goes further: with the bit-for-bit RANLUX port (D-058) and
  `--fortran-compat` a complete run is identical to the Fortran run of the same seed, which is proved for every
  shipped input (chapter 2).
* **D-003: the JSON input is a superset of the `.inp`**, with a documented mapping and the converter `inp2json`.
* **D-032: one number per spectrometer** (`core/SpectrometerId.h`; arm code 6 is rejected).
* **D-044: bitwise references are tied to a toolchain** (glibc libm, compilers; `reference/TOOLCHAIN.txt`).
* **D-055: finish the translation first, debug later**: non-blocking items are collected in the P15.X list
  (the expert list, chapter 6).

## 1.4 How the translation was done

* **Roles.** Three kinds of agents worked under a director (`project_control/`): the *fortran-reference* agent built
  and ran the Fortran program and the L1 harness and produced every reference (`/home/ayerbe/NewDisk/reference`), the
  *cpp-translator* agents wrote the C++ against those references, and the director reviewed every task, recorded the
  decisions (`DECISIONS.md`, D-001 ... D-070) and the evidence (`VALIDATION_LEDGER.md`), and ran `tools/check_all.sh` at
  every task close. The Fortran source was never modified; the C++ never wrote a reference for itself.
* **One subsystem per phase**, always in the same order: (F) the Fortran side first — reference tables of every routine
  of the subsystem, captures of real runs, A/B reference runs for L2, and a findings report
  (`reference/P<n>.F<k>_report.md`) listing every bug and oddity found while instrumenting the code; (C) the C++ port,
  routine by routine, each with its L1 test; then the integration into the driver with replays of the captures, a
  regression baseline and an L2 comparison in both modes; (R) the director's review and classification of the findings
  (class a/b/c/d, chapter 5).
* **Reading order of a port.** Read the whole Fortran routine and its callers, the COMMON blocks and includes; note units,
  1-based and column-major indexing, implicit typing (`i`-`n` integers), `real*4` vs `real*8`, `-fno-automatic`
  (static locals), `save`/`data` initialisation and GOTO flow; port the same algorithm with the same operation order and
  the **same number and order of random-number calls**; cite the origin; then prove it.
* **Proof levels** (chapter 2): L1 routine tables (bitwise in practice), replays of captured Fortran random streams
  through the whole chain, complete-run identity with the bit-for-bit RANLUX port, and the statistical L2 comparison
  that D-002 defines as the acceptance test. Every refactor is protected by C++-against-C++ output identity.
* **Fortran bugs** found on the way were reproduced first (so that the bitwise tests prove the port), then — after the
  user's decision D-031 — corrected by default behind the `fortran_compat` switch, each with a test that shows both
  behaviours (chapter 5 and 6).

## 1.5 History (phases)

The project ran in phases (`project_control/PLAN.md`, decisions in `DECISIONS.md`):

| phase | content |
|---|---|
| P0 | audit of the pre-existing C++ (found unfaithful, D-004), Fortran reference build and runs |
| P1 | validation infrastructure: L1 runner, `simc_compare`, summary JSON, `.inp` reader, architecture v1 |
| P2–P3 | core, set-up (dbase.f/init.f), target and energy loss |
| P4–P6 | H(e,e'p) generation, radiation, HMS and SHMS (COSY engine) |
| P7 | main loop, first L2 PASS, modernisation (P7.4: StageContext, kernel/ layers, Errors, Units, examples), corrections (P7.5, D-031) |
| P8–P13 | D/A(e,e'p), exclusive pion, kaon, Delta, SIDIS, rho; RANLUX port (P13.2) |
| P14 | SOS, HRS-L/R, calorimeter arms and π0 → 2γ |
| P15 | hardening: RANLUX settles, shipped-input suite, target-field tracking (P15.2b), clean-up and CMake package (P15.3), performance and `--threads` (P15.4), manuals (P15.D1/D2) |

The legacy pre-translation C++ (snapshot `ab6cd87` on `main`) was removed in P7.2 (`f9cea51`); it is in git history
and is not trusted (D-004).
