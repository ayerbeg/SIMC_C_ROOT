# SIMC++ driver (task P7.2, updated in P7.4)

The driver is the C++ counterpart of the `simc.f` main program (l.1-623): set-up, central event, event loop,
normalisation and outputs. It is built from the validated stage components (P2-P7.1) and does not change their
numerics. Code: `include/simc/driver/`, `src/driver/`, CLI `tools/simc_main.cpp` (binary `simc++`).

## 1. Running it

```
simc++ <file.inp | file.json> [-o OUTDIR] [--seed N] [--ngen N] [--pipeline overrides.json]
       [--run-name NAME] [--data-dir DIR] [--infiles-dir DIR] [--quiet] [--print-pipeline]
```

* A Fortran `.inp` is read exactly like `dbase.f` (`InpReader`, CTP rules, `extra_dbase_file`, float rounding of
  real literals, D-018).
* A JSON file has up to three sections: `parameters` (the `inp2json` output, every `regallvars` variable in `.inp`
  units, D-003), `pipeline` (plug-in choice, below) and `run` (`name`, `output_dir`). Examples:
  `data/config/h_elastic_norad.json`, `h_elastic_rad.json`, `acceptance_phsp_hms_shms.json`.
* `--ngen` and `--seed` override `ngen` and `random_seed` (the `.hist` then prints the overridden values).
  `random_seed = -1` means time() (dbase.f:119). `random_state_file` is ignored (RANLUX state, D-002).
* Outputs in OUTDIR, `<run>` = the input file stem unless `--run-name`/`run.name`:

| file | content | Fortran |
|---|---|---|
| `<run>.root` | tree `h10`, 46 Float_t branches, successes only (only if `Nntu > 0`) | `worksim/<run>.bin` → `make_root_tree` |
| `<run>.summary.json`, `.summary.txt` | run summary (docs/OUTPUT_FORMAT.md §2, schema-valid) | `.hist` + `.geni` |
| `<run>.gen.root` | the diagnostic histograms as TH1D | `.gen` |
| `<run>.hist`, `<run>.geni`, `<run>.gen` | the Fortran text reports, same gfortran layout | `outfiles/` |

`ngen = 0` writes nothing (simc.f:446).

## 2. Design

```
SimcParameters ──BuildRunSetup(target model's ITargetLossProvider)──> RunSetup (const for the stages, init.f/dbase.f)
PipelineConfig::Defaults(setup) ⊕ JSON "pipeline"
SimcRun owns (unique_ptr, created by name through Registry<I>):
   IRandom rng · ITarget target · IRadiator radiator · IBeam beam · IReaction reaction · IPhysicsModel model
   MonteCarlo (owns the ISpectrometer arms) · IOutputSink sinks[]
wiring: StageContext{setup, target, radiator} at init() — interfaces only, no casts
calculate_central → CentralResult {central, vertex0}
loop: EventState ev (ONE object for the run) ──> stages ──> RunTotals
finish: Normalise → RunReport → RunSummary → sinks (endRun; IReportSink::writeReport)
```

* **Stages and plug-ins.** Every stage is chosen by registry name in the `pipeline` JSON. A stage is a string or
  `{"type": name, <config>...}`; config keys are merged over the defaults. `sinks` replaces the list.
  ```json
  "pipeline": {
    "rng": {"type": "mt64", "seed": 12345}, "target": "hallc_cryo", "beam": "standard",
    "reaction": "hydrogen_elastic", "physics": "sigep", "radiator": "simc_rad",
    "spectrometers": {"hadron": "hms", "electron": {"type": "shms"}},
    "sinks": ["root_tree", "run_summary", "gen_hist", "fortran_hist"]
  }
  ```
  `PipelineConfig::Defaults` reproduces the Fortran choices from the flags: `phase_space` if `doing_phsp`
  (D-009: on hydrogen it still runs the elastic kinematics with sigcc = 1), `hydrogen_elastic` for H(e,e'p),
  `deuterium_eep` (D) / `nucleus_eep` (A > 2), `pion_production` (doing_pion, P9.2), `kaon_production` (doing_kaon, P10.2), `delta_production` (doing_delta, P11.1), `sigep` (H) / `deforest` (D, A) / `pion`
  (doing_pion) / `kaon` (doing_kaon) / `delta` (doing_delta; no cross section, D-048 F4), `simc_rad` (radc_init_ev runs for every event, even without radiation), the arms
  from `electron_arm`/`hadron_arm` (`MonteCarlo::DefaultArm`), `root_tree` only if `Nntu > 0`.
  `simc++ --print-pipeline` shows the resolved pipeline.
* **Wiring (P7.4).** `SimcRun::initialise` creates the stages by name, initialises the target and the radiator
  first and passes them to the other stages in the `StageContext` of `init()` (beam, reaction: `ctx.target`,
  `ctx.radiator`; the spectrometers get the context from `MonteCarlo`). What the driver reads back from the radiator
  (central%rad, Egamma_used, rad_proton_this_ev) goes through `IRadiator::currentValues()`, `photonEnergies()` and
  `hadronRadiates()`. The set-up's start-up losses (target_init, extreme_trip_thru_target) come from a dedicated
  instance of the pipeline target through `ITargetLossProvider`. The driver never casts a stage to a concrete class.
* **D/A(e,e'p) (P8.2).** `deuterium_eep` / `nucleus_eep` (NuclearEep) with a spectral-function component
  (`reaction.spectral_function`: `theory` or `benhar`, optional `file`). The choice reaches the set-up before it is built
  (`PipelineConfig::SpectralFunction`, `SetupOptions::sf_file`, `use_benhar_sf`) because the set-up reads the theory file
  (init.f:theory_init) and uses it for the generation limits. Data files: `--data-dir DIR`, default the source tree's
  `data/spectral` (compiled in as `SIMC_SPECTRAL_DATA_DIR`). The physics model is created before the reaction, which
  assembles the event weight with it (`IReaction::completeMain`, D-030).
* **Arms.** Two arms with the same type and the same merged config share one spectrometer object (one Fortran COMMON
  block of counters, e.g. HMS for both arms); otherwise each arm has its own object.
* **Quiet.** `--quiet` adds `"verbose": false` to every stage config (stage messages only).
* **Per-event state.** `EventState {main, vertex, orig, recon, success, pass_cuts, ntried}` is a single object
  passed to every stage on every try. This is deliberate: the Fortran keeps stale values in these structures across
  tries (orig after a failed generate_rad, main.SP/RECON/FP after an early montecarlo return, ...; D-024) and some of
  them reach the outputs. `SimcRun::resetForTry` applies the per-try resets of simc.f:174-183. All other Fortran
  state (COMMON blocks and SAVEd locals) is inside the stage objects (radiator: RadcState; spectrometers: stop
  counters and COSY maps; montecarlo: tmpfact; target: enerloss statics and th_can).
* **Run state.** `RunTotals` (counters, sums, contrib "found" ranges, slop ranges), `Normalisation`
  (simc.f:360-432) and `RunReport` (also the central event and the stop counters). No globals.
* **Reports.** `IReportSink` is an optional extension of `IOutputSink` for sinks that need the whole `RunReport`.
  `FortranHistSink` ("fortran_hist") writes `.hist`, `.geni` and `.gen` with a small gfortran format emulator
  (`FortranWriter.h`: A/F/E/I/L/X/T edits, list-directed records).

## 3. Fortran correspondence

| simc.f | C++ |
|---|---|
| l.53-101 min_max_init, dbase_read, radc_init, targetfac, luminosity | `BuildRunSetup` (P2) in `SimcRun::initialise` |
| l.85 calculate_central (l.1143-1309) | `driver::calculateCentral` (Central.h); central%rad = `IRadiator::currentValues()` |
| l.105-110 sgrnd(random_seed) | the `rng` stage (`mt64`, seeded with random_seed; or `ranlux`, the Fortran RANLUX stream itself, P13.2) |
| l.169-352 loop | `SimcRun::step` (one try) |
| l.203 generate | `reactions::generateEvent(beam, reaction, ...)` |
| l.211 montecarlo | `MonteCarlo::run` |
| l.216 complete_recon_ev, l.224 complete_main | `IReaction::completeRecon`, `reactions::completeMain` |
| l.226-243 pass_cuts, hard_cuts | `driver::PassCuts`, `SimcRun::step` (FORTRAN-BUG D-006 corrected unless fortran_compat) |
| l.250-287 inc() histograms, l.335 results_ntu_write | the sinks' `fill` (`gen_hist`, `fortran_hist`, `root_tree`) |
| l.295-333 ncontribute, sums, limits_update | `SimcRun::step`, `reactions::limitsUpdate` |
| l.341-347 nevent (ngen < 0 counts tries) | `SimcRun::step` |
| l.360-432 normfac, genvol, resolutions | `driver::Normalise` |
| l.462-620 .geni, .gen, .hist (report l.644-1142) | `FormatGeni`, `FortranHistSink::formatGen`, `FormatHist` |

## 4. Validation

* **Replay (bitwise).** `tests/l1/test_l1_simcrun.cpp` runs the complete driver on the four P4.F1 capture inputs
  with a `ScriptedRandom` holding the captured RANLUX stream. It requires the stream to be consumed exactly, the
  per-try bookkeeping, the normalisation, the 46 ntuple values (double, vs `.bin`) and the `h10` tree (Float_t,
  vs the Fortran `.root`) to be bit-identical, and `.hist` (except the time lines), `.geni`, `.gen` to be
  text-identical (CTest `l1_simcrun_replay_{norad,rad,phsp,phsp_tries}`).
* **CLI.** `cli_json_vs_inp_*`: the `.inp` and the JSON form of the same run give identical outputs.
* **L2.** `simc_compare` of 200k-event C++ runs (Mt64Random) against the Fortran A/B runs
  (reference/comparisons/CPP_vs_F/).

## 5. Not supported yet

Reactions other than H/D/A(e,e'p), exclusive pion (P9.2), kaon (P10.2) and Delta-region (e,e'p)pi (P11.1) production and
`doing_phsp` (rho, semi-inclusive: NotPorted at initialise), `random_state_file`, the semi/rho ntuple layouts. (The polarized-target field `using_tgt_field` is ported since P15.2b: optics::TargetField, owned by driver::MonteCarlo.)
(all spectrometer arms are ported: SOS P14.1, HRS P14.2, calorimeter + pi0 P14.3).
