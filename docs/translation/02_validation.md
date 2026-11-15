# 2. Validation framework

SIMC++ is validated at four levels against the Fortran, plus C++-against-C++ identity for refactors. Every result is
recorded in `/home/ayerbe/NewDisk/project_control/VALIDATION_LEDGER.md` (date, task, Fortran reference, C++ run,
commit, metrics, tolerance, verdict).

| level | question | Fortran side | C++ side |
|---|---|---|---|
| L1 tables | does one routine give the same numbers on fixed inputs? | harness drivers → `reference/functions/<routine>/` | `tests/l1/test_l1_*.cpp` (label `l1`) |
| replays | does the whole chain give the same numbers on the Fortran's own random stream? | capture programs (`ld --wrap` recorders) → `reference/functions/event_capture*`, `reference/runs/*_cap` | `l1_*_replay_*` (labels `l1;replay`) |
| RANLUX full runs | is a complete run identical to the Fortran run of the same seed? | `reference/runs/<tag>/` (.inp .hist .geni .gen .root) | `l1_simcrun_ranlux_*`, `regression_shipped_*` |
| L2 | is the run statistically the same (the acceptance test, D-002)? | `reference/runs/<tag>_A`, `_B` | `simc_compare` |
| identity | did a refactor change any output? | — | `regression_*`, `p153_golden.sh` |

## 2.1 The Fortran L1 harness and the tables

`/home/ayerbe/NewDisk/fortran_harness` (README.md there; conventions in `reference/functions/README.md`) links the
**unmodified** SIMC objects (archived into `lib/libsimcobjs.a`, built with the SIMC Makefile flags) with one driver per
routine, `<group>/drv_<routine>.f`. A driver reads `input.txt`, fills the arguments and COMMON fields **by column
name**, calls the routine and writes `output.txt`. The harness is run by the fortran-reference agent only; SIMC itself
is never modified or rebuilt.

Table format (`reference/functions/<routine>/`):

* `input.txt`, `output.txt`: whitespace-separated, first line = header with the column names, `#` comment lines; one row
  per point, same order in both files; `output.txt` echoes the input columns and appends the outputs;
* numbers in `%24.16E` / `ES24.16E3` (17 significant digits: text → double is exact); integer controls written as reals;
* structure fields named `<member>_<submember>` (`e_theta` = `vertex%e%theta`); `gen_*` columns document the generator;
* SIMC-internal units (MeV, MeV/c, MeV², rad, cm);
* `README.txt`: routine, file:lines, every COMMON/field read, units, grid, notes; `driver.log`: stdout incl. warnings;
  `check_*`: cross-checks against archived SIMC runs.

Variants: `<routine>_fixedrng/` (random numbers scripted through the columns `u1..uK` (grnd) and `g1..gK` (gauss1),
with `n_unif`, `n_gauss`, `rng_err` outputs; harness stubs `common/rng_stub_*.f` replace grnd/gauss1);
`distributions/`, `rad_distributions/`, `mc_distributions` (real RNG, summary statistics and histograms for routines whose
random structure is checked statistically); `<x>_state/` (the COMMON state a configuration leaves, e.g. `pion_state`).

## 2.2 The C++ L1 runners

`tests/l1/L1Runner.{h,cpp}` (`simcl1` namespace): `Table::load`, `Row` (columns by name), `ColumnSpec {name, rel_tol,
abs_tol}` (pass if |rel| ≤ rel_tol or |abs| ≤ abs_tol), `compare(in, ref, cols, adapter)`, `report`. The runner checks
that the echoed input columns of `output.txt` are bit-identical to `input.txt` (guards against a stale table) and
reports per column: max |rel|, max |abs|, the worst row, the number of bitwise-equal rows and NaN mismatches. Exit 0
pass, 1 fail, 77 table missing (CTest SKIP).

A routine needs an adapter: given one input row, evaluate the C++ and fill the output columns. Fixed-RNG tables feed
`ScriptedRandom` with the row's `u*`/`g*` columns (mode `Scripted` like the stubs, or `Algorithm` = the real gauss1 on the
scripted uniforms) and also compare the number of draws. The default tolerance is bitwise (or 1e-12 rel, D-010); a test
that needs more states why in its `.cmake` file. One `.cmake` per module (tests/l1/{material,radiation,hms,pion,...}.cmake),
included OPTIONAL by tests/l1/CMakeLists.txt, which labels every test `l1` (+ `replay` for names with `_replay_`).
Data-file tests (`l1_*_data_files`) check the md5 of the copies in `data/`.

## 2.3 Captures and replays

A capture is an unmodified SIMC program linked with `ld --wrap` recorders (e.g. `fortran_harness/event/cap_simc`):
the real routines run unchanged and the wrappers log every random number (`cap_rng.txt`: ntried, phase, kind grnd/gauss1,
value, number of grnd inside each gauss1) and the inputs/outputs of selected routines per try (`cap_generate.txt`,
`cap_complete_ev.txt`, `cap_montecarlo.txt`, ...). `montecarlo` cannot be wrapped (defined next to its caller in simc.o);
its outputs are recorded at the next wrapped call. Captures exist for H(e,e'p) (P4.F1, `event_capture/{norad,rad,phsp,
phsp_tries}`), D/A(e,e'p) (P8.F1), the pion collimator/decay calls (P6.F2), pion, kaon, Delta, SIDIS, rho, SOS/HRS, calo
and the target field (`reference/runs/p15f2_poltar_cap`). Files > 5 MB are `.txt.xz`.

Replays (tests/l1/test_l1_driver.cpp, test_l1_simcrun.cpp and the per-reaction test files):

* **per routine** (`l1_replay_event_capture_*`): every captured call of generate / complete_ev / ... re-evaluated
  from its recorded inputs;
* **montecarlo** (`l1_montecarlo_replay_*`): `MonteCarlo::run` on every captured try, with the captured montecarlo
  random numbers;
* **full run** (`l1_simcrun_replay_*`): `SimcRun` with `setRandom(ScriptedRandom)` fed with the whole captured stream,
  `--fortran-compat`: tree, per-try bookkeeping, normalisation, `.hist/.geni/.gen` must equal the Fortran.

## 2.4 RANLUX full-run identity

With the bit-for-bit RANLUX port (`"rng": {"type": "ranlux"}`, seeded from the `.inp` random_seed like simc.f:107) and
`--fortran-compat`, a SIMC++ run reproduces a Fortran run without any capture. `tests/l1/ranlux_fullrun.sh <simc++>
<simc_root_identity> <ref_run_dir> <run> <workdir>` runs it and requires: `.hist` identical except the BEGIN/END time
lines, `.geni` and `.gen` identical, the `h10` tree bitwise identical (`simc_root_identity`); then the corrected mode
must run and its `.geni` must agree with its own summary JSON. CTests `l1_simcrun_ranlux_*`.

The shipped-input suite (P15.2) applies the same to every Fortran reference of the shipped `infiles/` and of P15.F1:
`validation/scripts/p15_shipped_regression.sh` (one CTest `regression_shipped_<folder>` per `reference/runs/shipped_*`
and `p15f1_*` folder, labels `regression_shipped;slow`), plus `tests/fortran_diff/cli_shipped_inputs.sh` for the inputs
the Fortran itself does not run as shipped (legacy dotted names → the same STOP, spect_mode = 2, rho on A = 4, ...).
Policy (D-058): when an L2 comparison in compat mode flags something, the first check is a RANLUX same-seed run;
bitwise identity proves the flag is a fluctuation.

## 2.5 simc_compare (L2)

`validation/` is its own CMake project (`simc_compare`, `simc_mode_diff`); full documentation in
[validation/README.md](../../validation/README.md). Sources: simc_compare.cpp (CLI, loading, exit status),
CompareCommon.{h,cpp} (config, tree, records, normalisation), CompareRunLevel.cpp, CompareTreeLevel.cpp,
CompareReport.cpp (verdict, summary.txt/json, PDF), RunSummary.{h,cpp} (reads summary JSON or `.hist/.geni`), Stats.

```
simc_compare --a <run A> --b <run B> --out DIR [--config CFG.json] [--mode-a compat|corrected] [--mode-b ...]
             [--skip v1,v2] [--nbins N] [--prefer-hist] [--no-pdf] [--expect PASS|FAIL|IDENTICAL]
```

A run is a `.summary.json`, a Fortran `.hist` (+ `.geni`, `.root`), a `.root` or a directory.

* **Run level**: acceptance Ncontribute/Ntried (two-proportion z), Nco_no_rad_prot fraction, every stop counter
  (binomial pull) and a 2×K homogeneity χ² per spectrometer, the weighted yield per tried event (variance of the mean of
  the N per-try weights), `wtcontr`, `AVERAGE.sigcc`, `normfac` (info), and deterministic scalars (sigcc_central,
  charge, targetfac, luminosity, genvol, Ngen: equal within the print precision of the `.hist`, else a gross failure).
* **Tree level**, per common branch: `chi2_abs` (weighted histograms normalised per tried event: shape and
  normalisation), `chi2_ww` (Gagunashvili weighted/weighted), `chi2_uu` (unweighted), `ks_w` (unbinned weighted KS
  with effective sizes). Non-finite entries are counted and excluded; constant branches are compared as constants.
* **Verdict** (D-011): N counted tests, F with p < α = 0.01; PASS iff P(X ≥ F | Binomial(N, α)) ≥ 1e-3 and no gross
  failure (no deterministic mismatch, no p < 1e-5). The tests are correlated, so the count is approximate.
* **wtcontr sigma (D-053)**: `.hist` wtcontr sums only the `pass_cuts` events (simc.f:229-241, 311-313), so its
  relative error is that of the **pass_cuts subset** of the tree, not of the whole tree yield (which underestimated the
  sigma, e.g. SIDIS 0.415 instead of 0.58). The tool rebuilds the window from the printed SPedge and slop%MC; the
  hadron delta maximum depends on the mode: SPedge%e%delta%max for Fortran references and compat runs (D-006), the
  hadron limit for corrected runs. The mode comes from `run.fortran_compat` in the summary JSON (D-054), else from
  `--mode-a/--mode-b`; a `.hist` is compat. A warning is printed when the subset does not reproduce the run's wtcontr
  within 0.2 % (wrong mode).
* Expected outcomes: compat vs Fortran PASS; corrected vs Fortran may FAIL exactly on the quantities a correction
  changes (listed per reaction in FORTRAN_DIFFERENCES.md §7-17). There is no option to skip run-level quantities yet
  (P15.X).

The L2 scripts of each phase are in `validation/scripts/` (`p10_kaon_l2.sh`, `p11_delta_l2.sh`, `p12_sidis_l2.sh`,
`p13_rho_l2.sh`, `p14_*_l2.sh`, `p15_spect_mode_l2.sh`, `p15b_field_l2.sh`, `p154_l2.sh`); their outputs are in
`SIMC++_Project/simc_runs/<run_id>/` and `reference/comparisons/`.

## 2.6 Regression baselines (C++ vs C++)

* `tests/regression/output_identity.sh SIMCPP ROOT_IDENTITY INPUT SEED RUN_NAME BASELINE_DIR WORK_DIR [args]` runs
  `simc++ INPUT --seed SEED --run-name RUN_NAME -o WORK_DIR --quiet [args]` and requires every output to equal the
  baseline: `.root` trees and histograms bitwise (`simc_root_identity`, ROOT dates ignored), `.hist .geni .gen
  .summary.txt .summary.json` text-identical apart from time stamps, the same set of files. Families: H(e,e'p)
  rad/norad/phsp, D/Fe-56 (e,e'p), pion, kaon, Delta, SIDIS, rho, calorimeter, each in compat (`regression_identity_*`)
  and corrected (`regression_corrected_*`) mode; baselines in `SIMC++_Project/baselines/<prefix>_<family>_s4242/`,
  prefix `SIMC_*_BASELINE_COMMIT` (currently `55cf7a0_glibc243`).
* `validation/scripts/regen_baselines.sh [--commit REV] [--tag TAG] [--only SUBSTR]` builds REV in a clean worktree and
  writes **new** folders `<rev>_<tag>_<suffix>/` (never overwrites) with README.txt (command, toolchain line) and MD5SUMS.
  Then update the `*_COMMIT` cache defaults in tests/regression/CMakeLists.txt. Regenerate only when an output change is
  intended (a correction, a new summary key) or after a toolchain update (D-044).
* `validation/scripts/compare_baseline_sets.sh <simc_root_identity> <baselines dir> <new prefix> <map file>` proves that
  a regenerated set differs from the previous one only by allowed added summary keys (map file:
  `validation/scripts/p151_baseline_map.txt`).
* `validation/scripts/p153_golden.sh make|check`: a golden set of 20 inputs (every reaction family, every arm type, the
  target field) in both modes plus one `simc_compare`; `check` compares two sets file by file (326 files). Use it
  before/after any refactor or optimisation that should change nothing (P15.3, P15.4).

## 2.7 Adding a new Fortran reference, end to end

1. **Specify** what is needed and ask the director, who asks the fortran-reference agent: the routine (file:lines), the
   input grid (columns, units, ranges incl. edge cases), the COMMON state to set (or the `.inp` that sets it), the
   output columns, whether random numbers are involved (fixed-RNG table or capture), and for a run: the `.inp`
   parameters, seed, ngen, A and B runs for L2. Never fabricate a reference with the C++ itself.
2. The fortran-reference agent writes the driver and grid in the harness, regenerates byte-reproducibly, publishes
   `reference/functions/<routine>/` (README.txt with every COMMON field read) or `reference/runs/<tag>/`, and a report
   `reference/P<phase>.F<n>_report.md` with its findings.
3. **C++ test**: an adapter in the module's `tests/l1/test_l1_<module>.cpp` with the L1 runner, a CTest in
   `tests/l1/<module>.cmake` with `SKIP_RETURN_CODE 77` and the tolerance; for runs a RANLUX full-run test
   (`ranlux_fullrun.sh`) or a regression_shipped folder; for L2 a script in `validation/scripts/`.
4. **Evidence**: the run folder (`SIMC++_Project/tools/run_simc.sh ... --label <task>` or `tools/check_all.sh`), the metrics, and a line in
   VALIDATION_LEDGER.md (by the director). Findings go to DECISIONS.md and FORTRAN_DIFFERENCES.md (chapter 5); the catalogue of chapter 6 is regenerated from it.
