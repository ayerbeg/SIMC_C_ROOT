# 6. Coding standards, commits, technical debt

## 6.1 Coding standards

* **C++20**, `-Wall -Wextra -pedantic`, zero warnings. 4-space indentation, lines up to 120 characters, `#ifndef
  SIMC_<PATH>_H` include guards, the file path as the first comment line of every file.
* **Namespaces** `simc` and one sub-namespace per module (`simc::driver`, `simc::optics`, `simc::radiation`, ...).
  Classes `PascalCase`, functions and members `camelCase`, private members with a trailing `_`. Physics quantities keep
  their established short names (`main.SP.p.delta`, `targ.Eloss`), so that code and documents use the same words.
* **Ownership**: `std::unique_ptr` for stages, value types for events and settings, references for collaborators that
  outlive the user, raw pointers only for optional non-owning collaborators (`StageContext`). RAII for files.
* **No global mutable state** (chapter 4 §4.7). Constants are `constexpr` (`SimcConstants.h`, `Units.h`).
* **Errors**: `FortranStop`, `NotPorted`, `std::invalid_argument` (configuration), `std::logic_error` (misuse); never
  `exit()` in library code.
* **Kernels**: keep the operation order and the random-number order of a kernel exactly; do not "modernise" a kernel
  expression, refactor around it and prove that nothing changed with the golden set (`validation/scripts/p153_golden.sh`)
  and the replay tests. The numerical rules (no `std::pow(x, 2)`, which the unit test `lint_no_pow_square` enforces;
  `optics::powi` for integer powers; quad precision where required), the origin comments and the markers for
  behaviour that differs between the two modes are described in the [translation report](../../translation/index.md),
  chapters 3-5.
* **Interfaces** carry physics: return values instead of output arguments, `std::optional` for "written or not",
  0-based arrays (with a comment where a kernel keeps another convention).
* **Tests** for every change: a driver run for a new plug-in, a test on fixed inputs for a new kernel, a both-mode test
  for a new `fortran_compat` behaviour, an identity check for a refactor, the docs tests for new public headers.
* **Units**: internal units in all interfaces; the unit in a `///` comment on every member that has one.

## 6.2 Commits and working rules

* Branch `translation` of `SIMC++_Claude`. Never force-push, never rewrite history, never discard changes you did not
  make. Several agents may work in the same worktree: `git add <your paths>` only, never `git add -A`.
* Commit message: `[P<phase>.<task>] (<n>) <what>` (n = the commit number within the task), with the attribution lines
  required by the session at the end. Commit early and often.
* At the end of a task: `tools/check_all.sh --commit <hash> --label <task>` (chapter 2 §2.6) and quote its
  `simc_runs/index.tsv` line in the report.
* Reference data (`$SIMC_REFERENCE`) and the directory of shipped inputs are read-only; new references
  are produced by the reference process described in the translation report.
* Do not delete run folders in `SIMC++_Project/simc_runs/` or baselines; only reproducible build output and git-tracked
  files may be removed.
* Decisions are recorded in `project_control/DECISIONS.md` (append-only), evidence in `VALIDATION_LEDGER.md`.

## 6.3 Known technical debt

* `simc_compare` cannot skip run-level quantities (`--skip` acts on tree variables only): a `--skip-run` option is
  proposed for comparisons where a known mode difference changes a run-level number.
* `BuildRunSetup` has reaction-specific branches; it is split by responsibility (src/setup: SetupReaction, SetupTheory,
  SetupLimits), not per plug-in, so a new reaction that needs new set-up quantities edits src/setup.
* The arm codes still drive logic outside the spectrometer plug-ins (target-window materials in `hallc_cryo`, `.geni`
  blocks, default plug-in choice).
* Plug-ins are linked, not loaded at run time: an out-of-tree plug-in needs its own executable.
* The workers of `--threads N` are built sequentially (0.1-0.3 s each); `tgamma` remains the largest single cost of
  radiative runs (PERFORMANCE.md §4).
* `random_state_file` of the `.inp` is ignored (with a NOTE); the `ranlux` rng accepts a 25-integer state in JSON.
* `SIMC_REFERENCE_DIR`, `SIMC_FORTRAN_DIR` and `SIMC_BASELINE_DIR` default to absolute paths of the development machine;
  on another machine set them at configure time (the tests that need them are skipped otherwise). After moving an
  installation set `SIMC_DATA_DIR` or reinstall (chapter 3 §3.8).
* 48 untracked old baseline folders (767 MB) in `SIMC++_Project/baselines` are not referenced by any test (to be
  decided by the user).

Open physics questions about the original program's behaviour (the expert list, the acceptance mode, the Delta weight)
are listed in the [translation report](../../translation/06_catalogue.md).

## 6.4 Known loose ends

* The stale code texts noted while writing this manual (outdated "not ported" comments and messages, an old data-directory
  name) were corrected; runtime messages now speak of "legacy SIMC behaviour".
* Doxygen documentation of `include/` is partial (core and the stage interfaces); the coverage test `docs_doc_coverage`
  is disabled until it is completed.
* `src/physics/rho/RhoModel.cpp`: `kNotPortedTag` / `isNotPortedProblem` are dead code.
* `simc_compare` cannot exclude run-level quantities (a `--skip-run` option would help when comparing corrected runs).
* `inp2json` writes a `"generator"` string with an old version tag.
