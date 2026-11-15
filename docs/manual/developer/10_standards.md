# 10. Coding standards, commits, technical debt

## 10.1 Coding standards

* **C++20**, `-Wall -Wextra -pedantic`, zero warnings. 4-space indentation, lines up to 120 characters, `#ifndef
  SIMC_<PATH>_H` include guards, the file path as the first comment line of every file.
* **Namespaces** `simc` and one sub-namespace per module (`simc::driver`, `simc::optics`, `simc::radiation`, ...).
  Classes `PascalCase`, functions and members `camelCase`, private members with a trailing `_`; Fortran names are kept
  where they identify the Fortran quantity (`main.SP.p.delta`, `targ.Eloss`, `McHmsArgs`).
* **Ownership**: `std::unique_ptr` for stages, value types for events and settings, references for collaborators that
  outlive the user, raw pointers only for optional non-owning collaborators (`StageContext`). RAII for files.
* **No global mutable state** (chapter 8 §8.7). Constants are `constexpr` (`SimcConstants.h`, `Units.h`).
* **Errors**: `FortranStop`, `NotPorted`, `std::invalid_argument` (configuration), `std::logic_error` (misuse); never
  `exit()` in library code.
* **Faithful kernels**: cite every ported routine (`// Fortran: file.f:routine (l.a-b)`), keep the Fortran structure and
  operation order, follow chapter 5 (no `pow(x, 2)`, `powi` for integer powers, quad where the Fortran is quad), mark
  every difference (chapter 6). Do not "modernise" a kernel expression: refactor around it, and prove that nothing
  changed with the golden set and the replay tests.
* **Interfaces** carry physics, not Fortran-isms: return values instead of output arguments, `std::optional` for "written
  or not", 0-based arrays with a comment on the Fortran index.
* **Tests** for every change: an L1 table for a new port, a both-mode test for a correction, a driver run for a new
  plug-in, an identity check for a refactor; the docs tests for new public headers.
* **Units**: internal units in all interfaces; the unit in a `///` comment on every member that has one.

## 10.2 Commits and the working rules

* Branch `translation` of `SIMC++_Claude`. Never force-push, never rewrite history, never discard changes you did not
  make. Other agents may work in the same worktree: `git add <your paths>` only, never `git add -A`.
* Commit message: `[P<phase>.<task>] (<n>) <what>` (n = the commit number within the task), with the attribution lines
  required by the session at the end. Commit early and often.
* At the end of a task: `tools/check_all.sh --commit <hash> --label <task>` (chapter 2 §2.6) and quote its
  `simc_runs/index.tsv` line in the report.
* The Fortran tree `/home/ayerbe/NewDisk/simc_gfortran` and the references in `/home/ayerbe/NewDisk/reference` are
  read-only for C++ work; new references are requested from the fortran-reference agent through the director
  (chapter 7 §7.7).
* Do not delete run folders in `SIMC++_Project/simc_runs/` or baselines; deletion policy D-064 (only reproducible build
  output and git-tracked files).
* Decisions go to `project_control/DECISIONS.md` (append-only, by the director); evidence to `VALIDATION_LEDGER.md`.

## 10.3 Known technical debt and open items

**Open questions for the user / an expert** (the P15.X list in `project_control/PLAN.md`, D-055):

* the **expert list**: every class-b item of [FORTRAN_DIFFERENCES.md §4](../../FORTRAN_DIFFERENCES.md#4-class-b-ambiguous-fortran-behaviour-kept-questions-for-the-user-and-an-expert)
  (kept as Fortran in both modes), e.g. D-006 magnitude in spect_mode runs, D-057 rho kinematics (F5, F7, F8, F10-F12),
  D-061 F3/F6 (raster fry to SOS/HRS recon, HRS aerogel distance), D-062 G5/G8, D-063 (calorimeter windows reuse
  stale values), D-065 P1 (π0 transported as charged), D-067 T4/T5/T6/T13/T15 (target field);
* the **acceptance mode** (D-009/D-028, D-031 flag 1): `doing_phsp` on hydrogen keeps the elastic kinematics and can give
  0 events with radiation; Fortran behaviour kept (`FORTRAN-OPEN`);
* the **Delta weight** (D-048 F4): the doing_delta weight contains no cross section (phase space × flux).

**Tool and code items**:

* `simc_compare` cannot skip run-level quantities (`--skip` acts on tree variables only): a `--skip-run` option was
  proposed (D-070) for comparisons where a known correction changes a run-level number.
* Stale texts found in the code during P15.D1 (messages and comments, no numerical effect): see §10.4.
* The set-up `BuildRunSetup` is a port of dbase.f/init.f with reaction-specific branches; it is split by responsibility
  (P15.3), not per plug-in, so a new reaction that needs new set-up quantities edits src/setup.
* The arm codes still drive Fortran logic outside the spectrometer plug-ins (target-window materials in `hallc_cryo`,
  `.geni` blocks, default plug-in choice).
* Plug-ins are linked, not loaded at run time (D-069): an out-of-tree plug-in needs its own executable.
* `random_state_file` of the `.inp` is ignored (with a NOTE); the `ranlux` rng accepts a 25-integer state in JSON.

**Data-path notes** (D-069/D-070): data files are found through `SIMC_DATA_DIR`, the installed package, `<exe>/../share`,
then the source tree (chapter 3 §3.8). The `data/` copies are byte-identical to the Fortran files (md5 tests); after
moving an installation, set `SIMC_DATA_DIR` or reinstall. `SIMC_REFERENCE_DIR`, `SIMC_FORTRAN_DIR` and
`SIMC_BASELINE_DIR` default to absolute paths of the development machine; on another machine set them at configure time
(the tests that need them are skipped otherwise). 48 untracked old baseline folders (767 MB) in
`SIMC++_Project/baselines` are not referenced by any test (D-069 item 2; the user decides).

## 10.4 Inconsistencies noted while writing this manual (P15.D1)

Texts in the code that no longer match the behaviour (no numerical effect; to be fixed in a code task):

* include/simc/driver/Pipeline.h header comment: "other reactions are not ported yet (NotPorted)" and "rng mt64 (D-002:
  the RANLUX sequence is not reproduced)": every reaction is ported and `ranlux` exists.
* include/simc/driver/SimcRun.h `RunOptions::data_dir` comment names `SIMC_SPECTRAL_DATA_DIR`, replaced by
  `DataDir("spectral")` in P15.4.
* src/driver/SimcRun.cpp l.121-122: the `NotPorted` message for an empty reaction type says "only H/D/A(e,e'p) and
  doing_phsp; P9+".
* src/reactions/HydrogenElastic.cpp `init`: the message "only H(e,e'p) (doing_hyd_elast) and doing_phsp are ported" for
  pion/kaon/... flags (those reactions are other plug-ins now).
* include/simc/reactions/EventChain.h: `assembleWeight` "Throws NotPorted for the reactions whose cross-section branch is
  not ported (pion, kaon, Delta, rho, SIDIS)": those reactions assemble their own weight and never call it.
