# SIMC Fortran → SIMC++ translation report

SIMC++ is a C++20 / ROOT 6 translation of the Hall C Monte Carlo SIMC (Fortran, `simc_gfortran`). This report collects
everything about the relation between the two programs: how the translation was organised and proved, which C++ code
replaces which Fortran routine, the numerical rules that make the results bit-for-bit identical, the switch that
restores the exact Fortran behaviour, and the catalogue of every Fortran bug and quirk found on the way.

The SIMC++ framework itself (architecture, interfaces, extending) is described in the
[developer manual](../../manual/markdown/developer/index.md); running it in the [user manual](../../manual/markdown/user/index.md).

## Chapters

1. [Methodology](01_methodology.md): the Fortran SIMC in one page, how SIMC++ maps it, the design principles taken from
   the decisions (D-001, D-002, D-007, D-026, D-031, D-032, D-044), how each subsystem was translated, the phases.
2. [Validation framework](02_validation.md): the Fortran L1 harness and its tables, the C++ L1 runners, captures and
   replays, RANLUX full-run identity, the shipped-input suite, `simc_compare` (L2) and the wtcontr sigma (D-053),
   regression baselines, how to add a new Fortran reference end to end.
3. [Fortran → C++ routine map](03_routine_map.md): conventions and the curated map by subsystem, COMMON blocks; the
   generated table of all 172 routines is [fortran_routine_map.md](fortran_routine_map.md).
4. [Numerical conventions for bitwise fidelity](04_numerics.md): units, double / real*4 / REAL(16), powers, operand
   order and NaN signs, Fortran intrinsics, libm and toolchain dependence, random-number order.
5. [The Fortran compatibility layer](05_compat.md): classes and markers, how `--fortran-compat` reaches the code, how to
   add a correction.
6. [Catalogue of Fortran bugs, quirks and deviations](06_catalogue.md) (generated from
   [FORTRAN_DIFFERENCES.md](FORTRAN_DIFFERENCES.md)): every item with its id, Fortran file:line, class, effect and
   what SIMC++ does in each mode; the expert list; the cross-check with the decisions.

User-oriented chapters (maintained with the user manual):

* [Comparing SIMC++ with the Fortran](user_comparing_with_fortran.md)
* [Differences from the Fortran, for users](user_differences.md)

## Sources and tools

Paths: `$SIMC_FORTRAN` and `$SIMC_REFERENCE` are the Fortran tree and the Fortran reference data (e.g.
`~/src/simc_gfortran`, `~/simc_reference`; developer manual, chapter 2); `reference/...` in this report is relative to
`$SIMC_REFERENCE`; `run-archive:<label>` is a run folder of the project's run archive (`simc_runs/<timestamp>_<label>/`,
written by `tools/check_all.sh` and the run script); `project_control/` and `fortran_harness/` are directories of the
project next to the source tree and the reference data.

| what | where |
|---|---|
| Fortran source (read-only) | `$SIMC_FORTRAN` (e.g. `~/src/simc_gfortran`) |
| Fortran references (tables, captures, runs, findings reports, toolchain) | `$SIMC_REFERENCE` (e.g. `~/simc_reference`: `functions/README.md`, `runs/`, `P*.F*_report.md`, `TOOLCHAIN.txt`) |
| Fortran L1 harness | `fortran_harness/` of the project (README.md) |
| decisions, plan, evidence | `project_control/{DECISIONS.md, PLAN.md, VALIDATION_LEDGER.md}` of the project |
| the differences, per item and per reaction | [FORTRAN_DIFFERENCES.md](FORTRAN_DIFFERENCES.md) |
| generators | `docs/translation/tools/gen_fortran_map.py` (routine map), `gen_catalogue.py` (catalogue) |

Re-run both generators after a port or a new finding:

```bash
python3 -I docs/translation/tools/gen_fortran_map.py
python3 -I docs/translation/tools/gen_catalogue.py
```

The CTest `docs_translation` checks that no Fortran routine is left without a C++ citation or a reason, that the
catalogue is up to date with FORTRAN_DIFFERENCES.md, and that the verbatim code quotes match the code.
