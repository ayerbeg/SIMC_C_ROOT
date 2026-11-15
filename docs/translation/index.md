# SIMC Fortran → SIMC++ translation report

SIMC++ is a C++20 / ROOT 6 translation of the Hall C Monte Carlo SIMC (Fortran, `simc_gfortran`). This report collects
everything about the relation between the two programs: how the translation was organised and proved, which C++ code
replaces which Fortran routine, the numerical rules that make the results bit-for-bit identical, the switch that
restores the exact Fortran behaviour, and the catalogue of every Fortran bug and quirk found on the way.

The SIMC++ framework itself (architecture, interfaces, extending) is described in the
[developer manual](../manual/developer/index.md); running it in the [user manual](../manual/user/index.md).

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
   [FORTRAN_DIFFERENCES.md](../FORTRAN_DIFFERENCES.md)): every item with its id, Fortran file:line, class, effect and
   what SIMC++ does in each mode; the expert list; the cross-check with the decisions.

User-oriented chapters (maintained with the user manual):

* [Comparing SIMC++ with the Fortran](user_comparing_with_fortran.md)
* [Differences from the Fortran, for users](user_differences.md)

## Sources and tools

| what | where |
|---|---|
| Fortran source (read-only) | `/home/ayerbe/NewDisk/simc_gfortran` |
| Fortran references (tables, captures, runs, findings reports, toolchain) | `/home/ayerbe/NewDisk/reference` (`functions/README.md`, `runs/`, `P*.F*_report.md`, `TOOLCHAIN.txt`) |
| Fortran L1 harness | `/home/ayerbe/NewDisk/fortran_harness` (README.md) |
| decisions, plan, evidence | `/home/ayerbe/NewDisk/project_control/{DECISIONS.md, PLAN.md, VALIDATION_LEDGER.md}` |
| the differences, per item and per reaction | [FORTRAN_DIFFERENCES.md](../FORTRAN_DIFFERENCES.md) |
| generators | `docs/translation/tools/gen_fortran_map.py` (routine map), `gen_catalogue.py` (catalogue) |

Re-run both generators after a port or a new finding:

```bash
python3 -I docs/translation/tools/gen_fortran_map.py
python3 -I docs/translation/tools/gen_catalogue.py
```

The CTest `docs_translation` checks that no Fortran routine is left without a C++ citation or a reason, that the
catalogue is up to date with FORTRAN_DIFFERENCES.md, and that the verbatim code quotes match the code.
