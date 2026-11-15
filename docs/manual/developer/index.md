# SIMC++ developer manual

This manual is for people who maintain, extend or re-validate SIMC++, the C++20 / ROOT 6 translation of the Hall C
Monte Carlo SIMC (Fortran, `simc_gfortran`). Users who only want to run SIMC++ should read the user manual in
[`../user/`](../user/index.md) instead.

The chapters are the entry point. The detailed reference documents in `docs/` are linked from the chapters and kept
up to date with them:

| document | content |
|---|---|
| [ARCHITECTURE.md](../../ARCHITECTURE.md) | layers, stage table, plug-in wiring, per-event data flow, source layout (short reference) |
| [EXTENDING.md](../../EXTENDING.md) | every stage interface with its contract and the in-tree / out-of-tree plug-in mechanics |
| [DRIVER.md](../../DRIVER.md) | the `simc++` command line and `SimcRun` options |
| [PERFORMANCE.md](../../PERFORMANCE.md) | profile, bitwise-neutral optimisations, `--threads N` |
| [FORTRAN_DIFFERENCES.md](../../FORTRAN_DIFFERENCES.md) | every known difference from the Fortran (classes a, a-np, b, c, d, dev) with its test |
| [OUTPUT_FORMAT.md](../../OUTPUT_FORMAT.md), [INP_JSON_MAPPING.md](../../INP_JSON_MAPPING.md) | output files and the .inp → JSON mapping |
| [validation/README.md](../../../validation/README.md) | `simc_compare`: usage, statistics, verdict |

## Chapters

1. [Overview and history](01_overview.md): the Fortran SIMC, how SIMC++ maps it, and the design principles.
2. [Build system](02_build.md): CMake layout, targets, OBJECT libraries and registration, options, CTest labels,
   `check_all.sh`, the installed package.
3. [Architecture and data flow](03_architecture.md): driver, event structures, kernels and stages, the registry,
   every stage interface with its contract, the order of calls in one try, the normalisation, and the
   [public header reference](03_architecture.md#9-public-header-reference).
4. [Fortran → C++ map](04_fortran_map.md): conventions and the curated map. The generated table of every
   routine is in [fortran_routine_map.md](fortran_routine_map.md).
5. [Numerical conventions](05_numerics.md): units, float and double, REAL(16), `powi`, NaN signs, operand order, libm,
   random-number order.
6. [Fortran compatibility layer](06_compat.md): how `--fortran-compat` is passed through the code, and how to add a
   correction.
7. [Validation framework](07_validation.md): L1 tables, captures and replays, RANLUX full-run identity,
   `simc_compare`, regression baselines, the shipped-input suite, adding a new Fortran reference.
8. [Extending SIMC++](08_extending.md): a new reaction, spectrometer, physics model or sink, step by step, and the
   thread-safety rules.
9. [Parallelism and performance internals](09_parallel.md).
10. [Coding standards, commits, technical debt](10_standards.md).

## Conventions of this manual

* Paths are relative to the SIMC++ source tree `SIMC++_Claude/` unless they are absolute. The Fortran source is
  `/home/ayerbe/NewDisk/simc_gfortran` (read-only), the Fortran reference data are in `/home/ayerbe/NewDisk/reference`,
  the Fortran L1 harness is in `/home/ayerbe/NewDisk/fortran_harness`, and the project decisions (`D-xxx`) are in
  `/home/ayerbe/NewDisk/project_control/DECISIONS.md`.
* `simc.f:237` means line 237 of the Fortran file. `event.f:complete_ev` names a Fortran routine.
* A code block preceded by an HTML comment `<!-- from: path -->` is copied verbatim from that file, which is
  compiled and tested. The CTest `docs_dev_manual` checks that the copy still matches.
* The manual describes commit `translation` HEAD at the time of task P15.D1 (2026-11-14). When the code and this
  manual disagree, the code is right: please fix the manual.
