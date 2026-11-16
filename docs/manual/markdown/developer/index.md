# SIMC++ developer manual

This manual is for people who maintain and extend the SIMC++ framework: its architecture, data flow, interfaces,
build, tests, plug-ins, parallelism and coding standards. Users who only want to run SIMC++ should read the
[user manual](../user/index.md). How SIMC++ was translated from the Fortran SIMC and validated against it, with every
Fortran detail, is in the [translation report](../../../translation/markdown/index.md).

## Chapters

1. [Overview](01_overview.md): what SIMC++ is, the framework in one page, design principles, source layout.
2. [Build system and tests](02_build.md): CMake layout, targets, OBJECT libraries and registration, options, CTest
   labels, `check_all.sh`, the installed package.
3. [Architecture and data flow](03_architecture.md): driver, event structures, kernels and stages, the registry,
   every stage interface with its contract, the order of calls in one try, the normalisation, and the
   [public header reference](03_architecture.md#39-public-header-reference).
4. [Extending SIMC++](04_extending.md): a new reaction, spectrometer, physics model or sink, step by step, with compiled
   examples; state and thread-safety rules.
5. [Parallelism and performance internals](05_parallel.md).
6. [Coding standards, commits, technical debt](06_standards.md).

## Detailed references

| document | content |
|---|---|
| [ARCHITECTURE.md](../reference/ARCHITECTURE.md) | layers, stage table, plug-in wiring, per-event data flow, source layout (short reference) |
| [EXTENDING.md](../reference/EXTENDING.md) | every stage interface with its contract and the in-tree / out-of-tree plug-in mechanics |
| [DRIVER.md](../reference/DRIVER.md) | the `simc++` command line and `SimcRun` options |
| [PERFORMANCE.md](../reference/PERFORMANCE.md) | profile, bitwise-neutral optimisations, `--threads N` |
| [OUTPUT_FORMAT.md](../reference/OUTPUT_FORMAT.md), [INP_JSON_MAPPING.md](../reference/INP_JSON_MAPPING.md) | output files and the .inp → JSON mapping |
| [validation/README.md](../../../../validation/README.md) | `simc_compare`: usage, statistics, verdict |

## Conventions of this manual

* Paths are relative to the SIMC++ source tree `SIMC++_Claude/` unless they are absolute.
* A code block preceded by an HTML comment `<!-- from: path -->` is copied verbatim from that file, which is compiled
  and tested. The CTest `docs_dev_manual` checks that every copy still matches and that every header of `include/simc`
  is described here.
* When the code and this manual disagree, the code is right: please fix the manual.
