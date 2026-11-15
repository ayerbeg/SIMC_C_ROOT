# 2. Installation

## 2.1 Requirements

| item | version | note |
|---|---|---|
| C++ compiler | C++20 (GCC 16.2 and clang 22 are used by the project; both build without warnings with `-Wall -Wextra -pedantic`) | |
| CMake | ≥ 3.15 (≥ 3.29 adds a toolchain hint to failing bitwise tests) | `CMakeLists.txt`, `tests/toolchain_guard.sh` |
| ROOT | ≥ 6.40, built with C++20; components Core, Tree, Hist, RIO (the validation tool also needs Gpad, Graf, MathCore) | the project uses ROOT 6.40.04 (decision D-047) |
| nlohmann/json | ≥ 3.9 | downloaded automatically (FetchContent, v3.11.2) when it is not installed |
| Python 3 | any recent | only for the `docs` tests ([2.5](#25-running-the-tests)); the tests are skipped when Python is missing |

Runtime library path: put only `$ROOTSYS/lib` in `LD_LIBRARY_PATH`. On the development machine `/usr/local/lib` holds an
older ROOT and an older `libtbb.so.12`, and these break ROOT 6.40 (libImt). Keep that directory out of the path:

```bash
export ROOTSYS=/home/ayerbe/Builds/ROOT6_install        # your ROOT 6.40 installation
export LD_LIBRARY_PATH=$ROOTSYS/lib
```

## 2.2 Building

From the source tree (`SIMC++_Claude/`):

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=$ROOTSYS -DCMAKE_BUILD_TYPE=Release \
      -DSIMC_REFERENCE_DIR=/path/to/reference -DSIMC_FORTRAN_DIR=/path/to/simc_gfortran
cmake --build build -j8
```

The programs are in `build/bin/`:

| program | purpose |
|---|---|
| `simc++` | the Monte Carlo |
| `inp2json` | converts a Fortran `.inp` (and its `extra_dbase_file`) into a SIMC++ JSON configuration ([4.2](04_inputs.md#42-the-json-format)) |
| `hist2json` | converts a Fortran `.hist` + `.geni` into the SIMC++ run-summary JSON |
| `simc_root_identity` | compares two output ROOT files entry by entry and branch by branch (bitwise) |

CMake options and cache variables:

| option | default | meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | (none) | use `Release` (`-O3`). Do **not** add `-march=native` or `-ffast-math`: they change the floating-point results, and the bitwise tests then fail (`docs/PERFORMANCE.md` §4) |
| `SIMC_REFERENCE_DIR` | `/home/ayerbe/NewDisk/reference` | the Fortran reference data (`runs/`, `functions/`). Only the tests use it. Tests whose data are missing are skipped |
| `SIMC_FORTRAN_DIR` | `/home/ayerbe/NewDisk/simc_gfortran` | the Fortran tree, for its `infiles/`. Only the tests use it |
| `SIMC_BASELINE_DIR`, `SIMC_BASELINE_COMMIT` | project paths | the regression baselines (`regression` tests) |
| `SIMC_BUILD_EXAMPLES` | `ON` | builds the example plug-ins (`ideal_box`) into `simc++` |

The comparison tool `simc_compare` is a separate CMake project ([chapter 9](09_comparing_with_fortran.md)):

```bash
cmake -S validation -B build-val -DCMAKE_PREFIX_PATH=$ROOTSYS -DCMAKE_BUILD_TYPE=Release
cmake --build build-val -j8
```

## 2.3 Installing

```bash
cmake --install build --prefix $PREFIX
```

This installs:

* `$PREFIX/bin/` with the four programs;
* `$PREFIX/share/simc/data/` with the data files: COSY maps, spectral functions, pion/kaon/SIDIS model tables, the
  target-field map and example configurations. Each data directory has an `MD5SUMS` file;
* `$PREFIX/share/simc/docs/`, `$PREFIX/share/simc/examples/` and `$PREFIX/share/simc/README.md`;
* the headers, libraries and the CMake package `SIMC` (`find_package(SIMC)`), which you need only to write your own
  plug-ins (developer manual, `docs/EXTENDING.md`).

After `export PATH=$PREFIX/bin:$PATH`, the commands of this manual can be typed as `simc++ ...`. In a build tree, use
`build/bin/simc++ ...` instead; both find their data automatically ([2.4](#24-where-the-data-files-are-found)).

## 2.4 Where the data files are found

`simc++` reads its tables (COSY matrices, MAID/Saghai/CTEQ/fDSS tables, spectral functions, target-field map) from one data
directory. It uses the first of these that applies (`include/simc/core/DataPath.h`, decisions D-069/D-070):

1. the environment variable `SIMC_DATA_DIR`, if it is set and not empty;
2. for a program built against the installed CMake package (`find_package(SIMC)`): that package's
   `<prefix>/share/simc/data`;
3. for an installed program (`<prefix>/bin/simc++`): `<prefix>/share/simc/data`, if it exists;
4. the source tree the program was built from (`<source>/data`), for in-tree builds and the tests.

`simc++ --print-data-dir` prints the directory it will use and exits:

```
$ simc++ --print-data-dir
/home/ayerbe/NewDisk/SIMC++_Project/.tmp_P15.D2/install/share/simc/data
$ SIMC_DATA_DIR=/data/simc simc++ --print-data-dir
/data/simc
```

Two more ways to choose data files:

* `--data-dir DIR` sets the directory of the spectral-function and theory files of D/A(e,e'p), and the first place
  where a target-field map file (`tgt_field_file`) is looked for. For the field map, the current directory and
  `data/field` are tried next;
* the JSON keys `"data_dir"` (pion, kaon and SIDIS physics models) and `"matrix_dir"` (spectrometers) of a pipeline
  stage ([4.3](04_inputs.md#43-the-pipeline-section-choosing-plug-ins)).

A missing data file is an error at start-up (exit status 3), not in the middle of a run.

## 2.5 Running the tests

```bash
cd build
ctest -j8                      # everything (about an hour on 8 threads; the slow tests dominate)
ctest -L docs                  # the user-manual checks (seconds)
ctest -LE slow -j8             # everything except the long full-run tests
```

| label | what it tests | count |
|---|---|---|
| `unit` | unit tests of the C++ components (no Fortran data needed) | 75 |
| `l1` | routine-level tables: the C++ against the outputs of the Fortran routines, bitwise (`reference/functions/`) | 282 |
| `replay` | the complete program fed with the random numbers of a Fortran run reproduces it exactly (subset of `l1`) | 53 |
| `fortran_diff` | each documented difference with the Fortran, in both modes ([chapter 10](10_differences_from_fortran.md)) | 47 |
| `cli` | the command line: `.inp` vs JSON give identical outputs, rejected inputs, shipped inputs | 10 |
| `regression` | outputs of the standard runs (seed 4242) are unchanged with respect to the baselines, in both modes; `--threads` reproducibility | 30 |
| `regression_shipped` | full-run identity with the Fortran for every runnable shipped input (RANLUX, `--fortran-compat`) | 41 |
| `slow` | the long tests: `regression_shipped` and 12 long `l1` runs | 53 |
| `threads` | the parallel event loop | 4 |
| `examples`, `package` | the example plug-in and the installed CMake package | 4, 1 |
| `docs` | this manual: every input parameter and command-line option is documented, and every link resolves | 2 |

(Counts from `ctest -N -L <label>` on the commit of this manual; one test can carry several labels.)

The comparison tool has its own suite of 26 self-tests: `cd build-val && ctest`.

In the SIMC++ project, `tools/check_all.sh [--commit REV] [--label NAME]` does a clean build of both projects, runs both
test suites and logs the result under `SIMC++_Project/simc_runs/`. If other people may have uncommitted changes in the
worktree, use `--commit`. Put `TMPDIR` on a large disk.
