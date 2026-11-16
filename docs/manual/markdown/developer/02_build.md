# 2. Build system

## 2.1 Requirements and a build

C++20 compiler (gcc 16.2.1 on the reference machine), CMake ≥ 3.15 (≥ 3.29 for the toolchain hint of the L1 tests),
ROOT ≥ 6.40 (components Core, Tree, Hist, RIO), nlohmann/json ≥ 3.9 (fetched with FetchContent v3.11.2 when not
installed), Threads, Python 3 (only for the `docs` tests), `xz` (for the compressed captures of the replay tests).

The manuals use this layout (the same as the user manual's installation chapter):

| variable | example | content |
|---|---|---|
| `$SIMCPP_SRC` | `~/src/simcpp` | the SIMC++ source tree |
| `$SIMCPP_BUILD` | `~/build/simcpp` | a build directory |
| `$SIMCPP_PREFIX` | `~/opt/simcpp` | an installation (`bin/`, `lib/`, `include/`, `share/simc/data`) |
| `$SIMC_FORTRAN` | `~/src/simc_gfortran` | the Fortran SIMC tree (optional, validation only: its `infiles/` are test inputs) |
| `$SIMC_REFERENCE` | `~/simc_reference` | the Fortran reference data (optional, validation only) |
| `$SIMC_BASELINES` | `~/simc_baselines` | the C++ regression baselines (optional, validation only) |

```bash
export LD_LIBRARY_PATH=$ROOTSYS/lib          # the ROOT 6.40 install; no other ROOT or libtbb in the path
cmake -S $SIMCPP_SRC -B $SIMCPP_BUILD -DCMAKE_PREFIX_PATH=$ROOTSYS -DCMAKE_BUILD_TYPE=Release
cmake --build $SIMCPP_BUILD -j8
ctest --test-dir $SIMCPP_BUILD -j8 --output-on-failure -LE slow    # everything except the long tests
```

Without the validation data every test that needs it is skipped (exit 77, "Skipped"); with them, add
`-DSIMC_REFERENCE_DIR=$SIMC_REFERENCE -DSIMC_FORTRAN_DIR=$SIMC_FORTRAN -DSIMC_BASELINE_DIR=$SIMC_BASELINES` (or export
`SIMC_REFERENCE`, `SIMC_FORTRAN`, `SIMC_BASELINES` before the first configure: they are the defaults).

Flags: `-Wall -Wextra -pedantic` everywhere (the build has 0 warnings; keep it so), Release `-O3 -DNDEBUG`. **Do not
add `-march=native`, `-ffast-math` or `-ffp-contract=fast`**: FMA contraction and value-unsafe optimisations change the
floating-point results and break every bitwise reference test (the numerical rules are in the
[translation report](../../../translation/markdown/04_numerics.md)).

The comparison tool `simc_compare` is a separate CMake project in `validation/` (its own build directory and CTest
suite; the statistics are documented in [validation/README.md](../../../../validation/README.md)).

## 2.2 Layout and targets

The top-level `CMakeLists.txt` defines the ROOT-free core library and the tools, then adds one `CMakeLists.txt` per
module (so that parallel tasks never edit the same file):

| target | type | sources | content |
|---|---|---|---|
| `simcfw` | STATIC | src/core, src/input, src/output/{NtupleLayout,SummaryJson} | parameters, RunSetup, SpectrometerId, DataPath, `.inp` reader, ntuple layout, summary JSON (ROOT-free) |
| `simcfw_random` | STATIC | src/core/random | `Mt64Random` ("mt64"), `RanluxRandom` ("ranlux"), `ScriptedRandom`, `Gauss1` |
| `simcfw_root` | OBJECT | src/output/{RootTreeSink,RunSummarySink,GenHistSink} | sinks `root_tree`, `run_summary`, `gen_hist` |
| `simcfw_physics` | OBJECT | src/physics/{proton,nuclear} | models `sigep`, `deforest*`; spectral functions `theory`, `benhar` |
| `simcfw_material` | OBJECT | src/material | target `hallc_cryo`, material `simc_material` |
| `simcfw_kinematics`, `simcfw_setup` | STATIC | src/setup | `loren`, angles; `BuildRunSetup` (run set-up) |
| `simcfw_optics` | STATIC | src/optics | COSY engine, forward chains, `lfit`, target field |
| `simcfw_hms`, `_shms`, `_sos`, `_hrs`, `_calo` | OBJECT | src/spectrometers/* | the spectrometer plug-ins and kernels |
| `simcfw_pion`, `_kaon`, `_delta`, `_sidis`, `_rho` | OBJECT | src/physics/* | cross sections and their `IPhysicsModel` plug-ins |
| `simcfw_radiation` | OBJECT | src/radiation | `simc_rad` and the radiation kernel |
| `simcfw_reactions`, `simcfw_reactions_{pion,kaon,delta,sidis,rho}` | OBJECT | src/reactions | beam `standard`, the reactions, `EventChain` kernels |
| `simcfw_driver` | STATIC | src/driver | `SimcRun`, `MonteCarlo`, report, pipeline, command line |
| `simcfw_driver_sinks` | OBJECT | src/driver/ReportSink.cpp | sink `fortran_hist` |
| `simcpp` (output name `simc++`) | executable | tools/simc_main.cpp | `return simc::driver::CommandLineMain(argc, argv);` |
| `inp2json`, `hist2json` | executables | tools/ | `.inp` → JSON; `.hist` text report → summary JSON |
| `simc_root_identity` | executable | tests/regression/root_identity.cpp | bitwise comparison of trees/histograms (ignores ROOT dates) |
| `simc` / `SIMC::simc` | INTERFACE | — | everything a complete run links (the installed package's target) |

`SIMC_DRIVER_LIBS` (a cache variable set in src/driver/CMakeLists.txt) lists every library of a complete run. A test or
a program that runs the driver links `${SIMC_DRIVER_LIBS}`.

## 2.3 OBJECT libraries and registration

Plug-ins register themselves by name during static initialisation (`SIMC_REGISTER`, core/Registry.h, chapter 3). The
linker drops an object file of a static archive that nothing references, so **every library that contains a
`SIMC_REGISTER` line is an OBJECT library** and is linked directly by the executable. The exception is
`simcfw_random` (STATIC): `SimcRun::initialise` constructs `Mt64Random` / `RanluxRandom` directly when the registry
does not contain them, which pulls their object files (and so their registrations) in (SimcRun.cpp l.127-139).

The objects of an OBJECT library reach only targets that link it *directly*. The INTERFACE target `simc` therefore
passes them explicitly as `INTERFACE_SOURCES` (`$<TARGET_OBJECTS:...>`), in the build tree and in the installed package
(top-level CMakeLists.txt, "Installation"). This is what lets an out-of-tree executable that links `SIMC::simc` get
every built-in stage without whole-archive flags.

A new module: add `src/<module>/CMakeLists.txt` with an OBJECT library if it registers plug-ins (STATIC otherwise),
`add_subdirectory` it in the top-level file, and append it to `SIMC_DRIVER_LIBS` (chapter 4).

## 2.4 Options and cache variables

| variable | default | meaning |
|---|---|---|
| `SIMC_BUILD_EXAMPLES` | ON | build `examples/` and link the example plug-ins (`ideal_box`) into `simc++` |
| `SIMC_REFERENCE_DIR` | env `SIMC_REFERENCE`, else empty | reference data of the reference tests; empty or absent: those tests are skipped (exit 77) |
| `SIMC_FORTRAN_DIR` | env `SIMC_FORTRAN`, else empty | directory whose `infiles/` provide the shipped test inputs (read-only); empty: skipped |
| `SIMC_BASELINE_DIR` | env `SIMC_BASELINES`, else empty | C++ regression baselines; empty: skipped |
| `SIMC_*_BASELINE_COMMIT` | `55cf7a0_glibc243` | the baseline folder prefix per family (tests/regression/CMakeLists.txt) |

`src/core/DataPath.cpp` is compiled with `SIMC_SOURCE_DATA_ROOT=<source>/data`, the last fallback of the data directory
(chapter 3 §3.8).

## 2.5 Tests and CTest labels

Tests are grouped by directory and labelled:

| label | where | what | count at P15.D1 |
|---|---|---|---|
| `unit` | tests/unit, tests/fortran_diff | parameters, `.inp` reader, sinks, summary/hist2json, stage headers and registry, `lint_no_pow_square` | 75 |
| `l1` | tests/l1 | routine-level reference tables (bitwise or a stated tolerance), replays, full-run identity tests | 282 |
| `replay` | tests/l1 | subset of `l1`: the complete MonteCarlo / SimcRun driven by a recorded random stream (`ScriptedRandom`; names with `_replay_`) | 53 |
| `cli` | tests/l1, tests/fortran_diff, examples | the `simc++` command line (.inp vs JSON, refusals, `--pipeline`) | 10 |
| `fortran_diff` | tests/fortran_diff | one test per behaviour that differs between the default and the `fortran_compat` mode, both modes; input refusals | 47 |
| `regression` | tests/regression | before/after output identity against the C++ baselines (both modes, every family) | 30 |
| `regression_shipped` | tests/regression | full-run identity of `simc++` (`ranlux`, `--fortran-compat`) against the reference run of every shipped input (also `slow`) | 41 |
| `threads` | tests/regression | `--threads N` reproducibility, `--threads 1` = sequential, refusals | 4 |
| `examples` | examples | the example plug-ins through the complete driver | 5 |
| `package` | examples | install into a temporary prefix, out-of-tree plug-in build with `find_package(SIMC)`, run | 1 |
| `slow` | several | long tests (distributions, shipped suite) | 53 |
| `docs` | tests/user_manual, tests/developer_manual | the manuals: parameter reference, header coverage, verbatim snippets, Markdown links | 4 |

`ctest -L l1`, `ctest -LE slow`, `ctest -R replay_pion` select subsets. Tests return 77 (`SKIP_RETURN_CODE`) when
their reference input is missing. The L1 and fortran_diff tests run through `tests/toolchain_guard.sh`
(`CMAKE_TEST_LAUNCHER`): on a failure it prints the toolchain of the references (`reference/TOOLCHAIN.txt`) next to
the current one, because a libm update changes bits. What each reference test proves is described in the
[translation report](../../../translation/markdown/02_validation.md).

The `simc_compare` project has its own suite (`cmake -S validation -B vbuild && cmake --build vbuild && ctest
--test-dir vbuild`): A-vs-A, reference-vs-reference, sensitivity (rad vs norad, shifted `hsdelta`), JSON vs `.hist`,
empty trees, NaN entries, the wtcontr subset and mode detection.

## 2.6 `tools/check_all.sh`

The acceptance check of every change: a clean build of everything plus both test suites, logged.

```bash
TMPDIR=<large scratch dir> ROOTSYS=<ROOT install> tools/check_all.sh --commit <hash> --label <task>
```

* Paths: ROOT from `$ROOTSYS` (else `root-config --prefix`); the validation data from the environment variables
  `SIMC_REFERENCE`, `SIMC_FORTRAN`, `SIMC_BASELINES`, or, for those not set, from the untracked file `check_all.env` in
  the project directory (the parent of the source tree; shell syntax, one `NAME=path` per line). Nothing machine-specific
  is committed. Without the data check_all still runs, with those tests skipped; `git.txt` records which data were used.

* `--commit REV` builds a clean detached worktree of REV (use it whenever other agents have uncommitted work).
* It takes the same lock as `run_simc.sh` (`simc_runs/.lock`), so runs are serialised; check `ps` first.
* It writes `SIMC++_Project/simc_runs/<timestamp>_check_<label>/` (summary.md, git.txt with the toolchain
  fingerprint, build.log, tests.log, val_build.log, val_tests.log), removes the build trees, and appends one line to
  `simc_runs/index.tsv`: `run_id overall build run tests seconds config` (`build` and `tests` are
  `<main>/<validation>` exit codes; config `check_all@<commit>`). The exit code is 0 only if both builds and both
  suites pass. It runs every test including `slow` (about one hour with all validation data).
* Use a `TMPDIR` on a large disk (the worktree and the build trees are created there).

## 2.7 The installed package

```bash
cmake --install build --prefix $HOME/simc
```

installs `bin/` (`simc++`, `inp2json`, `hist2json`, `simc_root_identity`), `include/simc` (the public headers,
**without** the `kernel/` directories), `lib/simc` (static libraries) and `lib/simc/objects` (the objects of the
OBJECT libraries), `lib/cmake/SIMC` (`SIMCConfig.cmake`, `SIMCConfigVersion.cmake`, `SIMCTargets.cmake`),
`share/simc/data` (COSY maps, tables, spectral functions, field map), `share/simc/docs`, `share/simc/examples`, and
`share/simc/cmake/SimcInstalledDataRoot.cpp`.

`find_package(SIMC 2.0)` (version compatibility `SameMajorVersion`) finds ROOT 6.40, nlohmann_json and Threads, and
provides `SIMC::simc` and `SIMC_DATA_DIR`. `SIMCConfig.cmake` adds the interface source `SimcInstalledDataRoot.cpp`
to `SIMC::simc`, compiled into the consumer with `SIMC_INSTALLED_DATA_DIR=<prefix>/share/simc/data`: it registers the
package's data directory at start-up (`simc::detail::SetInstalledDataRoot`), so a consumer program finds the installed
data without `--data-dir`. The CTest `package_out_of_tree_plugin`
(`examples/out_of_tree_plugin/package_check.sh`) proves the whole chain on every `check_all`.
