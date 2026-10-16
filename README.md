# SIMC++

SIMC++ is a C++20 / ROOT 6 translation of SIMC, the Hall C Monte Carlo (Fortran, `simc_gfortran`), built as a
plug-in framework. The Fortran is the physics specification: the ported routines reproduce it bit for bit (L1 and
replay tests), and complete runs agree statistically with the Fortran (L2). Around these faithful kernels, every
stage of the event chain (target, beam, reaction, radiation, cross section, spectrometers, outputs) is a plug-in
chosen by name from a JSON pipeline.

Current scope: H(e,e'p) elastic, D(e,e'p) and A(e,e'p) (theory-file or Benhar spectral functions; task P8), exclusive
pion electroproduction (π+/π−/π0 on H, D and A, quasifree, coherent and π Δ final states, decay in flight; task P9; the
π0 runs with the calorimeter arm wait for P14), kaon electroproduction (K+Λ/Σ0/Σ− on H, D and A, quasifree and coherent,
decay in flight or survival weight; task P10), (e,e'p)π in the Delta region (`doing_delta`, H, and D/A in the corrected
mode; task P11), with and without radiation, and phase-space (`doing_phsp`) runs, HMS and SHMS arms. Other reactions and the
SOS/HRS/calorimeter arms are refused with a `NotPorted` error (later phases).

**Delta runs: the weight contains no cross section.** As in the Fortran (physics_delta.f:132, "Need new xsec model."),
the `doing_delta` event weight is phase space × virtual-photon flux × the proton CM→lab jacobian; the model cross section
is not applied. Delta yields are therefore not physical rates (docs/FORTRAN_DIFFERENCES.md, D-048 F4).

Documentation (`docs/`):

| file | content |
|---|---|
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | layers, stage interfaces, plug-in wiring, per-event data flow, Fortran correspondence, tests |
| [DRIVER.md](docs/DRIVER.md) | the `simc++` driver: command line, pipeline JSON, outputs, simc.f correspondence |
| [EXTENDING.md](docs/EXTENDING.md) | how to add a spectrometer, reaction, physics model, target or output sink |
| [INP_JSON_MAPPING.md](docs/INP_JSON_MAPPING.md) | Fortran `.inp` variables → JSON parameters |
| [OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md) | `.root` tree, summary JSON/text, `.hist/.geni/.gen` |

## Build

Requirements: CMake ≥ 3.15, a C++20 compiler, ROOT ≥ 6.40 built with C++20 (Core, Tree, Hist, RIO; the default is
`/home/ayerbe/Builds/ROOT6_install`, ROOT 6.40.04), nlohmann/json ≥ 3.9 (fetched automatically if absent). At run
time put only `$ROOTSYS/lib` in `LD_LIBRARY_PATH`: on this machine `/usr/local/lib` holds an old ROOT and an old
`libtbb.so.12` that break linking against ROOT 6.40 (libImt).

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/home/ayerbe/Builds/ROOT6_install -DCMAKE_BUILD_TYPE=Release \
      -DSIMC_REFERENCE_DIR=/home/ayerbe/NewDisk/reference -DSIMC_FORTRAN_DIR=/home/ayerbe/NewDisk/simc_gfortran
cmake --build build -j
cmake --install build --prefix <prefix>     # optional: programs, public headers, data, docs, examples
```

`SIMC_REFERENCE_DIR` (the Fortran reference data) and `SIMC_FORTRAN_DIR` (the Fortran tree, for its `infiles/`)
are only needed by the tests. Tests whose data are missing are skipped. `-DSIMC_BUILD_EXAMPLES=OFF` leaves the
example plug-ins out of `simc++`. The build is warning-free with `-Wall -Wextra -pedantic`.

Programs (in `build/bin/`): `simc++` (the Monte Carlo), `inp2json` (Fortran `.inp` → JSON), `hist2json` (Fortran
`.hist` → summary JSON), `simc_root_identity` (bitwise comparison of two output ROOT files).

## Run

With a Fortran input file, read exactly like `dbase.f` (CTP syntax, `extra_dbase_file`):

```bash
build/bin/simc++ /home/ayerbe/NewDisk/simc_gfortran/infiles/eep_h_elastic_rad_A.inp -o out --seed 4242
```

With a JSON configuration (`parameters` = the `inp2json` output, plus the optional `pipeline` and `run` sections):

```bash
build/bin/inp2json -o run.json my_run.inp          # or start from data/config/h_elastic_rad.json
build/bin/simc++ run.json -o out
```

Useful options: `--ngen N`, `--seed N`, `--run-name NAME`, `--pipeline FILE` (plug-in choice, for example
`examples/ideal_box/hadron_arm.json`), `--print-pipeline`, `--quiet`. The outputs in `out/` are
`<run>.root` (tree `h10`, the Fortran `make_root_tree` layout), `<run>.summary.json/.txt`, `<run>.gen.root` and
the Fortran text reports `<run>.hist`, `<run>.geni`, `<run>.gen`.

In the SIMC++ project, `tools/run_simc.sh [CONFIG] --label <id>` (in `SIMC++_Project/tools/`) builds, runs and
archives a run under `SIMC++_Project/simc_runs/`.

## Compare with the Fortran

The statistical comparison tool is a separate CMake project in `validation/`:

```bash
cmake -S validation -B build-val -DCMAKE_PREFIX_PATH=/home/ayerbe/Builds/ROOT6_install && cmake --build build-val -j
build-val/simc_compare --a /home/ayerbe/NewDisk/reference/runs/eep_h_elastic_rad_A \
                       --b out/eep_h_elastic_rad_A.summary.json --out cmp --config validation/config/compare_default.json
```

It writes `cmp/summary.txt`, `summary.json` and `compare.pdf`: the `.hist` normalisation pulls, the stop
counters, and the weighted χ²/KS test of every tree variable, with a PASS/FAIL verdict (see `validation/README.md`).
The Fortran reference runs are in `/home/ayerbe/NewDisk/reference/runs/<tag>/`.

## Tests

```bash
cd build && ctest -j8                 # everything
ctest -L unit | -L l1 | -L replay | -L cli | -L regression | -L examples
```

* `l1`: routine-level tables against the Fortran harness (`reference/functions/`).
* `replay`: the complete driver fed with captured Fortran random streams reproduces the Fortran outputs bit for bit.
* `regression`: the outputs of the standard runs (seed 4242; H(e,e'p), D/Fe(e,e'p), π+ H, π− D, K+Λ H, K+Σ− D, H(e,e'p)π Delta) must equal the baseline outputs of
  their reference commits (`SIMC++_Project/baselines/<commit>_*`), in both modes: every tree branch and entry bitwise, the text outputs apart from time stamps.
  Every refactor must keep these green.
* `tools/check_all.sh [--commit REV]`: a clean build of everything plus both test suites, logged in
  `SIMC++_Project/simc_runs/`.

## Layout

```
include/simc/   core/ stages/ driver/ input/ output/ setup/ material/ radiation/ reactions/ physics/ optics/
                spectrometers/      (*/kernel/ = private Fortran-mirror headers)
src/            the same modules
examples/       example plug-ins (ideal_box spectrometer) and their tests
tests/          unit/ l1/ regression/
tools/          simc++ main, inp2json, hist2json, check_all.sh
validation/     simc_compare (L2)
data/           COSY maps, spectral functions, pion and kaon model tables (copies of the Fortran files: matrices/simc_fortran,
                spectral, pion, kaon),
                tables, data/config/*.json
```
