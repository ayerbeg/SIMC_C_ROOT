# SIMC++

SIMC++ is a C++20 / ROOT 6 translation of SIMC, the Hall C coincidence Monte Carlo (Fortran, `simc_gfortran`). It reads
the Fortran `.inp` files unchanged, and it writes the Fortran `.hist`/`.geni`/`.gen` reports and the `h10` ROOT tree.
Every stage of the event chain is a plug-in that can be chosen from a JSON configuration. Every reaction, spectrometer and
the polarized-target field of the Fortran are ported. Clear Fortran bugs are corrected by default; `--fortran-compat`
restores the exact Fortran behaviour, and with the RANLUX generator SIMC++ then reproduces Fortran runs file for file.

## Manuals

* **[User manual](docs/manual/user/index.md)**: installation, quick start, inputs (complete parameter reference),
  reactions, spectrometers, outputs, run modes and reproducibility, performance, FAQ.
* **[Developer manual](docs/manual/developer/index.md)**: architecture, Fortran → C++ map, numerical conventions, validation
  framework, extending SIMC++.

Translation and validation (comparison with the Fortran SIMC, every difference): `docs/translation/`, with the user
guides [docs/translation/user_comparing_with_fortran.md](docs/translation/user_comparing_with_fortran.md) and
[docs/translation/user_differences.md](docs/translation/user_differences.md). Reference documents: [docs/INP_JSON_MAPPING.md](docs/INP_JSON_MAPPING.md), [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md),
[docs/FORTRAN_DIFFERENCES.md](docs/FORTRAN_DIFFERENCES.md), [docs/DRIVER.md](docs/DRIVER.md),
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), [docs/EXTENDING.md](docs/EXTENDING.md),
[docs/PERFORMANCE.md](docs/PERFORMANCE.md), [validation/README.md](validation/README.md).

**Delta-region runs (`doing_delta`) have no cross section in the weight**, as in the Fortran (`physics_delta.f:132`):
their yields are not physical rates.

## In short

```bash
export SIMCPP_SRC=~/src/simcpp SIMCPP_BUILD=~/build/simcpp ROOTSYS=~/opt/root LD_LIBRARY_PATH=~/opt/root/lib
cmake -S $SIMCPP_SRC -B $SIMCPP_BUILD -DCMAKE_PREFIX_PATH=$ROOTSYS -DCMAKE_BUILD_TYPE=Release
cmake --build $SIMCPP_BUILD -j8
$SIMCPP_BUILD/bin/simc++ $SIMCPP_SRC/docs/manual/user/examples/h_elastic_rad.inp -o out --ngen 5000
(cd $SIMCPP_BUILD && ctest -L docs)
```

The Fortran SIMC is not needed. Requirements: a C++20 compiler, CMake ≥ 3.15, ROOT ≥ 6.40 built with C++20, nlohmann/json ≥ 3.9 (fetched if absent).
See [chapter 2 of the user manual](docs/manual/user/02_installation.md) for the options, installation, data directory and
test labels.

## Layout

```
include/simc/, src/   the library (core, stages, driver, reactions, physics, optics, spectrometers, output)
tools/                simc++ main, inp2json, hist2json, check_all.sh, check_md_links.py
tests/                unit/ l1/ fortran_diff/ regression/ user_manual/
validation/           simc_compare (statistical comparison with the Fortran) and validation scripts
examples/             example plug-ins
data/                 COSY maps, model tables, spectral functions, field map, example configurations
docs/                 reference documents and the manuals (docs/manual/user, docs/manual/developer)
```
