# SIMC++ user manual

How to install and run SIMC++. No knowledge of the legacy Fortran SIMC is needed. The developer manual
([../developer/index.md](../developer/index.md)) covers maintaining and extending the code.

| chapter | content |
|---|---|
| [1. Introduction](01_introduction.md) | what SIMC++ simulates, run modes, limitations |
| [2. Installation](02_installation.md) | requirements, build, install, data directory (`SIMC_DATA_DIR`, `--print-data-dir`), tests |
| [3. Quick start](03_quick_start.md) | run an example input, the outputs, look at the tree in ROOT, the same run from JSON |
| [4. Inputs](04_inputs.md) | `.inp` syntax, old names, JSON, plug-in selection, **complete parameter reference** |
| [5. Reactions](05_reactions.md) | (e,e'p), pion, kaon, Delta, SIDIS, rho, phase space: what each computes and its limitations |
| [6. Spectrometers](06_spectrometers.md) | HMS, SHMS, SOS, HRS, calorimeter; spectrometer numbers; `spect_mode`; target field |
| [7. Outputs](07_outputs.md) | tree branches per reaction, `.hist`, `.geni`, `.gen`, summary JSON, exit status |
| [8. Running modes and reproducibility](08_running_and_reproducibility.md) | seeds, mt64 vs RANLUX, the legacy-compatible mode, reproducing legacy SIMC results exactly, `--threads` |
| [9. Performance tips](09_performance.md) | speed, threads, memory |
| [10. FAQ and troubleshooting](10_faq.md) | common stops, old inputs, data not found, NaN columns, ... |
| [Example inputs](examples/) | `h_elastic_rad.inp` (chapter 3), `legacy_names.inp` (chapter 4.1.3) |
| [Verified examples](VERIFIED_EXAMPLES.md) | transcript of every command of this manual, as run on the stated commit |

Reference documents: `docs/INP_JSON_MAPPING.md` (the `.inp` → JSON table), `docs/OUTPUT_FORMAT.md` (every output value),
`docs/DRIVER.md`, `docs/PERFORMANCE.md`. Comparison with the legacy Fortran SIMC and the list of differences:
`docs/translation/`.

Conventions: commands are shown as typed in a shell with `$SIMCPP_PREFIX/bin` on the `PATH` (directory layout: [2.1](02_installation.md#21-directory-layout-used-in-this-manual)). `$EX` is the directory of the example inputs (`$SIMCPP_PREFIX/share/simc/docs/manual/user/examples`).
`$SIMC_REFERENCE` (the project's reference runs of the legacy program) is needed only for the example of chapter 8.4.
