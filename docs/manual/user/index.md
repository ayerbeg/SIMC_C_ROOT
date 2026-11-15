# SIMC++ user manual

For physicists who want to run SIMC++, with or without SIMC experience. The developer manual is in
[../developer/index.md](../developer/index.md) (maintaining, extending and re-validating the code).

| chapter | content |
|---|---|
| [1. Introduction](01_introduction.md) | what SIMC++ is, relation to the Fortran SIMC, corrected vs `--fortran-compat` mode, what "validated" means, limitations |
| [2. Installation](02_installation.md) | requirements, build, install, data directory (`SIMC_DATA_DIR`, `--print-data-dir`), tests and their labels |
| [3. Quick start](03_quick_start.md) | run a shipped input, the outputs, look at the tree in ROOT, the same run from JSON |
| [4. Inputs](04_inputs.md) | `.inp` (CTP) syntax, legacy names, JSON, plug-in selection, **complete parameter reference** |
| [5. Reactions](05_reactions.md) | (e,e'p), pion, kaon, Delta, SIDIS, rho, phase space: what each computes and its limitations |
| [6. Spectrometers](06_spectrometers.md) | HMS, SHMS, SOS, HRS, calorimeter; spectrometer numbers; `spect_mode`; target field |
| [7. Outputs](07_outputs.md) | tree branches per reaction, `.hist`, `.geni`, `.gen`, summary JSON, exit status |
| [8. Running modes and reproducibility](08_running_and_reproducibility.md) | seeds, mt64 vs RANLUX, `--fortran-compat`, reproducing a Fortran run exactly, `--threads` |
| [9. Comparing with Fortran SIMC](09_comparing_with_fortran.md) | `simc_compare`, reading its verdict, expected differences |
| [10. Differences from the Fortran](10_differences_from_fortran.md) | every correction of the default mode and what it changes; the expert list |
| [11. Performance tips](11_performance.md) | speed, threads, memory |
| [12. FAQ and troubleshooting](12_faq.md) | common STOPs, legacy inputs, data not found, NaN columns, ... |
| [Verified examples](VERIFIED_EXAMPLES.md) | transcript of every command of this manual, as run on the stated commit |

Detailed references that the manual points to: `docs/INP_JSON_MAPPING.md` (the `.inp` → JSON table),
`docs/OUTPUT_FORMAT.md` (every output value with its Fortran source), `docs/FORTRAN_DIFFERENCES.md` (every difference
with the Fortran, with tests and measured effects), `docs/DRIVER.md`, `docs/PERFORMANCE.md`, `validation/README.md`
(`simc_compare` statistics).

Conventions: commands are shown as typed in a shell with `$PREFIX/bin` on the `PATH` (in a build tree use
`build/bin/simc++`). `$SIMC_FORTRAN` is a `simc_gfortran` tree, and `$SIMC_REF` the project's Fortran reference data
(`reference/`), which is needed only for the comparison examples. "Dxxx" are the project decisions
(`project_control/DECISIONS.md`), quoted for traceability.
