# src/physics/nuclear: nuclear spectral functions and the nuclear elastic model

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `theory` | `ISpectralFunction` | `nuclear::TheorySpectralFunction` | `SpectralFunctions.cpp` | spectral function of the theory files (h2/c12/fe56/au197.theory) for D/A(e,e'p) |
| `benhar` | `ISpectralFunction` | `nuclear::BenharSpectralFunction` | `SpectralFunctions.cpp` | Benhar tables (benharsf_*.dat) with a constant transparency |
| `nuc_elastic` | `IPhysicsModel` | `nuclear::NucElasticModel` | `NuclearElastic.cpp` | A(e,e')A elastic cross section with a charge form factor |
| (kernels) | -- | `nuclear::SfTable`, `TheorySfWeight` | `SpectralFunctionTables.cpp` | table reading and interpolation (Fortran sf_lookup.f, event.f) |
| (kernels) | -- | `nuclear::SigNucElastic`, `NucFormFactor`, `Fgauss` | `NuclearElastic.cpp` | the formulas of Fortran nuclear_elastic.f |

The headers are in `include/simc/physics/nuclear/`. The spectral functions are components of the reactions
`deuterium_eep` and `nucleus_eep` (`src/reactions/NuclearEep.cpp`); the data files are in `data/spectral/`.

## Physics notes

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_02_quasielastic_eep.pdf` | D(e,e'p) and A(e,e'p): the reactions `deuterium_eep`, `nucleus_eep`, the model `deforest` (in `src/physics/proton/`) and the spectral functions `theory` and `benhar` of this folder | `docs/physics/latex/02_quasielastic_eep.tex` |
| `SIMCpp_Physics_03_nuclear_elastic.pdf` | A(e,e')A nuclear elastic scattering: the reaction `nuclear_elastic` with the model `nuc_elastic` | `docs/physics/latex/03_nuclear_elastic.tex` |

These PDFs are **generated copies**: the build target `physics_pdf` writes them from the masters in `docs/physics/pdf/`,
so they are always identical to the masters (test `docs_physics_pdf_copies`). Do not edit or replace them by hand;
change the LaTeX source and rebuild. The common elements are in note 00 and the whole set is bound as the
*SIMC++ Physics Reference* (`docs/physics/pdf/`); see `docs/physics/README.md`.
