# src/physics/pion: pion electroproduction cross sections

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `pion` | `IPhysicsModel` | `pion::PionModel` | `PionModel.cpp` | exclusive pion cross section of the reaction `pion_production` (Fortran physics_pion.f peepi) |
| (kernel) | -- | `pion::peepi` | `PionCrossSections.cpp` | centre of mass to laboratory, photon flux, model blend |
| (kernels) | -- | `pion::sigParam2021`, `exclfit` | `PionCrossSections.cpp` | the 2021 fit of P. Bosted |
| (kernels) | -- | `pion::MaidTables`, `sigmaid` | `PionCrossSections.cpp` | MAID2007 tables below W = 2 GeV |
| (kernels) | -- | `pion::sigBlok`, `sigParam04`, `sigParam3000` | `PionCrossSections.cpp` | older fits (sigBlok is used by the Delta model) |
| (class) | -- | `pion::SigMultipole` | `SigMultipole.cpp` | multipole model (not called, kept for the L1 tests) |
| (function) | -- | `pion::pizeroDecay` | `PizeroDecay.cpp` | pi0 -> 2 gamma |

The headers are in `include/simc/physics/pion/`; the data files (MAID tables, Delta mass table) in `data/pion/`.

## Physics notes

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_04_pion.pdf` | Exclusive pion electroproduction: the reaction `pion_production` with the model `pion` | `docs/physics/latex/04_pion.tex` |

These PDFs are **generated copies**: the build target `physics_pdf` writes them from the masters in `docs/physics/pdf/`,
so they are always identical to the masters (test `docs_physics_pdf_copies`). Do not edit or replace them by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
