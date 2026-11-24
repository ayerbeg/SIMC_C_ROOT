# src/physics/sidis: semi-inclusive (SIDIS) cross-section model

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `sidis` | `IPhysicsModel` | `sidis::SidisModel` | `SidisModel.cpp` | the SIDIS cross section of the reaction `sidis_production` (reads `data/sidis/`) |
| (kernel) | -- | `sidis::peepiX` | `SidisCrossSection.cpp` | the model (Fortran semi_physics.f peepiX) |
| (kernels) | -- | `sidis::slacemcfit`, `sidis::rhadGlobal` | `SidisCrossSection.cpp` | SLAC EMC ratio and hadron-attenuation fit |
| (class) | -- | `sidis::Cteq5Pdf` | `Cteq5Pdf.cpp` | CTEQ5M parton distributions (cteq5/Ctq5Pdf.f) |
| (class) | -- | `sidis::FdssGrid` | `Fdss.cpp` | fDSS kaon fragmentation functions (fdss/fdss.f) |
| (functions) | -- | `sidis::f1f2in21`, `sf`, `resmod` | `F1F2In21.cpp` | inclusive structure functions (F1F2IN21_v1.0.f) |
| (function) | -- | `sidis::semiDilution` | `SemiDilution.cpp` | dilution factor (semi_dilution.f; not called, as in the Fortran) |

The headers are in `include/simc/physics/sidis/`; the reaction plug-in is in `src/reactions/sidis/`.

## Physics note

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_08_sidis.pdf` | semi-inclusive (e,e'h)X production: the reaction `sidis_production` with the model `sidis` | `docs/physics/latex/08_sidis.tex` |

This PDF is a **generated copy**: the build target `physics_pdf` writes it from the master in `docs/physics/pdf/`, so it
is always identical to the master (test `docs_physics_pdf_copies`). Do not edit or replace it by hand; change the LaTeX
source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss) are in note 00 and
the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see `docs/physics/README.md`.
