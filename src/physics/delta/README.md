# src/physics/delta: weight model of the Delta-region reaction

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `delta` | `IPhysicsModel` | `delta::DeltaModel` | `DeltaModel.cpp` | weight of `delta_production`: photon flux times the centre-of-mass jacobian, no cross section (Fortran physics_delta.f peedelta) |
| (kernel) | -- | `delta::peedelta` | `DeltaModel.cpp` | the formula; also writes the Blok fit (sigcm) |

The header is `include/simc/physics/delta/DeltaModel.h`.

## Physics notes

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_05_delta.pdf` | (e,e'p)pi in the Delta region: the reaction `delta_production` with the model `delta` | `docs/physics/latex/05_delta.tex` |

These PDFs are **generated copies**: the build target `physics_pdf` writes them from the masters in `docs/physics/pdf/`,
so they are always identical to the masters (test `docs_physics_pdf_copies`). Do not edit or replace them by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
