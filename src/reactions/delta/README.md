# src/reactions/delta: (e,e'p)pi in the Delta region

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `delta_production` | `IReaction` | `reactions::DeltaProduction` | `DeltaProduction.cpp` | H, D, A(e,e'p)pi with the proton detected: generation and completion of `PionProduction`, its own weight (Fortran event.f complete_main, doing_delta branch) |

The header is `include/simc/reactions/DeltaProduction.h`; the weight model is the plug-in `delta` of `src/physics/delta/`.

## Physics notes

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_05_delta.pdf` | (e,e'p)pi in the Delta region: the reaction `delta_production` with the model `delta` | `docs/physics/latex/05_delta.tex` |

These PDFs are **generated copies**: the build target `physics_pdf` writes them from the masters in `docs/physics/pdf/`,
so they are always identical to the masters (test `docs_physics_pdf_copies`). Do not edit or replace them by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
