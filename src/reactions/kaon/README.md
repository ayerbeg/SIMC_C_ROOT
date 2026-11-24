# src/reactions/kaon: the kaon-production reaction

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `kaon_production` | `IReaction` | `KaonProduction` | `KaonProduction.cpp` | N(e,e'K+)Y and coherent hypernuclei; derives from `PionProduction` (generation, vertex completion, central event, in `src/reactions/pion/`) and overrides `completeMain` (kaon cross section, target nucleon count, survival weight) |

The header is `include/simc/reactions/KaonProduction.h`; the cross-section model `kaon` is in `src/physics/kaon/`.

## Physics note

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_06_kaon.pdf` | kaon electroproduction: the reaction `kaon_production` with the model `kaon` | `docs/physics/latex/06_kaon.tex` |

This PDF is a **generated copy**: the build target `physics_pdf` writes it from the master in `docs/physics/pdf/`,
so it is always identical to its master (test `docs_physics_pdf_copies`). Do not edit or replace it by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
