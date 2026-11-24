# src/reactions/rho: the rho-production reaction

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `rho_production` | `IReaction` | `RhoProduction` | `RhoProduction.cpp` | p(e,e'rho0)p with one decay pion detected; derives from `PionProduction` (meson kinematics, in `src/reactions/pion/`) and implements the generation, the rho decay, the vertex completion and `completeMain` |

The header is `include/simc/reactions/RhoProduction.h`; the cross-section model `rho` and the generation and decay
kernels are in `src/physics/rho/`.

## Physics note

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_07_rho.pdf` | exclusive rho0 electroproduction: the reaction `rho_production` with the model `rho` | `docs/physics/latex/07_rho.tex` |

This PDF is a **generated copy**: the build target `physics_pdf` writes it from the master in `docs/physics/pdf/`,
so it is always identical to its master (test `docs_physics_pdf_copies`). Do not edit or replace it by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
