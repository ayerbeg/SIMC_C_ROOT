# src/physics/rho: exclusive rho0 electroproduction kernels

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `rho` | `IPhysicsModel` | `rho::RhoModel` | `RhoModel.cpp` | p(e,e'rho0)p cross section of the reaction `rho_production` (Fortran rho_physics.f:peerho) |
| (kernel) | -- | `rho::peerho` | `RhoPhysics.cpp` | flux times dsigma/dOmega* with the t slope, W and Q2 dependence and R = sigma_L/sigma_T |
| (kernel) | -- | `rho::generateRho` | `RhoPhysics.cpp` | Breit-Wigner mass, flat 4pi direction in the centre of mass, boost to the laboratory (generate_rho.f) |
| (kernel) | -- | `rho::rhoDecay` | `RhoPhysics.cpp` | rho0 -> pi pi decay and the slopes of the detected pion (rho_decay.f) |

The headers are in `include/simc/physics/rho/`; the reaction is in `src/reactions/rho/`.

## Physics note

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_07_rho.pdf` | exclusive rho0 electroproduction: the reaction `rho_production` with the model `rho`, the rho mass, decay and cross section | `docs/physics/latex/07_rho.tex` |

This PDF is a **generated copy**: the build target `physics_pdf` writes it from the master in `docs/physics/pdf/`,
so it is always identical to its master (test `docs_physics_pdf_copies`). Do not edit or replace it by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
