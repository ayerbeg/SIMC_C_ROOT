# src/physics/proton: proton cross-section models

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `sigep` | `IPhysicsModel` | `proton::SigEPModel` | `ProtonModels.cpp` | H(e,e'p) elastic cross section (Rosenbluth, Bosted form factors) |
| `deforest` | `IPhysicsModel` | `proton::DeForestModel` | `ProtonModels.cpp` | de Forest off-shell cross section, flavour from `deForest_flag` |
| `deforest_cc1` | `IPhysicsModel` | `proton::DeForestCC1` | `ProtonModels.cpp` | de Forest sigma_cc1 |
| `deforest_cc2` | `IPhysicsModel` | `proton::DeForestCC2` | `ProtonModels.cpp` | de Forest sigma_cc2 |
| `deforest_cc1_onshell` | `IPhysicsModel` | `proton::DeForestCC1OnShell` | `ProtonModels.cpp` | de Forest sigma_cc1 on shell |
| (kernels) | -- | `proton::sigMott`, `fofa_best_fit`, `sigep`, `deForest` | `ProtonCrossSections.cpp` | the formulas (Fortran physics_proton.f) |

The headers are in `include/simc/physics/proton/`.

## Physics notes

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_01_hydrogen_elastic.pdf` | H(e,e'p) elastic scattering: the reaction `hydrogen_elastic` with the model `sigep` | `docs/physics/latex/01_hydrogen_elastic.tex` |
| `SIMCpp_Physics_02_quasielastic_eep.pdf` | D(e,e'p) and A(e,e'p): the reactions `deuterium_eep`, `nucleus_eep` with the model `deforest` and the spectral functions of `src/physics/nuclear/` | `docs/physics/latex/02_quasielastic_eep.tex` |

These PDFs are **generated copies**: the build target `physics_pdf` writes them from the masters in `docs/physics/pdf/`,
so they are always identical to the masters (test `docs_physics_pdf_copies`). Do not edit or replace them by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
