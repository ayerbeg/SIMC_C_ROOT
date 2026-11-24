# src/physics/kaon: kaon electroproduction cross section

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `kaon` | `IPhysicsModel` | `kaon::KaonModel` | `KaonModel.cpp` | N(e,e'K+)Y cross section of the reaction `kaon_production`: the factorized model in the weight, the Saghai model for the `saghai` column, the kaon survival probability |
| (kernel) | -- | `kaon::peeK` | `KaonCrossSections.cpp` | five-fold laboratory cross section and survival probability (Fortran physics_kaon.f:peeK) |
| (kernel) | -- | `kaon::sigFactorized` | `KaonCrossSections.cpp` | factorized fit in Q2, W, t (physics_kaon.f:sig_factorized) |
| (kernel) | -- | `kaon::eekeek` | `KaonCrossSections.cpp` | Saghai amplitude model, Lambda (physics_kaon.f:eekeek) |
| (kernel) | -- | `kaon::eekeeks` | `KaonCrossSections.cpp` | Saghai amplitude model, Sigma (physics_kaon.f:eekeeks) |
| (kernel) | -- | `kaon::phspwght` | `KaonCrossSections.cpp` | phase-space weight routine, not called by SIMC (physics_kaon.f:phspwght) |
| (kernel) | -- | `kaon::fint` | `KaonCrossSections.cpp` | CERNLIB multilinear interpolation (cern/fint.f) |
| (tables) | -- | `kaon::SaghaiTables` | `KaonCrossSections.cpp` | the four data files of `data/kaon` (dbase.f) |

The headers are in `include/simc/physics/kaon/`; the reaction is in `src/reactions/kaon/`.

## Physics note

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_06_kaon.pdf` | kaon electroproduction: the reaction `kaon_production` with the model `kaon`, Lambda, Sigma0, Sigma- and coherent final states, kaon decay in flight | `docs/physics/latex/06_kaon.tex` |

This PDF is a **generated copy**: the build target `physics_pdf` writes it from the master in `docs/physics/pdf/`,
so it is always identical to its master (test `docs_physics_pdf_copies`). Do not edit or replace it by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
