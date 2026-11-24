# src/reactions/pion: exclusive pion production

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `pion_production` | `IReaction` | `reactions::PionProduction` | `PionProduction.cpp` | H, D, A(e,e'pi+-), (e,e'pi0): generation, event completion, weight (Fortran event.f, doing_pion branches); base class of the kaon, Delta, SIDIS and rho reactions |
| (class) | -- | `reactions::DeltaMassTable` | `PionProduction.cpp` | Delta(1232) recoil mass for which_pion 2/3 (delta_relativistic_bw.inp) |

The header is `include/simc/reactions/PionProduction.h`; the cross section is the plug-in `pion` of `src/physics/pion/`.

## Physics notes

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_04_pion.pdf` | Exclusive pion electroproduction: the reaction `pion_production` with the model `pion` | `docs/physics/latex/04_pion.tex` |

These PDFs are **generated copies**: the build target `physics_pdf` writes them from the masters in `docs/physics/pdf/`,
so they are always identical to the masters (test `docs_physics_pdf_copies`). Do not edit or replace them by hand;
change the LaTeX source and rebuild. The common elements (architecture, notation, weighting, radiation, energy loss)
are in note 00 and the whole set is bound as the *SIMC++ Physics Reference* (`docs/physics/pdf/`); see
`docs/physics/README.md`.
