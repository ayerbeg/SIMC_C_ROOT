# src/reactions/sidis: semi-inclusive (SIDIS) reaction

## The plug-ins and classes of this folder

| registered name | interface | class | file | role |
|---|---|---|---|---|
| `sidis_production` | `IReaction` | `reactions::SidisProduction` | `SidisProduction.cpp` | (e,e'pi)X / (e,e'K)X on H, D, A: generation of both energies and directions, vertex completion, weight (Fortran event.f, doing_semi branches) |

The class derives from `reactions::PionProduction` (`src/reactions/pion/`) and `reactions::EepReactionBase`
(`src/reactions/`); the header is `include/simc/reactions/SidisProduction.h`. The cross section is the physics model
`sidis` of `src/physics/sidis/`.

## Physics note

| PDF in this folder | note | LaTeX source |
|---|---|---|
| `SIMCpp_Physics_08_sidis.pdf` | semi-inclusive (e,e'h)X production: the reaction `sidis_production` with the model `sidis` | `docs/physics/latex/08_sidis.tex` |

This PDF is a **generated copy**: the build target `physics_pdf` writes it from the master in `docs/physics/pdf/`, so it
is always identical to the master (test `docs_physics_pdf_copies`). Do not edit or replace it by hand; change the LaTeX
source and rebuild. The common elements are in note 00 and the whole set is bound as the *SIMC++ Physics Reference*
(`docs/physics/pdf/`); see `docs/physics/README.md`. Phase-space runs (`doing_phsp`) are described in note 09
(`docs/physics/pdf/SIMCpp_Physics_09_phase_space.pdf`).
