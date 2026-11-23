# src/physics/proton: physics note

`SIMCpp_Physics_01_hydrogen_elastic.pdf` is the physics note of the H(e,e'p) elastic reaction (reaction
`hydrogen_elastic`, model `sigep`, implemented in this folder by `ProtonCrossSections.cpp` and `ProtonModels.cpp`): the
generated variables and the generation volume, the elastic kinematics, the jacobian, the Rosenbluth cross section with
the Bosted form-factor fit, the radiative tails, the inputs and tree columns, the Fortran-to-C++ map and the Fortran
findings. The de Forest cross sections of the same files are described in note 02 (`src/physics/nuclear/`).

This PDF is a generated copy. The source is `docs/physics/latex/01_hydrogen_elastic.tex` and the master PDF is
`docs/physics/pdf/SIMCpp_Physics_01_hydrogen_elastic.pdf`; the build target `physics_pdf` writes both, so they are
always identical (test `docs_physics_pdf_copies`). Do not edit or replace this file by hand. The common elements
(frames, weighting, radiation, energy loss) are in `docs/physics/pdf/SIMCpp_Physics_00_common.pdf`; see
`docs/physics/README.md`.
