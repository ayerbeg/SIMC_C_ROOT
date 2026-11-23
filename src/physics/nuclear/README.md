# src/physics/nuclear: physics notes

Two physics notes belong to the code of this folder:

* `SIMCpp_Physics_02_quasielastic_eep.pdf`: D(e,e'p) and A(e,e'p) quasi-elastic knock-out (reactions `deuterium_eep`,
  `nucleus_eep`, model `deforest`): generation and limits, the proton-energy solution and jacobian for deuterium, the
  de Forest cross sections, the spectral functions (theory files, Benhar tables, implemented here in
  `SpectralFunctions.cpp` and `SpectralFunctionTables.cpp`), transparency, radiation, inputs, tree columns, the
  Fortran-to-C++ map and the Fortran findings;
* `SIMCpp_Physics_03_nuclear_elastic.pdf`: A(e,e')A nuclear elastic scattering (reaction `nuclear_elastic`, model
  `nuc_elastic`, implemented here in `NuclearElastic.cpp`): kinematics, the charge form factor, radiation, outputs and
  the upstream status of the reaction.

These PDFs are generated copies. The sources are `docs/physics/latex/02_quasielastic_eep.tex` and
`docs/physics/latex/03_nuclear_elastic.tex`, the master PDFs are in `docs/physics/pdf/`; the build target
`physics_pdf` writes masters and copies together, so they are always identical (test `docs_physics_pdf_copies`). Do not
edit or replace these files by hand. The common elements (frames, weighting, radiation, energy loss) are in
`docs/physics/pdf/SIMCpp_Physics_00_common.pdf`; see `docs/physics/README.md`.
