# SIMC++ physics notes

One short PDF per reaction of SIMC++, for physicists who use the program: what the reaction simulates, which
variables are generated and over which ranges, how the event is completed, the jacobians, the cross-section model
with every formula the code evaluates, how radiation, energy loss and decays apply, the inputs and tree columns of the
reaction, the Fortran-to-C++ map and the Fortran findings that touch it. A common note collects what all reactions
share and is referenced by the others.

Every displayed formula is transcribed from the Fortran SIMC (`simc_gfortran`, upstream commit `59ccdf6`, the
specification SIMC++ follows) and cites the file and line it comes from; the build refuses an equation without a
citation.

| note | PDF (master) | copy | LaTeX source |
|---|---|---|---|
| 00 Common elements: frames, kinematic variables, event chain, generation limits, weighting and normalisation, target effects, energy loss and multiple scattering, Coulomb correction, radiative corrections | [pdf/SIMCpp_Physics_00_common.pdf](pdf/SIMCpp_Physics_00_common.pdf) | (none) | [latex/00_common.tex](latex/00_common.tex) |
| 01 H(e,e'p) elastic (`hydrogen_elastic`, `sigep`) | [pdf/SIMCpp_Physics_01_hydrogen_elastic.pdf](pdf/SIMCpp_Physics_01_hydrogen_elastic.pdf) | [src/physics/proton/](../../src/physics/proton/README.md) | [latex/01_hydrogen_elastic.tex](latex/01_hydrogen_elastic.tex) |
| 02 D(e,e'p) and A(e,e'p) quasi-elastic (`deuterium_eep`, `nucleus_eep`, `deforest`) | [pdf/SIMCpp_Physics_02_quasielastic_eep.pdf](pdf/SIMCpp_Physics_02_quasielastic_eep.pdf) | [src/physics/nuclear/](../../src/physics/nuclear/README.md) | [latex/02_quasielastic_eep.tex](latex/02_quasielastic_eep.tex) |
| 03 A(e,e')A nuclear elastic (`nuclear_elastic`, `nuc_elastic`) | [pdf/SIMCpp_Physics_03_nuclear_elastic.pdf](pdf/SIMCpp_Physics_03_nuclear_elastic.pdf) | [src/physics/nuclear/](../../src/physics/nuclear/README.md) | [latex/03_nuclear_elastic.tex](latex/03_nuclear_elastic.tex) |

The notes for pion production, kaon production, the Delta region, SIDIS, rho production and phase space follow the
same model (later tasks); their PDFs will be copied to `src/physics/<reaction>/`.

## Layout

```
docs/physics/
├── README.md              this file
├── CMakeLists.txt         target physics_pdf; tests docs_physics_pdf_build, docs_physics_pdf_copies (label docs)
├── build_physics_pdf.py   lint (code citations), xelatex passes, log checks, publication of the PDFs and their copies
├── latex/
│   ├── notes.json         the list of notes: source, PDF name, the source-tree folders that receive a copy
│   ├── simcnote.cls       the document class (look of the SIMC++ books: docs/common/latex/simcfig.sty, simcfonts.sty)
│   ├── frag_*.tex         fragments shared by several notes (the (e,e'p) tree layout)
│   └── NN_<name>.tex      one note per file; 00_common.tex is built first
└── pdf/                   the published PDFs (masters)
```

## Writing a note

Copy `latex/01_hydrogen_elastic.tex` and keep its nine sections. Conventions of `simcnote.cls`:

* displayed equations only with `\begin{feq}{event.f:518} ... \end{feq}` (one line) or `\begin{feqa}{...} ... \end{feqa}`
  (an `align` block); the argument is the code citation printed under the equation. `equation`, `align`, `\[ \]` and
  `$$` are refused by the build;
* `\fsrc{file.f:line}` cites the Fortran, `\csrc{dir/File.cpp: Function}` the C++, `\finding{D-033 F3}` a finding of
  `docs/translation/markdown/FORTRAN_DIFFERENCES.md`;
* `\eqr{eq:x}`, `\secr{sec:x}`, `\tabr{tab:x}` refer inside the note; `\simccommon` in the preamble loads the labels of
  the common note, referred to with `\ceqr`, `\csecr`, `\ctabr`;
* boxes: `simcnote`, `simcwarning`, `simclegacy` (a Fortran behaviour), `simcfinding` (a candidate finding not yet in
  FORTRAN_DIFFERENCES.md).

Add the note to `latex/notes.json` with its `copies` (the reaction's `src/physics/<reaction>/` folder) and add or update
the `README.md` of that folder (it must name the PDF; `docs_physics_pdf_copies` checks it).

## Building

```bash
cmake --build <build> --target physics_pdf        # builds and publishes docs/physics/pdf/*.pdf + the copies
ctest --test-dir <build> -R docs_physics_pdf      # the build in its own directory + the copy check
```

or by hand, without publishing:

```bash
python3 -I docs/physics/build_physics_pdf.py --source-root . --build-dir /tmp/simcpp_physics [--note 02_quasielastic_eep.tex]
```

Requirements: XeLaTeX with the packages listed at the top of `latex/simcnote.cls` (TeX Live), the fonts of
`docs/common/latex/simcfonts.sty` and TeX Gyre Pagella Math, and `pdfinfo` (optional, for the page count).

The copies in `src/physics/<reaction>/` are written only by `physics_pdf` from the master in `pdf/`, so the two can
never differ; never edit or replace a copy by hand. Commit the masters and the copies together.
