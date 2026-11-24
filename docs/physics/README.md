# SIMC++ physics notes and the SIMC++ Physics Reference

One physics note per reaction of SIMC++, for physicists who use or extend the program. A note describes the reaction
as SIMC++ builds it (the reaction and physics-model plug-ins, how the pipeline selects them, the stages they use) and
then every formula the code evaluates: generated variables and generation volume, event completion, jacobians, the
cross-section model, radiation, inputs and outputs, validation, and the Fortran origin with its findings. Every
displayed equation names the C++ function that computes it first and the Fortran SIMC line (upstream `59ccdf6`) it
was translated from second; the build refuses an equation without both citations.

The same sources give three products:

| product | where | how it is made |
|---|---|---|
| one standalone PDF per note (masters) | [pdf/](pdf/) | `physics_pdf` |
| the bound book **SIMC++ Physics Reference** (note 00 as Part I, the reactions as the chapters of Part II) | [pdf/SIMCpp_Physics_Reference.pdf](pdf/SIMCpp_Physics_Reference.pdf) | `physics_pdf`, from [latex/physics_reference.tex](latex/physics_reference.tex) |
| a copy of each reaction note in every plug-in folder of its reaction | `src/physics/<reaction>/` (and `src/reactions/<reaction>/` where it exists) | `physics_pdf`, byte-identical to the master |

| note | master PDF | copies | source |
|---|---|---|---|
| 00 Common elements: plug-in architecture, event chain, notation and glossary, frames, kinematics, generation and weighting, target effects, Coulomb correction, radiative corrections, Mott cross section and proton form factors | [SIMCpp_Physics_00_common.pdf](pdf/SIMCpp_Physics_00_common.pdf) | (none) | [00_common.tex](latex/00_common.tex) |
| 01 H(e,e'p) elastic (`hydrogen_elastic`, `sigep`) | [SIMCpp_Physics_01_hydrogen_elastic.pdf](pdf/SIMCpp_Physics_01_hydrogen_elastic.pdf) | [src/physics/proton](../../src/physics/proton/README.md) | [01_hydrogen_elastic.tex](latex/01_hydrogen_elastic.tex) |
| 02 D(e,e'p) and A(e,e'p) (`deuterium_eep`, `nucleus_eep`, `deforest`, spectral functions `theory`, `benhar`) | [SIMCpp_Physics_02_quasielastic_eep.pdf](pdf/SIMCpp_Physics_02_quasielastic_eep.pdf) | [src/physics/proton](../../src/physics/proton/README.md), [src/physics/nuclear](../../src/physics/nuclear/README.md) | [02_quasielastic_eep.tex](latex/02_quasielastic_eep.tex) |
| 03 A(e,e')A nuclear elastic (`nuclear_elastic`, `nuc_elastic`) | [SIMCpp_Physics_03_nuclear_elastic.pdf](pdf/SIMCpp_Physics_03_nuclear_elastic.pdf) | [src/physics/nuclear](../../src/physics/nuclear/README.md) | [03_nuclear_elastic.tex](latex/03_nuclear_elastic.tex) |

Notes 04-09 (pion, kaon, Delta, SIDIS, rho, phase space) follow the same template; they are added to
`latex/notes.json` when they are converted (see "Converting a draft" below).

The user manual (chapter 5, Reactions) links each reaction to its note; the developer manual (chapter 4) and
`docs/manual/markdown/reference/EXTENDING.md` require a physics note for every new reaction or physics-model plug-in.

## Layout

```
docs/physics/
├── README.md              this file
├── CMakeLists.txt         targets physics_pdf (publish), physics_pdf_build (dependency of `pdf`); tests
│                          docs_physics_pdf_build, docs_physics_pdf_copies, docs_physics_template_lint (label docs)
├── build_physics_pdf.py   lint, xelatex passes and log checks, the book, publication of the masters and the copies
├── latex/
│   ├── notes.json         the notes in order (source, PDF name, copy folders) and the book
│   ├── simcnote.cls       the class: standalone note (article) or, with [book], the book (class of the manuals)
│   ├── physics_reference.tex   the book: front matter, Part I = 00_common, Part II = the reaction notes
│   ├── template_reaction.tex   the template of a new note (nine sections, instructions as comments)
│   ├── frag_eep_columns.tex    a fragment shared by notes 01-03 (the (e,e'p) tree columns)
│   └── NN_<name>.tex      one note per file; 00_common.tex is built first
└── pdf/                   the published PDFs (masters and the book)
```

## Single source

A note file is a complete LaTeX document whose preamble contains only `\documentclass{simcnote}` and, for a reaction
note, `\simccommon`; everything else (`\simcnotemeta`, title, text) comes after `\begin{document}`. The standalone build
compiles the file as it is. The book is `\documentclass[book]{simcnote}`: `\simcinputnote{<file>}` inputs the same file
without its preamble, and `\simcnotemeta` then opens a chapter instead of a title band. Labels are written `\lbl{...}`
and stored with the note number (`N01-eq:...`), so the labels of all notes coexist in the book; the standalone reaction
notes read the labels of note 00 from `00_common.aux` (xr-hyper).

## Writing a note

Start from [latex/template_reaction.tex](latex/template_reaction.tex); the comments at its top list the steps. Rules:

* sections, in this order: 1 The reaction in SIMC++ (process; components table, one row per component; selection and
  JSON override; event flow), 2 Kinematics, generated variables, generation volume, 3 Completing the event,
  4 Jacobians, 5 Cross-section model (what a model must provide; the model's formulas), 6 Radiation, energy loss,
  decays, 7 Inputs and outputs, 8 Validation in SIMC++, 9 Fortran origin and Fortran findings;
* displays: `\begin{feq}{<C++ file: function>}{<Fortran file:lines>}` (or `feqa` for an align block); the Fortran
  citation may be `none` for a SIMC++-only formula;
* labels `\lbl{eq:x}`; references `\eqr`, `\secr`, `\tabr`; to note 00 `\ceqr`, `\csecr`, `\ctabr`;
* literature: `\begin{simcrefs} \simcref{key} ... \end{simcrefs}`, cited with `\citer{key}`;
* tables: one item per row (one input, one tree column, one routine, one finding, one test); `L{w}` columns whose
  widths sum to `\linewidth`; file paths in cells with `\file{...}` (breakable);
* text: physics notation and the symbols and terms of note 00 (tables Symbols and Glossary); the Fortran names only in
  tables and citations; finding ids (`D-0xx`, `U-..`, `Pn.Fn`) only in the `simcfindings` table of section 9;
* the reaction prime is `\pr`: `H(e,e\pr p)` (a text apostrophe; `$'$` detaches in sans-serif titles);
* boxes: `simcnote`, `simcwarning`, `simclegacy`, `simcfindingbox` (a candidate finding not yet in
  FORTRAN_DIFFERENCES.md).

Then add the note to `latex/notes.json` with its `copies`, add or update the `README.md` of each copy folder (it must
name the PDF), and publish.

## Converting a draft written for the first template (notes 04-09)

The conversion is mechanical:

1. move `\simcnotemeta{...}` after `\begin{document}`; keep `\simccommon` in the preamble;
2. replace `\title`-like uses of `$'$` by `\pr` (e.g. `H(e,e\pr p)`);
3. every `\begin{feq}{<Fortran>}` becomes `\begin{feq}{<C++ file: function>}{<Fortran>}` (same for `feqa`): look the
   function up in the `// Fortran: file.f:routine` comments of the C++ sources;
4. `\label{` becomes `\lbl{`; `\ref{C-x}` / `\ceqr{x}` keep working with `\ceqr`, `\csecr`, `\ctabr` (drop any `C-`);
5. `thebibliography`/`\bibitem{k}` become `simcrefs`/`\simcref{k}`; `\cite{k}` becomes `\citer{k}`;
6. reorder the sections as above (add "The reaction in SIMC++" and "Validation in SIMC++"; the Fortran map and the
   findings go to section 9);
7. split every table row that lists several items; move finding ids out of the text into the section-9 table
   (`simcfindings` environment) and describe the fact in words where it was;
8. `python3 -I docs/physics/build_physics_pdf.py --source-root . --lint-file NN_<name>.tex` lists what is left.

## Building

```bash
cmake --build <build> --target physics_pdf        # builds and publishes docs/physics/pdf/*.pdf, the book and the copies
cmake --build <build> --target pdf                # the manuals; also builds the physics notes into <build>/docs
ctest --test-dir <build> -R docs_physics          # the build in its own directory, the copy check, the template lint
python3 -I docs/physics/build_physics_pdf.py --source-root . --build-dir /tmp/simcpp_physics [--note 02_quasielastic_eep.tex]
```

Requirements: XeLaTeX with the packages listed at the top of `latex/simcnote.cls` (TeX Live), the fonts of
`docs/common/latex/simcfonts.sty` and TeX Gyre Pagella Math, and `pdfinfo` (optional, for the page count).

The copies in the source folders are written only by `physics_pdf` from the masters in `pdf/`, so they cannot differ;
never edit or replace a copy by hand, and commit masters and copies together.
