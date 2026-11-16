# The SIMC++ documentation tree

Three documents, each with a **Markdown edition** (readable in the repository browser) and a **LaTeX edition** (the PDF
books). Both editions are sources: edit the one you need, and keep the other in step when the change matters to its
readers. The Doxygen API reference is generated from the headers.

```
docs/
├── README.md                 this map
├── common/                   shared by every edition
│   ├── latex/                simcpp.cls (book class), simcfig.sty (palette, TikZ styles), simcfonts.sty (typefaces),
│   │                         simcpp.ist (index style), conventions.tex (front matter of all books),
│   │                         books.json (manifest for tools/md2latex.py), CMakeLists.txt (target pdf, test docs_pdf_build)
│   └── figures/              simcpp_logo.png
├── manual/                   the user manual and the developer manual
│   ├── markdown/
│   │   ├── user/             user manual (index.md, chapters 01-10, VERIFIED_EXAMPLES.md, examples/*.inp)
│   │   ├── developer/        developer manual (index.md, chapters 01-06)
│   │   └── reference/        framework references: ARCHITECTURE.md, DRIVER.md, EXTENDING.md, OUTPUT_FORMAT.md,
│   │                         INP_JSON_MAPPING.md, PERFORMANCE.md, summary.schema.json
│   ├── latex/                user_manual.tex, developer_manual.tex (main files), user/*.tex, developer/*.tex
│   │                         (chapters), front/*.tex (PDF-only front matter)
│   └── figures/              TikZ sources *.tex and their PNG exports *.png (used by the Markdown edition)
├── translation/              the translation and validation report (Fortran SIMC -> SIMC++)
│   ├── markdown/             index.md, chapters 01-06, appendices, and FORTRAN_DIFFERENCES.md (the living list of
│   │                         findings, source of the catalogue 06_catalogue.md)
│   ├── latex/                translation_report.tex (main file), chapters/*.tex, front/*.tex
│   ├── figures/              TikZ sources and PNG exports (compat_switch, validation_ladder; also used by the manuals)
│   └── tools/                gen_catalogue.py (06_catalogue.md from FORTRAN_DIFFERENCES.md),
│                             gen_fortran_map.py (fortran_routine_map.md)
└── doxygen/                  Doxygen configuration: mainpage.md, groups.dox, CMakeLists.txt (target doc),
                              check_doc_coverage.py, run_doxygen.py; the Doxyfile.in is at the top of the source tree
```

The reference documents in `manual/markdown/reference/` are linked from the developer manual (chapters 3-5 summarise
them; the references keep the full detail: interface contracts, JSON pipeline schema, output formats, the `.inp` to
JSON mapping, performance measurements) and are not part of the PDF books. `summary.schema.json` is also read by tests.

## Editing the Markdown edition

Plain GitHub Markdown. The subset understood by the converter is listed at the top of `tools/md2latex.py`. Check the
links with `python3 -I tools/check_md_links.py docs README.md`, or run `ctest -L docs` in the build directory.
Figures: a Markdown image `![caption](../figures/x.png)` points to the PNG export of `figures/x.tex`.

## Editing the LaTeX edition

The `.tex` files in `manual/latex/` and `translation/latex/` are ordinary, hand-editable LaTeX. Each book has a main
file (`user_manual.tex`, `developer_manual.tex`, `translation_report.tex`) with the parts and the `\input` list of the
chapters; edit the structure there. The chapters use the macros of `common/latex/simcpp.cls` (`\code`, `\simcabstract`,
the `simctable`, `simcshell`, `simccode`, `simcnote` ... environments): copy from an existing chapter.

Two files are produced at build time and are never edited: `simcpp_version.tex` (version and date on the cover) and
`translation_report_dashboard.tex` (the "translation at a glance" page, whose numbers are read from the catalogue).

### Re-converting a chapter from Markdown (optional)

`tools/md2latex.py` created the chapters once. After a large change of a Markdown chapter it can re-create the LaTeX
chapter:

```bash
python3 -I tools/md2latex.py --source-root . --book user_manual --write-sources --chapter 04_inputs.md
```

This **overwrites** the `.tex` file and therefore any hand edit made to it: look at `git diff` and re-apply the edits.
Without `--chapter` every chapter of the book is rewritten; `--skeleton` also rewrites the main file. The chapter
`translation/latex/chapters/06_catalogue.tex` mirrors a generated Markdown file and is always re-converted after
`gen_catalogue.py` (see the top of `translation/markdown/FORTRAN_DIFFERENCES.md`).

## Building the PDF books

```bash
cmake --build <build> --target pdf          # -> <build>/docs/*.pdf
ctest --test-dir <build> -R docs_pdf_build  # the same build as a test
```

Both call `tools/build_docs_pdf.py`, which (1) compiles every figure `docs/{manual,translation}/figures/*.tex` alone
into a PDF (`preview.sty` crops it), (2) writes the build-time files above (`tools/md2latex.py --generated`), and (3) runs
`xelatex`, `makeindex -s docs/common/latex/simcpp.ist`, then `xelatex` until the references are stable, and fails on
a LaTeX error, an undefined reference or a missing glyph. Requirements: XeLaTeX and makeindex (TeX Live, with the
packages listed at the top of `simcpp.cls`, and `preview`), the fonts TeX Gyre Pagella, TeX Gyre Heros and
DejaVu Sans Mono, and `pdftoppm` for the PNG export.

By hand, for example for the user manual, from the source tree:

```bash
D=$PWD/docs; OUT=/tmp/simcpp_pdf; mkdir -p $OUT/figures
export TEXINPUTS=$D/manual/latex:$D/common/latex//:$D/common/figures:$OUT/figures:$OUT:
for f in $D/manual/figures/*.tex $D/translation/figures/*.tex; do xelatex -output-directory $OUT/figures $f; done
python3 -I tools/md2latex.py --source-root . --book user_manual --generated $OUT     # optional
cd $OUT && xelatex user_manual && makeindex -s $D/common/latex/simcpp.ist user_manual && xelatex user_manual && xelatex user_manual
```

After changing a figure, re-export its PNG for the Markdown edition and commit both files:
`python3 -I tools/build_docs_pdf.py --source-root . --build-dir <dir> --figures-only --png`.

## The API reference

`cmake --build <build> --target doc` runs Doxygen (configuration `Doxyfile.in` + `docs/doxygen/`); the HTML is in
`<build>/docs/doxygen/html/index.html`. With `SIMC_DOXYGEN_MANUALS=ON` (default) the Markdown editions are rendered
as pages too.
