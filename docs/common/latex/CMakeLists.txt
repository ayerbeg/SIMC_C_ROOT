# docs/common/latex - the PDF editions of the user manual, the developer manual and the translation report
# (tasks P15.D5, P15.D6). See docs/README.md.
#
# The LaTeX editions are committed sources: docs/manual/latex/{user_manual,developer_manual}.tex with their chapters in
# user/ and developer/, docs/translation/latex/translation_report.tex with chapters/; the shared style is here
# (simcpp.cls, simcfig.sty, simcfonts.sty, simcpp.ist, conventions.tex), the logo in docs/common/figures/, the TikZ
# figures in docs/{manual,translation}/figures/*.tex. tools/build_docs_pdf.py compiles the figures, writes the
# build-time inputs (version, dashboard; tools/md2latex.py --generated), runs xelatex + makeindex and checks the logs.
# tools/md2latex.py --write-sources is an optional helper that re-converts chapters from the Markdown edition.
#
#   target `pdf`:            <build>/docs/SIMCpp_User_Manual.pdf, SIMCpp_Developer_Manual.pdf,
#                            SIMC_Fortran_to_SIMCpp_Translation_Report.pdf (intermediate files in <build>/docs/latex).
#                            Without xelatex the target only prints a message.
#   CTest docs_pdf_build     (label docs): the same build in its own directory; fails on a LaTeX error, an undefined
#                            reference, a missing file or glyph; SKIP (77) without xelatex.
#   Re-export the figure PNGs after changing a figure (writes next to the .tex sources):
#       python3 -I tools/build_docs_pdf.py --source-root . --build-dir build/docs/latex --figures-only --png
find_program(SIMC_XELATEX xelatex)
find_package(Python3 COMPONENTS Interpreter)
set(_simc_pdf_script ${CMAKE_SOURCE_DIR}/tools/build_docs_pdf.py)
if(SIMC_XELATEX AND Python3_Interpreter_FOUND)
    add_custom_target(pdf
        COMMAND ${Python3_EXECUTABLE} -I ${_simc_pdf_script} --source-root ${CMAKE_SOURCE_DIR}
                --build-dir ${CMAKE_BINARY_DIR}/docs/latex --pdf-dir ${CMAKE_BINARY_DIR}/docs
        COMMENT "PDF manuals (xelatex): ${CMAKE_BINARY_DIR}/docs/*.pdf"
        VERBATIM)
else()
    add_custom_target(pdf
        COMMAND ${CMAKE_COMMAND} -E echo "pdf: xelatex (or Python 3) not found: the PDF manuals are not built"
        VERBATIM)
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_pdf_build
             COMMAND ${Python3_EXECUTABLE} -I ${_simc_pdf_script} --source-root ${CMAKE_SOURCE_DIR}
                     --build-dir ${CMAKE_CURRENT_BINARY_DIR}/ctest_pdf)
    set_tests_properties(docs_pdf_build PROPERTIES LABELS docs SKIP_RETURN_CODE 77 TIMEOUT 900)
endif()
