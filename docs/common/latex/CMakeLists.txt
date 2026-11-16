# docs/latex — the PDF editions of the user manual, the developer manual and the translation report (task P15.D5).
#
# The Markdown files (docs/manual/user, docs/manual/developer, docs/translation) are the source of the text;
# tools/md2latex.py converts them, docs/latex/simcpp.cls gives the design, docs/latex/figures/*.tex are the TikZ
# figures (their PNG exports in docs/figures/ are used by the Markdown view). Everything is driven by
# tools/build_docs_pdf.py (figures, conversion, xelatex + makeindex passes, log checks).
#
#   target `pdf`:            <build>/docs/SIMCpp_User_Manual.pdf, SIMCpp_Developer_Manual.pdf,
#                            SIMC_Fortran_to_SIMCpp_Translation_Report.pdf (intermediate files in <build>/docs/latex).
#                            Without xelatex the target only prints a message.
#   CTest docs_pdf_build     (label docs): the same build in its own directory; fails on a converter error, a LaTeX
#                            error, an undefined reference or a missing glyph; SKIP (77) without xelatex.
#   Re-export the figure PNGs after changing a figure (writes into the source tree):
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
