# docs/physics - the SIMC++ physics notes and the bound book "SIMC++ Physics Reference" (docs/physics/README.md).
#
#   target `physics_pdf`:        builds every note of docs/physics/latex/notes.json and the book, and PUBLISHES them:
#                                docs/physics/pdf/<output>.pdf (masters, the book included) and the generated copies of
#                                each note in its plug-in folders ("copies" in notes.json: src/physics/proton/, ...).
#                                Commit the masters and the copies together. Intermediate files: <build>/docs/physics.
#   target `physics_pdf_build`:  the same build without publishing (PDFs in <build>/docs); a dependency of the target
#                                `pdf` (docs/common/latex), so `pdf` also produces the physics notes and the book.
#   CTest docs_physics_pdf_build (label docs): the build in its own directory, nothing published; fails on a LaTeX
#                                error, an undefined reference or glyph, or a display without its two code citations;
#                                SKIP (77) without xelatex.
#   CTest docs_physics_pdf_copies (label docs): every committed copy is byte-identical to its master in
#                                docs/physics/pdf, every copy folder has a README.md naming the PDF, the book is present.
#   CTest docs_physics_template_lint (label docs): the template template_reaction.tex passes the lint.
find_program(SIMC_XELATEX_PHYS xelatex)
find_package(Python3 COMPONENTS Interpreter)
set(_simc_phys_script ${CMAKE_CURRENT_SOURCE_DIR}/build_physics_pdf.py)
if(SIMC_XELATEX_PHYS AND Python3_Interpreter_FOUND)
    add_custom_target(physics_pdf
        COMMAND ${Python3_EXECUTABLE} -I ${_simc_phys_script} --source-root ${CMAKE_SOURCE_DIR}
                --build-dir ${CMAKE_CURRENT_BINARY_DIR} --publish
        COMMENT "SIMC++ physics notes + Physics Reference (xelatex): docs/physics/pdf/ and the copies in src/physics/<reaction>/"
        VERBATIM)
    add_custom_target(physics_pdf_build
        COMMAND ${Python3_EXECUTABLE} -I ${_simc_phys_script} --source-root ${CMAKE_SOURCE_DIR}
                --build-dir ${CMAKE_CURRENT_BINARY_DIR}/unpublished --pdf-dir ${CMAKE_BINARY_DIR}/docs
        COMMENT "SIMC++ physics notes + Physics Reference (xelatex): ${CMAKE_BINARY_DIR}/docs/SIMCpp_Physics_*.pdf"
        VERBATIM)
else()
    add_custom_target(physics_pdf
        COMMAND ${CMAKE_COMMAND} -E echo "physics_pdf: xelatex (or Python 3) not found: the physics notes are not built"
        VERBATIM)
    add_custom_target(physics_pdf_build
        COMMAND ${CMAKE_COMMAND} -E echo "physics_pdf_build: xelatex (or Python 3) not found"
        VERBATIM)
endif()
if(TARGET pdf)
    add_dependencies(pdf physics_pdf_build)
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_physics_pdf_build
             COMMAND ${Python3_EXECUTABLE} -I ${_simc_phys_script} --source-root ${CMAKE_SOURCE_DIR}
                     --build-dir ${CMAKE_CURRENT_BINARY_DIR}/ctest_pdf)
    set_tests_properties(docs_physics_pdf_build PROPERTIES LABELS docs SKIP_RETURN_CODE 77 TIMEOUT 900)
    add_test(NAME docs_physics_pdf_copies
             COMMAND ${Python3_EXECUTABLE} -I ${_simc_phys_script} --source-root ${CMAKE_SOURCE_DIR} --check-copies)
    set_tests_properties(docs_physics_pdf_copies PROPERTIES LABELS docs TIMEOUT 60)
    add_test(NAME docs_physics_template_lint
             COMMAND ${Python3_EXECUTABLE} -I ${_simc_phys_script} --source-root ${CMAKE_SOURCE_DIR}
                     --lint-file template_reaction.tex)
    set_tests_properties(docs_physics_template_lint PROPERTIES LABELS docs TIMEOUT 60)
endif()
