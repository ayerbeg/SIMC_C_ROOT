# docs/physics - the SIMC++ physics notes (task PHYS.1; docs/physics/README.md).
#
#   target `physics_pdf`:  builds every note of docs/physics/latex/notes.json (intermediate files in
#                          <build>/docs/physics) and PUBLISHES the PDFs: docs/physics/pdf/<output>.pdf (the masters) and
#                          the generated copies in the reaction folders listed under "copies" (src/physics/proton/,
#                          src/physics/nuclear/, ...). Commit the masters and the copies together.
#                          Without xelatex the target only prints a message.
#   CTest docs_physics_pdf_build (label docs): the same build in its own directory, nothing published; fails on a LaTeX
#                          error, an undefined reference or glyph, or a displayed equation without a code citation;
#                          SKIP (77) without xelatex.
#   CTest docs_physics_pdf_copies (label docs): every committed copy is byte-identical to its master in docs/physics/pdf.
find_program(SIMC_XELATEX_PHYS xelatex)
find_package(Python3 COMPONENTS Interpreter)
set(_simc_phys_script ${CMAKE_CURRENT_SOURCE_DIR}/build_physics_pdf.py)
if(SIMC_XELATEX_PHYS AND Python3_Interpreter_FOUND)
    add_custom_target(physics_pdf
        COMMAND ${Python3_EXECUTABLE} -I ${_simc_phys_script} --source-root ${CMAKE_SOURCE_DIR}
                --build-dir ${CMAKE_CURRENT_BINARY_DIR} --publish
        COMMENT "SIMC++ physics notes (xelatex): docs/physics/pdf/*.pdf and the copies in src/physics/<reaction>/"
        VERBATIM)
else()
    add_custom_target(physics_pdf
        COMMAND ${CMAKE_COMMAND} -E echo "physics_pdf: xelatex (or Python 3) not found: the physics notes are not built"
        VERBATIM)
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_physics_pdf_build
             COMMAND ${Python3_EXECUTABLE} -I ${_simc_phys_script} --source-root ${CMAKE_SOURCE_DIR}
                     --build-dir ${CMAKE_CURRENT_BINARY_DIR}/ctest_pdf)
    set_tests_properties(docs_physics_pdf_build PROPERTIES LABELS docs SKIP_RETURN_CODE 77 TIMEOUT 600)
    add_test(NAME docs_physics_pdf_copies
             COMMAND ${Python3_EXECUTABLE} -I ${_simc_phys_script} --source-root ${CMAKE_SOURCE_DIR} --check-copies)
    set_tests_properties(docs_physics_pdf_copies PROPERTIES LABELS docs TIMEOUT 60)
endif()
