# docs/doxygen — the Doxygen API reference of SIMC++ (task P15.D3, decision D-072).
#
#   target `doc`:          runs Doxygen on the configured Doxyfile (Doxyfile.in at the top of the source tree, filled in by
#                          configure_file: no machine-specific path in the repository, D-073); HTML in
#                          <build>/docs/doxygen/html/index.html. Without Doxygen the target only prints a message.
#   CTest docs_doc_coverage (label docs): check_doc_coverage.py — every public entity of include/ has a doc comment, every
#                          documented function its @param / @return (works without Doxygen).
#   CTest docs_doxygen     (label docs): run_doxygen.py — runs Doxygen and requires 0 warnings for include/; skipped
#                          (exit 77) when Doxygen is not installed.
# Options: SIMC_DOXYGEN_MANUALS (ON) also renders the Markdown manuals (docs/manual, docs/translation and the reference
# documents of docs/) as pages, so that the main page's links to them resolve inside the HTML.
find_package(Doxygen QUIET OPTIONAL_COMPONENTS dot)
option(SIMC_DOXYGEN_MANUALS "Render the Markdown manuals as pages of the Doxygen reference" ON)

set(SIMC_DOXYGEN_OUT "${CMAKE_CURRENT_BINARY_DIR}")
set(SIMC_DOXYGEN_MANUAL_INPUTS "")
if(SIMC_DOXYGEN_MANUALS)
    foreach(_d manual/user manual/developer translation ARCHITECTURE.md EXTENDING.md DRIVER.md PERFORMANCE.md
               OUTPUT_FORMAT.md INP_JSON_MAPPING.md FORTRAN_DIFFERENCES.md)
        string(APPEND SIMC_DOXYGEN_MANUAL_INPUTS " \"${CMAKE_SOURCE_DIR}/docs/${_d}\"")
    endforeach()
endif()
if(DOXYGEN_DOT_EXECUTABLE)
    set(SIMC_DOXYGEN_HAVE_DOT YES)
    get_filename_component(SIMC_DOXYGEN_DOT_PATH "${DOXYGEN_DOT_EXECUTABLE}" DIRECTORY)
else()
    set(SIMC_DOXYGEN_HAVE_DOT NO)
    set(SIMC_DOXYGEN_DOT_PATH "")
endif()
configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

if(DOXYGEN_EXECUTABLE)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        COMMAND ${CMAKE_COMMAND} -E echo "SIMC++ API reference: ${CMAKE_CURRENT_BINARY_DIR}/html/index.html (warnings: ${CMAKE_CURRENT_BINARY_DIR}/doxygen_warnings.log)"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Doxygen API reference of SIMC++"
        VERBATIM)
    set(_doxygen_exe "${DOXYGEN_EXECUTABLE}")
else()
    add_custom_target(doc
        COMMAND ${CMAKE_COMMAND} -E echo "doxygen not found: the API reference is not built (install doxygen and graphviz, then re-run cmake)"
        VERBATIM)
    set(_doxygen_exe "doxygen-not-found")
endif()

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_doc_coverage
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_CURRENT_SOURCE_DIR}/check_doc_coverage.py ${CMAKE_SOURCE_DIR}/include)
    add_test(NAME docs_doxygen
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_CURRENT_SOURCE_DIR}/run_doxygen.py ${_doxygen_exe}
                     ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doxygen_warnings.log
                     ${CMAKE_SOURCE_DIR}/include)
    set_tests_properties(docs_doc_coverage docs_doxygen PROPERTIES LABELS docs TIMEOUT 900)
    set_tests_properties(docs_doxygen PROPERTIES SKIP_RETURN_CODE 77)
else()
    message(STATUS "Python 3 not found: the API documentation tests (docs/doxygen) are not registered")
endif()
