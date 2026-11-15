# tests/developer_manual — checks of the developer manual docs/manual/developer/ (task P15.D1), label "docs".
#   docs_dev_manual:       every header of include/simc is covered, every stage-interface class is described in chapter 3,
#                          every `<!-- from: path -->` snippet is a verbatim copy of compiled code, the generated Fortran
#                          routine map has no uncited routine (check_dev_manual.py);
#   docs_dev_manual_links: every relative Markdown link (file and #anchor) of the developer manual and of the detailed
#                          references it links (docs/ARCHITECTURE.md, EXTENDING.md, PERFORMANCE.md, DRIVER.md,
#                          FORTRAN_DIFFERENCES.md) resolves (tools/check_md_links.py, shared with the user manual).
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_dev_manual
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_CURRENT_SOURCE_DIR}/check_dev_manual.py ${CMAKE_SOURCE_DIR})
    add_test(NAME docs_dev_manual_links
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_SOURCE_DIR}/tools/check_md_links.py
                     ${CMAKE_SOURCE_DIR}/docs/manual/developer ${CMAKE_SOURCE_DIR}/docs/ARCHITECTURE.md
                     ${CMAKE_SOURCE_DIR}/docs/EXTENDING.md ${CMAKE_SOURCE_DIR}/docs/PERFORMANCE.md
                     ${CMAKE_SOURCE_DIR}/docs/DRIVER.md ${CMAKE_SOURCE_DIR}/docs/FORTRAN_DIFFERENCES.md)
    set_tests_properties(docs_dev_manual docs_dev_manual_links PROPERTIES LABELS docs TIMEOUT 120)
else()
    message(STATUS "Python 3 not found: the docs tests (tests/developer_manual) are not registered")
endif()
