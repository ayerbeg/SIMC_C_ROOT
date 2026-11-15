# tests/developer_manual — checks of the developer manual docs/manual/developer/ and of the translation report
# docs/translation/ (task P15.D1), label "docs".
#   docs_dev_manual:             every header of include/simc is described in the developer manual, every stage-interface
#                                class in its chapter 3, every `<!-- from: path -->` snippet is a verbatim copy of compiled
#                                code (check_dev_manual.py);
#   docs_dev_manual_links:       every relative Markdown link (file and #anchor) of the developer manual and of the detailed
#                                references it links (docs/ARCHITECTURE.md, EXTENDING.md, PERFORMANCE.md, DRIVER.md)
#                                resolves (tools/check_md_links.py, shared with the user manual);
#   docs_translation:            the translation report: verbatim snippets, no Fortran routine without a C++ citation or a
#                                reason (fortran_routine_map.md), the catalogue 06_catalogue.md up to date with
#                                docs/FORTRAN_DIFFERENCES.md (tools/gen_catalogue.py);
#   docs_translation_links:      the links of docs/translation/ and docs/FORTRAN_DIFFERENCES.md.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_dev_manual
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_CURRENT_SOURCE_DIR}/check_dev_manual.py ${CMAKE_SOURCE_DIR} developer)
    add_test(NAME docs_dev_manual_links
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_SOURCE_DIR}/tools/check_md_links.py
                     ${CMAKE_SOURCE_DIR}/docs/manual/developer ${CMAKE_SOURCE_DIR}/docs/ARCHITECTURE.md
                     ${CMAKE_SOURCE_DIR}/docs/EXTENDING.md ${CMAKE_SOURCE_DIR}/docs/PERFORMANCE.md
                     ${CMAKE_SOURCE_DIR}/docs/DRIVER.md)
    add_test(NAME docs_translation
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_CURRENT_SOURCE_DIR}/check_dev_manual.py ${CMAKE_SOURCE_DIR} translation)
    add_test(NAME docs_translation_links
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_SOURCE_DIR}/tools/check_md_links.py
                     ${CMAKE_SOURCE_DIR}/docs/translation ${CMAKE_SOURCE_DIR}/docs/FORTRAN_DIFFERENCES.md)
    set_tests_properties(docs_dev_manual docs_dev_manual_links docs_translation docs_translation_links
                         PROPERTIES LABELS docs TIMEOUT 120)
else()
    message(STATUS "Python 3 not found: the docs tests (tests/developer_manual) are not registered")
endif()
