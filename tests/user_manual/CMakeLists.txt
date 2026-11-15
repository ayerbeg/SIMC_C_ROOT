# tests/user_manual — checks of the user manual docs/manual/user/ (task P15.D2), label "docs".
#   docs_user_manual_parameters: every SimcParameters input variable has a row in the parameter reference
#                                (04_inputs.md) and every `simc++ --help` option is documented in the manual;
#   docs_user_manual_links:      every relative Markdown link (file and #anchor) of the user manual and README.md
#                                resolves, and those of the two user guides in docs/translation/ (tools/check_md_links.py, reusable for any Markdown file or directory).
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_user_manual_parameters
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_CURRENT_SOURCE_DIR}/check_parameters.py $<TARGET_FILE:inp2json>
                     $<TARGET_FILE:simcpp> ${CMAKE_SOURCE_DIR}/docs/manual/user/04_inputs.md
                     ${CMAKE_SOURCE_DIR}/docs/manual/user)
    add_test(NAME docs_user_manual_links
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_SOURCE_DIR}/tools/check_md_links.py
                     ${CMAKE_SOURCE_DIR}/docs/manual/user ${CMAKE_SOURCE_DIR}/README.md
                     ${CMAKE_SOURCE_DIR}/docs/translation/user_comparing_with_fortran.md
                     ${CMAKE_SOURCE_DIR}/docs/translation/user_differences.md)
    set_tests_properties(docs_user_manual_parameters docs_user_manual_links PROPERTIES LABELS docs TIMEOUT 120)
else()
    message(STATUS "Python 3 not found: the docs tests (tests/user_manual) are not registered")
endif()
