# tests/user_manual — checks of the user manual docs/manual/markdown/user/ (task P15.D2), label "docs".
#   docs_user_manual_parameters: every SimcParameters input variable has a row in the parameter reference
#                                (04_inputs.md) and every `simc++ --help` option is documented in the manual;
#   docs_user_manual_links:      every relative Markdown link (file and #anchor) of the user manual and README.md
#                                resolves, and those of the two user guides in docs/translation/markdown/ (tools/check_md_links.py, reusable for any Markdown file or directory);
#   docs_no_machine_paths:       no "/home/" or "NewDisk" in any file under docs/ or in README.md (D-073).
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME docs_user_manual_parameters
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_CURRENT_SOURCE_DIR}/check_parameters.py $<TARGET_FILE:inp2json>
                     $<TARGET_FILE:simcpp> ${CMAKE_SOURCE_DIR}/docs/manual/markdown/user/04_inputs.md
                     ${CMAKE_SOURCE_DIR}/docs/manual/markdown/user)
    add_test(NAME docs_user_manual_links
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_SOURCE_DIR}/tools/check_md_links.py
                     ${CMAKE_SOURCE_DIR}/docs/manual/markdown/user ${CMAKE_SOURCE_DIR}/README.md
                     ${CMAKE_SOURCE_DIR}/docs/translation/markdown/user_comparing_with_fortran.md
                     ${CMAKE_SOURCE_DIR}/docs/translation/markdown/user_differences.md)
    # D-073: no machine-specific paths ("/home/", "NewDisk") in any file under docs/ or in README.md
    add_test(NAME docs_no_machine_paths
             COMMAND ${Python3_EXECUTABLE} -I ${CMAKE_CURRENT_SOURCE_DIR}/check_no_machine_paths.py
                     ${CMAKE_SOURCE_DIR}/docs ${CMAKE_SOURCE_DIR}/README.md)
    set_tests_properties(docs_user_manual_parameters docs_user_manual_links docs_no_machine_paths
                         PROPERTIES LABELS docs TIMEOUT 120)
else()
    message(STATUS "Python 3 not found: the docs tests (tests/user_manual) are not registered")
endif()
