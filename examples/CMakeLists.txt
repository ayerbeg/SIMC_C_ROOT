# examples/ — example plug-ins (docs/EXTENDING.md). Each one is an OBJECT library (so that its SIMC_REGISTER line is
# linked) appended to SIMC_EXAMPLE_PLUGINS, which simc++ links.

# ISpectrometer "ideal_box": perfect resolution, rectangular acceptance (examples/ideal_box/IdealBoxSpectrometer.cpp)
add_library(simcfw_example_ideal_box OBJECT ideal_box/IdealBoxSpectrometer.cpp)
target_link_libraries(simcfw_example_ideal_box PUBLIC simcfw)
target_compile_options(simcfw_example_ideal_box PRIVATE -Wall -Wextra -pedantic)
list(APPEND SIMC_EXAMPLE_PLUGINS simcfw_example_ideal_box)
set(SIMC_EXAMPLE_PLUGINS ${SIMC_EXAMPLE_PLUGINS} PARENT_SCOPE)

# Test: the complete driver with ideal_box on the hadron arm, then on both arms (label "examples").
add_executable(test_example_ideal_box ideal_box/test_ideal_box.cpp)
target_link_libraries(test_example_ideal_box PRIVATE ${SIMC_DRIVER_LIBS} simcfw_example_ideal_box)
target_compile_definitions(test_example_ideal_box PRIVATE SIMC_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
set(_inp ${SIMC_FORTRAN_DIR}/infiles/eep_h_elastic_rad_A.inp)
foreach(_cfg hadron_arm both_arms)
    add_test(NAME example_ideal_box_${_cfg}
             COMMAND test_example_ideal_box ${_inp} ${CMAKE_CURRENT_SOURCE_DIR}/ideal_box/${_cfg}.json
                     ${CMAKE_CURRENT_BINARY_DIR}/ideal_box_${_cfg})
    set_tests_properties(example_ideal_box_${_cfg} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600 LABELS "examples")
endforeach()
# The same through the simc++ command line: the plug-in is chosen purely by the --pipeline JSON.
add_test(NAME example_ideal_box_cli
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ideal_box/cli_check.sh $<TARGET_FILE:simcpp> ${_inp}
                 ${CMAKE_CURRENT_SOURCE_DIR}/ideal_box/hadron_arm.json ${CMAKE_CURRENT_BINARY_DIR}/ideal_box_cli)
set_tests_properties(example_ideal_box_cli PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600 LABELS "examples;cli")

# P15.3: the installable package + an out-of-tree plug-in (examples/out_of_tree_plugin): cmake --install into a
# temporary prefix, find_package(SIMC) build of the example, a run that uses its plug-in (labels "examples;package").
add_test(NAME package_out_of_tree_plugin
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/out_of_tree_plugin/package_check.sh ${CMAKE_BINARY_DIR}
                 ${CMAKE_CURRENT_SOURCE_DIR}/out_of_tree_plugin ${CMAKE_CURRENT_BINARY_DIR}/out_of_tree_plugin ${_inp}
                 ${ROOT_DIR})
set_tests_properties(package_out_of_tree_plugin PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900 LABELS "examples;package")
