# tests/l1/setup.cmake — task P2.1 (run setup + kinematics utilities). Included from tests/l1/CMakeLists.txt.
set(_fn ${SIMC_REFERENCE_DIR}/functions)
set(_runs ${SIMC_REFERENCE_DIR}/runs)

# L1: physics_angles, spectrometer_angles, loren (|rel| <= 1e-12), target_init non-eloss columns.
add_executable(test_l1_setup ${CMAKE_CURRENT_SOURCE_DIR}/test_l1_setup.cpp)
target_link_libraries(test_l1_setup PRIVATE simcfw_l1runner simcfw_setup)
foreach(r physics_angles spectrometer_angles loren)
    add_test(NAME l1_${r} COMMAND test_l1_setup ${r} ${_fn}/${r} 1e-12)
    set_tests_properties(l1_${r} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()
add_test(NAME l1_target_init_nonloss COMMAND test_l1_setup target_init_nonloss ${_fn}/target_init 1e-12 ${SIMC_FORTRAN_DIR})
set_tests_properties(l1_target_init_nonloss PROPERTIES SKIP_RETURN_CODE 77)
