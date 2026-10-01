# tests/l1/driver.cmake — task P7.1: simc.f:montecarlo (simc::driver::MonteCarlo) replayed against the real-call
# records of the P4.F1 captures (reference/functions/event_capture/<cfg>, reference/runs/p4f1_<cfg>_cap).
# Included (OPTIONAL) by tests/l1/CMakeLists.txt. Everything bitwise; needs `xz` for the .txt.xz capture files.
add_executable(test_l1_driver test_l1_driver.cpp)
target_link_libraries(test_l1_driver PRIVATE simcfw_driver simcfw_setup simcfw_hms simcfw_shms simcfw_optics
    simcfw_material)
target_compile_options(test_l1_driver PRIVATE -O2)
foreach(_c norad rad phsp phsp_tries)
    add_test(NAME l1_montecarlo_replay_${_c} COMMAND test_l1_driver replay
        ${SIMC_REFERENCE_DIR}/functions/event_capture/${_c} ${SIMC_REFERENCE_DIR}/runs/p4f1_${_c}_cap ${_c}
        ${SIMC_FORTRAN_DIR})
    set_tests_properties(l1_montecarlo_replay_${_c} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)
endforeach()
