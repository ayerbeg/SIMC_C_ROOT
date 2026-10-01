# tests/l1/reactions.cmake — task P4.1: event generation (event.f generate, complete_ev, complete_recon_ev,
# complete_main; jacobians.f transform_to_cm) against the P4.F1 tables and capture runs
# (reference/P4.F1_report.md). Included (OPTIONAL) by tests/l1/CMakeLists.txt. Tolerance 1e-12 (bitwise expected).
add_executable(test_l1_reactions test_l1_reactions.cpp)
target_link_libraries(test_l1_reactions PRIVATE simc_l1runner simcfw_reactions simcfw_radiation simcfw_setup
                      simcfw_material simcfw_kinematics simcfw_physics)

set(_fne ${SIMC_REFERENCE_DIR}/functions)
set(_runs ${SIMC_REFERENCE_DIR}/runs)
foreach(_t generate_fixedrng complete_ev complete_recon_ev complete_main)
    add_test(NAME l1_${_t} COMMAND test_l1_reactions table ${_t} ${_fne}/${_t} 1e-12 ${_runs} ${SIMC_FORTRAN_DIR})
    set_tests_properties(l1_${_t} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()
add_test(NAME l1_transform_to_cm COMMAND test_l1_reactions transform_to_cm ${_fne}/transform_to_cm 1e-12)
set_tests_properties(l1_transform_to_cm PROPERTIES SKIP_RETURN_CODE 77)

# replay of every captured real call (generate with the captured RANLUX stream, complete_ev, complete_recon_ev,
# complete_main) for the four capture runs
foreach(_c norad rad phsp phsp_tries)
    add_test(NAME l1_event_capture_${_c}
             COMMAND test_l1_reactions capture ${_fne}/event_capture/${_c} ${_runs}/p4f1_${_c}_cap/p4f1_${_c}_cap.inp
                     1e-12 ${SIMC_FORTRAN_DIR})
    set_tests_properties(l1_event_capture_${_c} PROPERTIES SKIP_RETURN_CODE 77 LABELS slow)
endforeach()

add_test(NAME reactions_selftest COMMAND test_l1_reactions selftest)
