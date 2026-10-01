# tests/l1/radiation.cmake — task P5.1: radiative corrections (radc.f, brem.f, init.f radc_init*) against the P5.F1
# tables (reference/functions/*, reference/P5.F1_report.md). Included (OPTIONAL) by tests/l1/CMakeLists.txt.
# Deterministic tables: |rel| <= 1e-12 (D-010), NaN/Inf rows must match. Fixed-RNG tables: consumption counts exact.
add_executable(test_l1_radiation test_l1_radiation.cpp)
target_link_libraries(test_l1_radiation PRIVATE simc_l1runner simcfw_radiation simcfw_reactions simcfw_setup simcfw_material
                      simcfw_kinematics ROOT::MathCore)

set(_fnr ${SIMC_REFERENCE_DIR}/functions)
set(_rstate ${_fnr}/radc_state)
foreach(_t gamma spen spence schwinger lambda_dave brem bremos inter inter_prime srad extrad_friedrich extrad_phi
           radc_init radc_init_ev basicrad_init_ev peaked_rad_weight peaked_rad_weight_r1x3 basicrad_fixedrng)
    add_test(NAME l1_${_t} COMMAND test_l1_radiation ${_t} ${_fnr}/${_t} 1e-12)
    set_tests_properties(l1_${_t} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()

# generate_rad with scripted grnd/gauss1 and the full SIMC initialisation per configuration (22 .inp files);
# rows gen_src=1 are the real-stream replay of the first 1000 calls of the P1.2 rad_A run.
add_test(NAME l1_generate_rad_fixedrng
         COMMAND test_l1_radiation generate_rad_fixedrng ${_fnr}/generate_rad_fixedrng 1e-12 ${_rstate} ${SIMC_FORTRAN_DIR})
# every generate_rad call (21258) recorded in the capture run
add_test(NAME l1_generate_rad_capture
         COMMAND test_l1_radiation generate_rad_capture
                 ${SIMC_REFERENCE_DIR}/runs/p5f1_radA_capture/capture_generate_rad.txt 1e-12 ${_rstate} ${SIMC_FORTRAN_DIR})
# BuildRunSetup radiation quantities + radc_init_ev at the central vertex for the 22 configurations
add_test(NAME l1_radc_state COMMAND test_l1_radiation radc_state ${_rstate} 1e-12 ${SIMC_FORTRAN_DIR})
set_tests_properties(l1_generate_rad_fixedrng l1_generate_rad_capture l1_radc_state PROPERTIES SKIP_RETURN_CODE 77)

# real-RNG distributions (Mt64Random, 1e6 samples): quantiles, binned KS and chi2 (p >= 1e-3)
foreach(_id 1 2 3 4 5 6 7)
    add_test(NAME l1_rad_distributions_${_id}
             COMMAND test_l1_radiation rad_distributions ${_fnr}/rad_distributions ${_id} ${_rstate} ${SIMC_FORTRAN_DIR})
    set_tests_properties(l1_rad_distributions_${_id} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()

add_test(NAME radiation_selftest COMMAND test_l1_radiation selftest)
