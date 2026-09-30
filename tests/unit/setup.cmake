# tests/unit/setup.cmake — tasks P2.1/P2.2: RunSetup (BuildRunSetup with the real HallCCryoTarget provider,
# InpReader with CTP float literals) vs the Fortran .hist echo of every derived quantity.
# Included from tests/unit/CMakeLists.txt.
add_executable(test_setup_hist ${CMAKE_CURRENT_SOURCE_DIR}/test_setup_hist.cpp)
target_link_libraries(test_setup_hist PRIVATE simcfw_setup)

function(simc_setup_hist name stem)
    add_test(NAME setup_hist_${name} COMMAND test_setup_hist ${SIMC_REFERENCE_DIR}/runs/${stem}.inp
             ${SIMC_REFERENCE_DIR}/runs/${stem}.hist --data-dir ${SIMC_FORTRAN_DIR})
    set_tests_properties(setup_hist_${name} PROPERTIES SKIP_RETURN_CODE 77)
endfunction()

simc_setup_hist(norad_A eep_h_elastic_norad_A/eep_h_elastic_norad_A)
simc_setup_hist(rad_A eep_h_elastic_rad_A/eep_h_elastic_rad_A)
simc_setup_hist(acc_phsp_A acc_phsp_hms_shms_A/acc_phsp_hms_shms_A)
simc_setup_hist(eep_h_test_pctnames eep_h_test_pctnames/eep_h_test_pctnames)
simc_setup_hist(hms_shms_validation hms_shms_validation/test_hms_shms_validation)
