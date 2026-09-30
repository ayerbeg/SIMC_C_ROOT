# tests/unit/setup.cmake — task P2.1: RunSetup (BuildRunSetup) vs the Fortran .hist echo of every
# derived quantity. Included from tests/unit/CMakeLists.txt.
add_executable(test_setup_hist ${CMAKE_CURRENT_SOURCE_DIR}/test_setup_hist.cpp)
target_link_libraries(test_setup_hist PRIVATE simcfw_setup)

# simc_setup_hist(<name> <run dir/stem> <ctp_float 0|1> [cfg]) — with cfg, the energy-loss numbers come
# from the Fortran harness tables target_init/ and extreme_trip_thru_target/ row cfg. ctp_float=1 rounds
# the real inputs to float as CTP does (CTP/thTestParse.c:391-394); see test_setup_hist.cpp.
function(simc_setup_hist name stem ctp_float)
    set(_args ${SIMC_REFERENCE_DIR}/runs/${stem}.inp ${SIMC_REFERENCE_DIR}/runs/${stem}.hist
              --data-dir ${SIMC_FORTRAN_DIR} --ctp-float ${ctp_float})
    if(ARGC GREATER 3)
        list(APPEND _args --tables ${SIMC_REFERENCE_DIR}/functions --cfg ${ARGV3})
    endif()
    add_test(NAME setup_hist_${name} COMMAND test_setup_hist ${_args})
    set_tests_properties(setup_hist_${name} PROPERTIES SKIP_RETURN_CODE 77)
endfunction()

simc_setup_hist(norad_A eep_h_elastic_norad_A/eep_h_elastic_norad_A 1)
simc_setup_hist(norad_A_tables eep_h_elastic_norad_A/eep_h_elastic_norad_A 1 1)
simc_setup_hist(rad_A eep_h_elastic_rad_A/eep_h_elastic_rad_A 1)
simc_setup_hist(rad_A_tables eep_h_elastic_rad_A/eep_h_elastic_rad_A 1 1)
simc_setup_hist(acc_phsp_A acc_phsp_hms_shms_A/acc_phsp_hms_shms_A 1)
simc_setup_hist(acc_phsp_A_tables acc_phsp_hms_shms_A/acc_phsp_hms_shms_A 1 1)
simc_setup_hist(eep_h_test_pctnames eep_h_test_pctnames/eep_h_test_pctnames 1)
simc_setup_hist(hms_shms_validation hms_shms_validation/test_hms_shms_validation 1)
# KNOWN FAILURE (xfail, reported in P2.1): the same chain with the InpReader values as they are (full
# double). It fails (spec%e%P printed 8360.9890 vs Fortran 8360.9893; edge%e%E%max 10243.610 vs 10243.611)
# until src/input/InpReader reproduces the CTP float literals. WILL_FAIL makes CTest report it as passed
# while it fails; once the InpReader is fixed it turns red: then drop WILL_FAIL and --ctp-float.
simc_setup_hist(rad_A_tables_inpreader_double_XFAIL eep_h_elastic_rad_A/eep_h_elastic_rad_A 0 1)
set_tests_properties(setup_hist_rad_A_tables_inpreader_double_XFAIL PROPERTIES WILL_FAIL TRUE)
