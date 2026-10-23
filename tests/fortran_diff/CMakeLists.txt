# tests/fortran_diff — task P7.5 (D-031, D-032): the documented differences with the Fortran (docs/FORTRAN_DIFFERENCES.md).
# For every corrected clear Fortran bug (class a) one test shows BOTH behaviours: fortran_compat = true gives the
# Fortran result, the default (corrected) gives the corrected one. Label "fortran_diff" (+ "unit"/"cli").
add_executable(test_fortran_differences test_fortran_differences.cpp ${CMAKE_SOURCE_DIR}/tests/l1/L1Runner.cpp)
target_include_directories(test_fortran_differences PRIVATE ${CMAKE_SOURCE_DIR}/tests/l1)
target_link_libraries(test_fortran_differences PRIVATE ${SIMC_DRIVER_LIBS})
target_compile_definitions(test_fortran_differences PRIVATE SIMC_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# One CTest per item id (the ids of docs/FORTRAN_DIFFERENCES.md); the test program takes the id as argument.
set(SIMC_FORTRAN_DIFF_ITEMS
    D006_pass_cuts
    D006_kaon_decay_mass
    D011_D018_report_formats
    D014_deforest_cos_phi_warning
    D015_rotate_raddeg
    D017_extreme_trip_hadron_angle
    D021_hms_coll_slit_counters
    D028_spectrometer_angles_nan
    D032_spectrometer_id
    D033_F3_benhar_sf_diff
    D033_F5_pmper_stale_pm
    D033_F6_sf_lookup_stale
    D033_F7_heavy_mrec_guard
    D037_F1_shms_coll_slit_counters
    D037_F2_coll_stops_geni
    D037_F4_coll_ranged_out
    D037_F7_project_decay_drift
    D037_F9_decay_non_pi_k
    D038_F1_delta_pp_stop
    D038_F2_pi0_neutron_on_h
    D038_F5_central_sigcc
    D038_F15_hist_which_pion_11
    D043_F13_central_sigcc_kaon
    D048_F1_delta_DA_stop
    D048_F2_delta_tail1_egamma
    D048_F3_delta_sigblok_args
    D048_F6_delta_warning_text
    D048_F9_central_sigcc_delta
    D051_F3_F9_semi_pfermi_pmper
    D051_F12_F14_semi_hist_titles
)
foreach(_item ${SIMC_FORTRAN_DIFF_ITEMS})
    add_test(NAME fortran_diff_${_item}
             COMMAND test_fortran_differences ${_item} ${SIMC_REFERENCE_DIR} ${SIMC_FORTRAN_DIR})
    set_tests_properties(fortran_diff_${_item} PROPERTIES SKIP_RETURN_CODE 77 LABELS "fortran_diff;unit")
endforeach()

add_test(NAME cli_arm_code_6_rejected
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cli_arm_code_rejected.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:inp2json>
                 ${CMAKE_CURRENT_BINARY_DIR}/cli_arm_code_rejected)
set_tests_properties(cli_arm_code_6_rejected PROPERTIES LABELS "fortran_diff;cli")

# D-037 F9 (P6.4): doing_decay with a non-pi/K hadron: corrected = rejected at start-up (exit 2), compat = FortranStop (3)
add_test(NAME cli_decay_non_pik
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cli_decay_non_pik.sh $<TARGET_FILE:simcpp> ${SIMC_REFERENCE_DIR}
                 ${CMAKE_CURRENT_BINARY_DIR}/cli_decay_non_pik)
set_tests_properties(cli_decay_non_pik PROPERTIES LABELS "fortran_diff;cli" SKIP_RETURN_CODE 77)
