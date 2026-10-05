# tests/l1/eep.cmake — tasks P8.1/P8.2: D(e,e'p) and A(e,e'p) against the P8.F1 tables (reference/P8.F1_report.md).
# Included (OPTIONAL) by tests/l1/CMakeLists.txt. fortran_compat mode, bitwise (tolerance 0).
# The spectral-function files are the repository copies in data/spectral (byte-identical to the Fortran files,
# checked by l1_spectral_files).
add_executable(test_l1_eep test_l1_eep.cpp)
target_link_libraries(test_l1_eep PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})
target_compile_options(test_l1_eep PRIVATE -O2)

set(_fn ${SIMC_REFERENCE_DIR}/functions)
set(_sp ${CMAKE_SOURCE_DIR}/data/spectral)
set(_inp ${_fn}/eep_state/inp)
add_test(NAME l1_spectral_files COMMAND test_l1_eep files ${_sp} ${SIMC_FORTRAN_DIR})
add_test(NAME l1_theory_init COMMAND test_l1_eep theory_init ${_fn}/theory_init ${_inp} ${_sp} 0)
add_test(NAME l1_sf_lookup_init COMMAND test_l1_eep sf_lookup_init ${_fn}/sf_lookup_init ${_inp} ${_sp} 0)
add_test(NAME l1_sf_lookup COMMAND test_l1_eep sf_lookup ${_fn}/sf_lookup ${_sp} 0)
add_test(NAME l1_generate_em_fixedrng COMMAND test_l1_eep generate_em ${_fn}/generate_em_fixedrng ${_sp} 0)
set_tests_properties(l1_spectral_files l1_theory_init l1_sf_lookup_init l1_sf_lookup l1_generate_em_fixedrng
    PROPERTIES SKIP_RETURN_CODE 77)

# P8.2: the start-up state and the event routines of D(e,e'p)/A(e,e'p)
add_test(NAME l1_eep_state COMMAND test_l1_eep eep_state ${_fn}/eep_state ${_inp} ${_sp} 0)
add_test(NAME l1_generate_eepA_fixedrng COMMAND test_l1_eep generate ${_fn}/generate_eepA_fixedrng ${_inp} ${_sp} 0)
add_test(NAME l1_complete_ev_eepA COMMAND test_l1_eep complete_ev ${_fn}/complete_ev_eepA ${_inp} ${_sp} 0)
add_test(NAME l1_complete_recon_ev_eepA COMMAND test_l1_eep complete_recon_ev ${_fn}/complete_recon_ev_eepA ${_inp} ${_sp} 0)
add_test(NAME l1_complete_main_eepA COMMAND test_l1_eep complete_main ${_fn}/complete_main_eepA ${_inp} ${_sp} 0)
set_tests_properties(l1_eep_state l1_generate_eepA_fixedrng l1_complete_ev_eepA l1_complete_recon_ev_eepA
    l1_complete_main_eepA PROPERTIES SKIP_RETURN_CODE 77)
# replay of every captured real call (generate with the captured RANLUX stream and persistent state, complete_ev,
# complete_recon_ev, complete_main) of the D and Fe-56 capture runs
foreach(_c "d;p8f1_d_rad_cap" "fe56;p8f1_fe56_rad_cap")
    list(GET _c 0 _tag)
    list(GET _c 1 _run)
    add_test(NAME l1_event_capture_eepA_${_tag}
             COMMAND test_l1_eep capture ${_fn}/event_capture_eepA/cap_${_tag} ${SIMC_REFERENCE_DIR}/runs/${_run}/${_run}.inp
                     ${_sp} 0)
    set_tests_properties(l1_event_capture_eepA_${_tag} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900)
endforeach()
