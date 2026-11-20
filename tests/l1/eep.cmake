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

# U.1: upstream nuclear_elastic.f (sig_nuc_elastic, NUC_FORM_FACTOR, Fgauss) against the U.F1 tables, bitwise.
add_executable(test_l1_nuc_elastic test_l1_nuc_elastic.cpp)
target_link_libraries(test_l1_nuc_elastic PRIVATE simcfw_l1runner simcfw_physics simcfw_material simcfw_random)
add_test(NAME l1_nuc_elastic COMMAND test_l1_nuc_elastic ${SIMC_REFERENCE_DIR}/functions/nuclear_elastic)
set_tests_properties(l1_nuc_elastic PROPERTIES SKIP_RETURN_CODE 77 LABELS l1)

# U.1: full-run identity (RANLUX + --fortran-compat) against the upstream-SIMC references of task U.F1.
# A(e,e')A nuclear elastic (doing_nuc_elast): the 20d46fc build (last upstream commit with the complete path, finding U-A1);
# its ntuple is the 46-column (e,e'p) layout, i.e. the SIMC++ legacy layout (LAYOUT=legacy), compared bitwise in full.
foreach(_r "up_ee_calcium48_cap1000_20d46fc;upf1_ca48_cap" "up_ee_calcium48_20d46fc;ee_calcium48")
    list(GET _r 0 _dir)
    list(GET _r 1 _run)
    add_test(NAME l1_simcrun_ranlux_${_dir}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ranlux_fullrun.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${SIMC_REFERENCE_DIR}/runs/${_dir} ${_run} ${CMAKE_CURRENT_BINARY_DIR}/ranlux_${_dir})
    set_tests_properties(l1_simcrun_ranlux_${_dir} PROPERTIES LABELS "l1;replay" SKIP_RETURN_CODE 77 TIMEOUT 600
                         ENVIRONMENT "LAYOUT=legacy")
endforeach()
# U.2 (D-082): the upstream 59ccdf6 tree layout (the default) against the 59ccdf6 runs of task U.F1, the WHOLE h10 tree
# bitwise (every column), .hist/.geni/.gen text-identical; RANLUX + --fortran-compat. ranlux_fullrun.sh selects
# --tree-layout upstream from the reference .hist. Short runs (< 1 min) in the default set, the others "slow".
# U.R: H which_pion 3 (20d46fc: runs in compat mode now) and the kaon decay with Nntu = 0 (1fae2d7 transp.f:158 Mpi;
# .hist/.geni/.gen only, the kaon tree STOPs at 59ccdf6, U-A2). Not here: pion C-12 which_pion 2/3 (39d7de0 1/Z, U-A3,
# old weight kept pending JeffersonLab/simc_gfortran#75).
foreach(_r "up_upf1_rho_h_k1_rad_A_s;upf1_rho_h_k1_rad_A_s;fast"
           "up_upf1_pion_c12_wp0_A;upf1_pion_c12_wp0_A;fast"
           "up_eep_h_elastic_rad_A;eep_h_elastic_rad_A;slow" "up_eep_h_elastic_norad_A;eep_h_elastic_norad_A;slow"
           "up_eep_d_rad_A;eep_d_rad_A;slow" "up_eep_c12_rad_A;eep_c12_rad_A;slow"
           "up_upf1_eep_cu64_bh_rad_A;upf1_eep_cu64_bh_rad_A;slow"
           "up_pion_h_k1_rad_dec0_A;pion_h_k1_rad_dec0_A;slow" "up_delta_h_d1_rad_A;delta_h_d1_rad_A;slow"
           "up_p14f2_excl_h_ng2_A;p14f2_excl_h_ng2_A;slow" "up_sidis_h_k1_pip_rad_A;sidis_h_k1_pip_rad_A;slow"
           "up_u_nps_5p639;u_nps_5p639;slow"
           "up_upf1_pion_h_wp3_A;upf1_pion_h_wp3_A;slow"
           "up_upf1_kaon_h_k1_rad_dec1_nontu;upf1_kaon_h_k1_rad_dec1_nontu;slow")
    list(GET _r 0 _dir)
    list(GET _r 1 _run)
    list(GET _r 2 _speed)
    add_test(NAME l1_simcrun_ranlux_${_dir}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ranlux_fullrun.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${SIMC_REFERENCE_DIR}/runs/${_dir} ${_run} ${CMAKE_CURRENT_BINARY_DIR}/ranlux_${_dir})
    if(_speed STREQUAL "slow")
        set_tests_properties(l1_simcrun_ranlux_${_dir} PROPERTIES LABELS "l1;replay;upstream_layout;slow"
                             SKIP_RETURN_CODE 77 TIMEOUT 2400)
    else()
        set_tests_properties(l1_simcrun_ranlux_${_dir} PROPERTIES LABELS "l1;replay;upstream_layout"
                             SKIP_RETURN_CODE 77 TIMEOUT 600)
    endif()
endforeach()
