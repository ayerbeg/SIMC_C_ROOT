# tests/l1/calo.cmake — task P14.3: the calorimeter arm (ISpectrometer "calo", calo/mc_calo.f) and pi0 -> 2 gamma
# (driver::PizeroCalorimeter, pizero_decay) against the P14.F2 tables reference/functions/calo/ (fortran_compat, bitwise).
# Included (OPTIONAL) by tests/l1/CMakeLists.txt after its label loop: labels are set here.
add_executable(test_l1_calo test_l1_calo.cpp)
target_link_libraries(test_l1_calo PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})
target_compile_options(test_l1_calo PRIVATE -O2)
set(_cfn ${SIMC_REFERENCE_DIR}/functions/calo)
add_test(NAME calo_selftest COMMAND test_l1_calo selftest)
add_test(NAME l1_mc_calo COMMAND test_l1_calo mc ${_cfn}/mc_calo 0)
add_test(NAME l1_mc_calo_stops COMMAND test_l1_calo stops ${_cfn}/mc_calo_stops)
add_test(NAME l1_mc_calo_recon COMMAND test_l1_calo recon ${_cfn}/mc_calo_recon 0)
add_test(NAME l1_calo_frame_rotation COMMAND test_l1_calo rotation ${_cfn}/calo_frame_rotation)
add_test(NAME l1_pizero_decay_nps_fixedrng COMMAND test_l1_pion pizero_decay ${_cfn}/pizero_decay_nps_fixedrng 0)
# the P9.F1 event drivers on the NPS pi0 configurations (cfg -> reference/runs/<run>/<run>.inp)
set(_cruns ${SIMC_REFERENCE_DIR}/runs)
set(_csp ${CMAKE_SOURCE_DIR}/data/spectral)
set(_cpd ${CMAKE_SOURCE_DIR}/data/pion)
add_test(NAME l1_pi0calo_state COMMAND test_l1_pion_event pion_state ${_cfn}/pi0calo_state ${_cruns} ${_csp} ${_cpd} 0 calo)
foreach(_t generate:generate_fixedrng_pi0calo complete_ev:complete_ev_pi0calo complete_recon_ev:complete_recon_ev_pi0calo
        complete_main:complete_main_pi0calo)
    string(REPLACE ":" ";" _p ${_t})
    list(GET _p 0 _mode)
    list(GET _p 1 _tab)
    add_test(NAME l1_${_tab} COMMAND test_l1_pion_event ${_mode} ${_cfn}/${_tab} ${_cruns} ${_csp} ${_cpd} 0 calo)
    list(APPEND _calo_ev l1_${_tab})
endforeach()
set_tests_properties(calo_selftest PROPERTIES LABELS "l1;unit")
set_tests_properties(l1_mc_calo l1_mc_calo_stops l1_mc_calo_recon l1_calo_frame_rotation l1_pizero_decay_nps_fixedrng
    l1_pi0calo_state ${_calo_ev} PROPERTIES LABELS l1 SKIP_RETURN_CODE 77 TIMEOUT 900)
# full simc++ runs of the two P14.F2 captures (1000 events each; excl H pi0 ng2 and H(e,e'pi0)X) with RANLUX + compat:
# .hist/.geni/.gen and the h10 tree identical to the Fortran (the strongest form of the capture replay)
foreach(_r p14f2_excl_h_ng2_cap p14f2_sidis_h_cap)
    add_test(NAME l1_simcrun_ranlux_${_r}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ranlux_fullrun.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${SIMC_REFERENCE_DIR}/runs/${_r} ${_r} ${CMAKE_CURRENT_BINARY_DIR}/ranlux_${_r})
    set_tests_properties(l1_simcrun_ranlux_${_r} PROPERTIES LABELS "l1;replay" SKIP_RETURN_CODE 77 TIMEOUT 900)
endforeach()
# the 100k P14.F2 L2 run with the calorimeter as electron arm (1.8e6 tries)
add_test(NAME l1_simcrun_ranlux_p14f2_eep_caloe_A
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ranlux_fullrun.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                 ${SIMC_REFERENCE_DIR}/runs/p14f2_eep_caloe_A p14f2_eep_caloe_A ${CMAKE_CURRENT_BINARY_DIR}/ranlux_p14f2_eep_caloe_A)
set_tests_properties(l1_simcrun_ranlux_p14f2_eep_caloe_A PROPERTIES LABELS "l1;replay;slow" SKIP_RETURN_CODE 77 TIMEOUT 1800)
