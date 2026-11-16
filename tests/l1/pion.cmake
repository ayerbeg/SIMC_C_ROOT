# tests/l1/pion.cmake — task P9.1: exclusive pion cross sections (physics_pion.f), transform_to_cm on pion kinematics
# (jacobians.f) and pizero_decay.f against the P9.F1 tables (reference/P9.F1_report.md). Included (OPTIONAL) by
# tests/l1/CMakeLists.txt AFTER its label loop, so the labels are set here. fortran_compat mode, |rel| <= 1e-12 (the
# reports give the bitwise counts). The MAID tables and MULTIPOL.file are the repository copies in data/pion
# (byte-identical to the Fortran files, checked by l1_pion_data_files).
add_executable(test_l1_pion test_l1_pion.cpp)
target_link_libraries(test_l1_pion PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})

set(_pfn ${SIMC_REFERENCE_DIR}/functions)
set(_pdata ${CMAKE_SOURCE_DIR}/data/pion)
set(_ptol 1e-12)
add_test(NAME l1_pion_data_files COMMAND test_l1_pion files ${_pdata} ${SIMC_FORTRAN_DIR})
foreach(_r sig_blok sig_param04 sig_param_3000 sig_param_2021 exclfit)
    add_test(NAME l1_${_r} COMMAND test_l1_pion ${_r} ${_pfn}/${_r} ${_ptol})
endforeach()
add_test(NAME l1_sigmaid COMMAND test_l1_pion sigmaid ${_pfn}/sigmaid ${_pdata} ${_ptol})
add_test(NAME l1_sig_multipole COMMAND test_l1_pion sig_multipole ${_pfn}/sig_multipole ${_pdata} ${_ptol})
add_test(NAME l1_transform_to_cm_pion COMMAND test_l1_pion transform_to_cm ${_pfn}/transform_to_cm_pion ${_ptol})
add_test(NAME l1_peepi COMMAND test_l1_pion peepi ${_pfn}/peepi ${_pfn}/pion_state ${_pdata} ${_ptol})
add_test(NAME l1_pizero_decay_fixedrng COMMAND test_l1_pion pizero_decay ${_pfn}/pizero_decay_fixedrng ${_ptol})
add_test(NAME l1_pion_plugin COMMAND test_l1_pion plugin ${_pfn}/peepi ${_pfn}/pion_state ${_pdata})
set(_pion_l1 l1_pion_data_files l1_sig_blok l1_sig_param04 l1_sig_param_3000 l1_sig_param_2021 l1_exclfit l1_sigmaid
    l1_sig_multipole l1_transform_to_cm_pion l1_peepi l1_pizero_decay_fixedrng l1_pion_plugin)
set_tests_properties(${_pion_l1} PROPERTIES SKIP_RETURN_CODE 77 LABELS l1)

# D-038 F10 (FORTRAN-BUG, corrected by default): both modes on the peepi table rows (docs/translation/markdown/FORTRAN_DIFFERENCES.md)
add_test(NAME fortran_diff_D038_F10_sigcm2_stale
         COMMAND test_l1_pion fortran_diff_D038_F10 ${_pfn}/peepi ${_pfn}/pion_state ${_pdata})
set_tests_properties(fortran_diff_D038_F10_sigcm2_stale PROPERTIES SKIP_RETURN_CODE 77 LABELS "fortran_diff;unit")

# task P9.2: the pion reaction (reactions::PionProduction) against the P9.F1 event tables, fortran_compat, bitwise
add_executable(test_l1_pion_event test_l1_pion_event.cpp)
target_link_libraries(test_l1_pion_event PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})
target_compile_options(test_l1_pion_event PRIVATE -O2)
set(_pinp ${_pfn}/pion_state/inp)
set(_psp ${CMAKE_SOURCE_DIR}/data/spectral)
add_test(NAME l1_pion_state COMMAND test_l1_pion_event pion_state ${_pfn}/pion_state ${_pinp} ${_psp} ${_pdata} 0)
add_test(NAME l1_generate_pion_fixedrng
         COMMAND test_l1_pion_event generate ${_pfn}/generate_pion_fixedrng ${_pinp} ${_psp} ${_pdata} 0)
add_test(NAME l1_complete_ev_pion COMMAND test_l1_pion_event complete_ev ${_pfn}/complete_ev_pion ${_pinp} ${_psp} ${_pdata} 0)
add_test(NAME l1_complete_recon_ev_pion
         COMMAND test_l1_pion_event complete_recon_ev ${_pfn}/complete_recon_ev_pion ${_pinp} ${_psp} ${_pdata} 0)
add_test(NAME l1_complete_main_pion
         COMMAND test_l1_pion_event complete_main ${_pfn}/complete_main_pion ${_pinp} ${_psp} ${_pdata} 0)
set(_pion_ev l1_pion_state l1_generate_pion_fixedrng l1_complete_ev_pion l1_complete_recon_ev_pion l1_complete_main_pion)
foreach(_c "pip_h;p9f1_pip_h_dec1_cap" "pim_d;p9f1_pim_d_dec1_cap")
    list(GET _c 0 _tag)
    list(GET _c 1 _run)
    add_test(NAME l1_replay_event_capture_pion_${_tag}
             COMMAND test_l1_pion_event capture ${_pfn}/event_capture_pion/cap_${_tag}
                     ${SIMC_REFERENCE_DIR}/runs/${_run}/${_run}.inp ${_psp} ${_pdata} 0)
    list(APPEND _pion_ev l1_replay_event_capture_pion_${_tag})
endforeach()
set_tests_properties(${_pion_ev} PROPERTIES SKIP_RETURN_CODE 77 LABELS l1 TIMEOUT 900)
set_tests_properties(l1_replay_event_capture_pion_pip_h l1_replay_event_capture_pion_pim_d PROPERTIES LABELS "l1;replay")
