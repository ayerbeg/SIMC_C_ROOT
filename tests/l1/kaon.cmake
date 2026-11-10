# tests/l1/kaon.cmake — task P10.1: kaon cross sections (physics_kaon.f, cern/fint.f, dbase.f:624-679) against the P10.F1
# tables (reference/P10.F1_report.md). Included (OPTIONAL) by tests/l1/CMakeLists.txt after its label loop. fortran_compat
# mode, BITWISE (tolerance 0) unless stated; the *_fintabi tables are the reference of the corrected fint path (D-043 F1).
add_executable(test_l1_kaon test_l1_kaon.cpp)
target_link_libraries(test_l1_kaon PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})

set(_kfn ${SIMC_REFERENCE_DIR}/functions)
set(_kdata ${CMAKE_SOURCE_DIR}/data/kaon)
add_test(NAME l1_kaon_data_files COMMAND test_l1_kaon files ${_kdata} ${SIMC_FORTRAN_DIR})
add_test(NAME l1_saghai_tables COMMAND test_l1_kaon saghai_tables ${_kfn}/saghai_tables ${_kdata})
add_test(NAME l1_sig_factorized COMMAND test_l1_kaon sig_factorized ${_kfn}/sig_factorized 0)
add_test(NAME l1_eekeek COMMAND test_l1_kaon eekeek ${_kfn}/eekeek ${_kdata} 0 0)
add_test(NAME l1_eekeek_fintabi COMMAND test_l1_kaon eekeek ${_kfn}/eekeek_fintabi ${_kdata} 0 1)
add_test(NAME l1_eekeeks COMMAND test_l1_kaon eekeeks ${_kfn}/eekeeks ${_kdata} 0 0)
add_test(NAME l1_eekeeks_fintabi COMMAND test_l1_kaon eekeeks ${_kfn}/eekeeks_fintabi ${_kdata} 0 1)
add_test(NAME l1_phspwght COMMAND test_l1_kaon phspwght ${_kfn}/phspwght ${_kdata} 0 0)
add_test(NAME l1_phspwght_fintabi COMMAND test_l1_kaon phspwght ${_kfn}/phspwght_fintabi ${_kdata} 0 1)
add_test(NAME l1_peek COMMAND test_l1_kaon peek ${_kfn}/peek ${_kfn}/kaon_state ${_kdata} 0 0)
add_test(NAME l1_peek_fintabi COMMAND test_l1_kaon peek ${_kfn}/peek_fintabi ${_kfn}/kaon_state ${_kdata} 0 1)
add_test(NAME l1_kaon_plugin COMMAND test_l1_kaon plugin ${_kfn}/peek ${_kfn}/kaon_state ${_kdata})
add_test(NAME l1_transform_to_cm_kaon COMMAND test_l1_pion transform_to_cm ${_kfn}/transform_to_cm_kaon 0)
set(_kaon_l1 l1_kaon_data_files l1_saghai_tables l1_sig_factorized l1_eekeek l1_eekeek_fintabi l1_eekeeks
    l1_eekeeks_fintabi l1_phspwght l1_phspwght_fintabi l1_peek l1_peek_fintabi l1_kaon_plugin l1_transform_to_cm_kaon)
set_tests_properties(${_kaon_l1} PROPERTIES SKIP_RETURN_CODE 77 LABELS l1)
# D-043 F1/F2/F3/F5/F6 (FORTRAN-BUG, corrected by default): both modes
add_test(NAME fortran_diff_D043_kaon_model
         COMMAND test_l1_kaon fortran_diff_D043 ${_kfn}/peek ${_kfn}/kaon_state ${_kdata} ${_kfn}/eekeek ${_kfn}/saghai_tables)
set_tests_properties(fortran_diff_D043_kaon_model PROPERTIES SKIP_RETURN_CODE 77 LABELS "fortran_diff;unit")

# task P10.2: the kaon reaction (reactions::KaonProduction + "kaon" model) against the P10.F1 event tables, fortran_compat,
# bitwise; the same runner as the pion tables (test_l1_pion_event ... kaon <data/kaon>)
set(_kinp ${_kfn}/kaon_state/inp)
set(_ksp ${CMAKE_SOURCE_DIR}/data/spectral)
set(_kpd ${CMAKE_SOURCE_DIR}/data/pion)
add_test(NAME l1_kaon_state COMMAND test_l1_pion_event pion_state ${_kfn}/kaon_state ${_kinp} ${_ksp} ${_kpd} 0 kaon ${_kdata})
add_test(NAME l1_generate_kaon_fixedrng
         COMMAND test_l1_pion_event generate ${_kfn}/generate_kaon_fixedrng ${_kinp} ${_ksp} ${_kpd} 0 kaon ${_kdata})
add_test(NAME l1_complete_ev_kaon
         COMMAND test_l1_pion_event complete_ev ${_kfn}/complete_ev_kaon ${_kinp} ${_ksp} ${_kpd} 0 kaon ${_kdata})
add_test(NAME l1_complete_recon_ev_kaon
         COMMAND test_l1_pion_event complete_recon_ev ${_kfn}/complete_recon_ev_kaon ${_kinp} ${_ksp} ${_kpd} 0 kaon ${_kdata})
add_test(NAME l1_complete_main_kaon
         COMMAND test_l1_pion_event complete_main ${_kfn}/complete_main_kaon ${_kinp} ${_ksp} ${_kpd} 0 kaon ${_kdata})
set(_kaon_ev l1_kaon_state l1_generate_kaon_fixedrng l1_complete_ev_kaon l1_complete_recon_ev_kaon l1_complete_main_kaon)
foreach(_c "kl_h;p10f1_kl_h_dec1_cap" "ksm_d;p10f1_ksm_d_dec1_cap")
    list(GET _c 0 _tag)
    list(GET _c 1 _run)
    add_test(NAME l1_replay_event_capture_kaon_${_tag}
             COMMAND test_l1_pion_event capture ${_kfn}/event_capture_kaon/cap_${_tag}
                     ${SIMC_REFERENCE_DIR}/runs/${_run}/${_run}.inp ${_ksp} ${_kpd} 0 kaon ${_kdata})
    list(APPEND _kaon_ev l1_replay_event_capture_kaon_${_tag})
endforeach()
set_tests_properties(${_kaon_ev} PROPERTIES SKIP_RETURN_CODE 77 LABELS l1 TIMEOUT 900)
set_tests_properties(l1_replay_event_capture_kaon_kl_h l1_replay_event_capture_kaon_ksm_d PROPERTIES LABELS "l1;replay")
# the kaon project/transp decay calls and mc_hms/mc_shms with kaon decay (same runner as the P9.F1 pion tables)
add_test(NAME l1_decay_kaon_fixedrng COMMAND test_l1_coll_decay pdecay ${_kfn}/decay_kaon_fixedrng)
add_test(NAME l1_mc_decay_kaon_fixedrng COMMAND test_l1_coll_decay mcdecay ${_kfn}/mc_decay_kaon_fixedrng)
set_tests_properties(l1_decay_kaon_fixedrng l1_mc_decay_kaon_fixedrng PROPERTIES SKIP_RETURN_CODE 77 LABELS l1)

# task P10.2: the full-run replay of the two P10.F1 kaon captures (K+ Lambda H and K+ Sigma- D, K1, rad, decay ON in the
# SHMS; reference/functions/event_capture_kaon/cap_{kl_h,ksm_d}, reference/runs/p10f1_{kl_h,ksm_d}_dec1_cap): 55-column
# tree, .bin, bookkeeping, normalisation, .hist/.geni/.gen text; and every montecarlo call (kaon decay draws).
foreach(_c "kl_h;p10f1_kl_h_dec1_cap" "ksm_d;p10f1_ksm_d_dec1_cap")
    list(GET _c 0 _tag)
    list(GET _c 1 _run)
    set(_out ${CMAKE_CURRENT_BINARY_DIR}/simcrun_replay_kaon_${_tag})
    file(MAKE_DIRECTORY ${_out})
    add_test(NAME l1_simcrun_replay_kaon_${_tag} COMMAND test_l1_simcrun replay
        ${_kfn}/event_capture_kaon/cap_${_tag} ${SIMC_REFERENCE_DIR}/runs/${_run} ${_tag} ${_out} ${_run})
    add_test(NAME l1_montecarlo_replay_kaon_${_tag} COMMAND test_l1_driver replay
        ${_kfn}/event_capture_kaon/cap_${_tag} ${SIMC_REFERENCE_DIR}/runs/${_run} ${_tag} ${_ksp} ${_run})
    set_tests_properties(l1_simcrun_replay_kaon_${_tag} l1_montecarlo_replay_kaon_${_tag}
        PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900 LABELS "l1;replay")
endforeach()
# P15.1 (D-058): the P10 L2 flag on kaon_h_k1_norad_dec0_A closed by a same-seed run: RANLUX + --fortran-compat reproduce the
# 200000-event Fortran run exactly (.hist except the time lines, .geni, .gen, h10 tree)
add_test(NAME l1_simcrun_ranlux_kaon_h_k1_norad_dec0_A
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ranlux_fullrun.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                 ${SIMC_REFERENCE_DIR}/runs/kaon_h_k1_norad_dec0_A kaon_h_k1_norad_dec0_A
                 ${CMAKE_CURRENT_BINARY_DIR}/ranlux_kaon_h_k1_norad_dec0_A)
set_tests_properties(l1_simcrun_ranlux_kaon_h_k1_norad_dec0_A PROPERTIES LABELS "l1;replay;slow" SKIP_RETURN_CODE 77 TIMEOUT 1800)
