# tests/l1/delta.cmake — task P11.1: the Delta region (e,e'p)pi (reactions::DeltaProduction + the "delta" model,
# physics_delta.f:peedelta) against the P11.F1 tables (reference/P11.F1_report.md, reference/functions/README.md "Delta
# region"). fortran_compat mode, BITWISE (tolerance 0). Same runners as the pion/kaon tables (test_l1_pion_event ... delta).
# Included (OPTIONAL) by tests/l1/CMakeLists.txt after the pion/kaon files.
set(_dfn ${SIMC_REFERENCE_DIR}/functions)
set(_dinp ${_dfn}/delta_state/inp)
set(_dsp ${CMAKE_SOURCE_DIR}/data/spectral)
set(_dpd ${CMAKE_SOURCE_DIR}/data/pion)
add_test(NAME l1_delta_state COMMAND test_l1_pion_event pion_state ${_dfn}/delta_state ${_dinp} ${_dsp} ${_dpd} 0 delta)
add_test(NAME l1_generate_delta_fixedrng
         COMMAND test_l1_pion_event generate ${_dfn}/generate_delta_fixedrng ${_dinp} ${_dsp} ${_dpd} 0 delta)
add_test(NAME l1_complete_ev_delta
         COMMAND test_l1_pion_event complete_ev ${_dfn}/complete_ev_delta ${_dinp} ${_dsp} ${_dpd} 0 delta)
add_test(NAME l1_complete_recon_ev_delta
         COMMAND test_l1_pion_event complete_recon_ev ${_dfn}/complete_recon_ev_delta ${_dinp} ${_dsp} ${_dpd} 0 delta)
add_test(NAME l1_complete_main_delta
         COMMAND test_l1_pion_event complete_main ${_dfn}/complete_main_delta ${_dinp} ${_dsp} ${_dpd} 0 delta)
add_test(NAME l1_peedelta COMMAND test_l1_pion_event peedelta ${_dfn}/peedelta ${_dinp} ${_dsp} ${_dpd} 0 delta)
add_test(NAME l1_transform_to_cm_delta COMMAND test_l1_pion transform_to_cm ${_dfn}/transform_to_cm_delta 0)
set(_delta_l1 l1_delta_state l1_generate_delta_fixedrng l1_complete_ev_delta l1_complete_recon_ev_delta
    l1_complete_main_delta l1_peedelta l1_transform_to_cm_delta)
# the two P11.F1 captures (H(e,e'p)pi rad, D1 = SHMS e- + HMS p, D2 = arms swapped): every generate / complete_ev /
# complete_recon_ev / complete_main call, the full SimcRun (53-column tree, .bin values, bookkeeping, .hist/.geni/.gen)
# and every montecarlo call
foreach(_c "d1;p11f1_dh_d1_rad_cap" "d2;p11f1_dh_d2_rad_cap")
    list(GET _c 0 _tag)
    list(GET _c 1 _run)
    add_test(NAME l1_replay_event_capture_delta_${_tag}
             COMMAND test_l1_pion_event capture ${_dfn}/event_capture_delta/cap_${_tag}
                     ${SIMC_REFERENCE_DIR}/runs/${_run}/${_run}.inp ${_dsp} ${_dpd} 0 delta)
    set(_out ${CMAKE_CURRENT_BINARY_DIR}/simcrun_replay_delta_${_tag})
    file(MAKE_DIRECTORY ${_out})
    add_test(NAME l1_simcrun_replay_delta_${_tag} COMMAND test_l1_simcrun replay
        ${_dfn}/event_capture_delta/cap_${_tag} ${SIMC_REFERENCE_DIR}/runs/${_run} ${_tag} ${_out} ${_run})
    add_test(NAME l1_montecarlo_replay_delta_${_tag} COMMAND test_l1_driver replay
        ${_dfn}/event_capture_delta/cap_${_tag} ${SIMC_REFERENCE_DIR}/runs/${_run} ${_tag} ${_dsp} ${_run})
    list(APPEND _delta_l1 l1_replay_event_capture_delta_${_tag} l1_simcrun_replay_delta_${_tag}
         l1_montecarlo_replay_delta_${_tag})
endforeach()
set_tests_properties(${_delta_l1} PROPERTIES SKIP_RETURN_CODE 77 LABELS l1 TIMEOUT 900)
foreach(_tag d1 d2)
    set_tests_properties(l1_replay_event_capture_delta_${_tag} l1_simcrun_replay_delta_${_tag}
        l1_montecarlo_replay_delta_${_tag} PROPERTIES LABELS "l1;replay")
endforeach()
