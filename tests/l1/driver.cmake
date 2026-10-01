# tests/l1/driver.cmake — task P7.1: simc.f:montecarlo (simc::driver::MonteCarlo) replayed against the real-call
# records of the P4.F1 captures (reference/functions/event_capture/<cfg>, reference/runs/p4f1_<cfg>_cap).
# Included (OPTIONAL) by tests/l1/CMakeLists.txt. Everything bitwise; needs `xz` for the .txt.xz capture files.
add_executable(test_l1_driver test_l1_driver.cpp)
target_link_libraries(test_l1_driver PRIVATE simcfw_driver simcfw_setup simcfw_hms simcfw_shms simcfw_optics
    simcfw_material)
target_compile_options(test_l1_driver PRIVATE -O2)
foreach(_c norad rad phsp phsp_tries)
    add_test(NAME l1_montecarlo_replay_${_c} COMMAND test_l1_driver replay
        ${SIMC_REFERENCE_DIR}/functions/event_capture/${_c} ${SIMC_REFERENCE_DIR}/runs/p4f1_${_c}_cap ${_c}
        ${SIMC_FORTRAN_DIR})
    set_tests_properties(l1_montecarlo_replay_${_c} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)
endforeach()

# task P7.2, VALIDATION A: the complete driver (SimcRun) replays the four capture runs with the captured RANLUX
# stream; tree, .bin values, per-try bookkeeping, normalisation, .hist/.geni/.gen text must equal the Fortran.
add_executable(test_l1_simcrun test_l1_simcrun.cpp)
target_link_libraries(test_l1_simcrun PRIVATE ${SIMC_DRIVER_LIBS})
target_compile_options(test_l1_simcrun PRIVATE -O2)
target_compile_definitions(test_l1_simcrun PRIVATE SIMC_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
foreach(_c norad rad phsp phsp_tries)
    set(_out ${CMAKE_CURRENT_BINARY_DIR}/simcrun_replay_${_c})
    file(MAKE_DIRECTORY ${_out})
    add_test(NAME l1_simcrun_replay_${_c} COMMAND test_l1_simcrun replay
        ${SIMC_REFERENCE_DIR}/functions/event_capture/${_c} ${SIMC_REFERENCE_DIR}/runs/p4f1_${_c}_cap ${_c} ${_out})
    set_tests_properties(l1_simcrun_replay_${_c} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900)
endforeach()

# task P7.2: the simc++ CLI, .inp vs the equivalent JSON (data/config) give identical outputs
foreach(_p "norad;eep_h_elastic_norad_A;h_elastic_norad" "rad;eep_h_elastic_rad_A;h_elastic_rad"
           "phsp;acc_phsp_hms_shms_A;acceptance_phsp_hms_shms")
    list(GET _p 0 _tag)
    list(GET _p 1 _ref)
    list(GET _p 2 _cfg)
    add_test(NAME cli_json_vs_inp_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cli_json_vs_inp.sh $<TARGET_FILE:simcpp>
                     ${SIMC_REFERENCE_DIR}/runs/${_ref}/${_ref}.inp ${CMAKE_SOURCE_DIR}/data/config/${_cfg}.json
                     ${CMAKE_CURRENT_BINARY_DIR}/cli_json_vs_inp_${_tag} 2000)
    set_tests_properties(cli_json_vs_inp_${_tag} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)
endforeach()
