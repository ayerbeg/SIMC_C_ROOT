# tests/l1/colldecay.cmake — task P6.4: HMS/SHMS collimator (pi/mu) and decay in flight against the P6.F2 tables
# (reference/functions/{mc_hms_coll,mc_shms_coll,mc_hms_decay,mc_shms_decay,transp_decay_hms,transp_decay_shms,
# project_decay,pion_coll_absorb}), all BITWISE in fortran_compat mode, and the bitwise replay of every mc_hms/mc_shms
# call recorded in the two P6.F2 capture runs (runs/p6f2_pi_{hms,shms}_cap; the call files are kept by the harness:
# ${SIMC_HARNESS_DIR}/spectrometers_coll_decay/work/p6f2_pi_*_cap/cap_calls.txt). Included by tests/l1/CMakeLists.txt.
add_executable(test_l1_coll_decay test_l1_coll_decay.cpp)
target_link_libraries(test_l1_coll_decay PRIVATE simcfw_l1runner simcfw_hms simcfw_shms simcfw_optics simcfw_material
    simcfw_random)
target_compile_options(test_l1_coll_decay PRIVATE -O2 -Wall -Wextra -pedantic)

if(NOT DEFINED SIMC_HARNESS_DIR)
    get_filename_component(SIMC_HARNESS_DIR "${SIMC_REFERENCE_DIR}/../fortran_harness" ABSOLUTE)
endif()
set(_fcd ${SIMC_REFERENCE_DIR}/functions)
set(_cap ${SIMC_HARNESS_DIR}/spectrometers_coll_decay/work)
add_test(NAME l1_mc_hms_coll COMMAND test_l1_coll_decay mc ${_fcd}/mc_hms_coll hms)
add_test(NAME l1_mc_hms_decay COMMAND test_l1_coll_decay mc ${_fcd}/mc_hms_decay hms)
add_test(NAME l1_mc_shms_coll COMMAND test_l1_coll_decay mc ${_fcd}/mc_shms_coll shms)
add_test(NAME l1_mc_shms_decay COMMAND test_l1_coll_decay mc ${_fcd}/mc_shms_decay shms)
add_test(NAME l1_transp_decay_hms COMMAND test_l1_coll_decay tp ${_fcd}/transp_decay_hms)
add_test(NAME l1_transp_decay_shms COMMAND test_l1_coll_decay tp ${_fcd}/transp_decay_shms)
add_test(NAME l1_project_decay COMMAND test_l1_coll_decay tp ${_fcd}/project_decay)
add_test(NAME l1_pion_coll_absorb COMMAND test_l1_coll_decay pca ${_fcd}/pion_coll_absorb)
add_test(NAME l1_replay_p6f2_pi_hms_cap COMMAND test_l1_coll_decay replay ${_cap}/p6f2_pi_hms_cap/cap_calls.txt)
add_test(NAME l1_replay_p6f2_pi_shms_cap COMMAND test_l1_coll_decay replay ${_cap}/p6f2_pi_shms_cap/cap_calls.txt)
set_tests_properties(l1_mc_hms_coll l1_mc_hms_decay l1_mc_shms_coll l1_mc_shms_decay l1_transp_decay_hms
    l1_transp_decay_shms l1_project_decay l1_pion_coll_absorb l1_replay_p6f2_pi_hms_cap l1_replay_p6f2_pi_shms_cap
    PROPERTIES SKIP_RETURN_CODE 77)
