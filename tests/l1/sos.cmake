# tests/l1/sos.cmake — task P14.1: ISpectrometer "sos" (mc_sos, mc_sos_hut, mc_sos_recon on the COSY engine).
# Included (OPTIONAL) by tests/l1/CMakeLists.txt.
#   sos_selftest                : maps, RNG consumption, counters, recon sanity, zero guard, sieve, decay, plug-in (no table)
#   l1_{mc_sos,mc_sos_smear,mc_sos_decay,transp_sos,transp_sos_chain,recon_sos,hut_sos} : the P14.F1 ray tables
#                                 (reference/functions/sos/), bitwise
#   l1_simcrun_ranlux_sos_eep_h_test_pctnames : the full simc++ run of reference/runs/eep_h_test_pctnames (HMS e-,
#                                 SOS p, 3874 tries) with the RANLUX port and --fortran-compat: .hist (except the time
#                                 lines), .geni, .gen and the h10 tree bitwise equal to the Fortran outputs.
add_executable(test_l1_sos test_l1_sos.cpp)
target_link_libraries(test_l1_sos PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})
target_compile_options(test_l1_sos PRIVATE -O2)

set(_fns ${SIMC_REFERENCE_DIR}/functions/sos)
add_test(NAME sos_selftest COMMAND test_l1_sos selftest)
# P14.F1 tables, all compared BITWISE (tolerance 0): every column, stop reasons row by row, counters, RNG consumption
foreach(_t "mc;mc_sos" "mc;mc_sos_smear" "mc;mc_sos_decay" "transp;transp_sos" "chain;transp_sos_chain" "recon;recon_sos"
        "hut;hut_sos")
    list(GET _t 0 _mode)
    list(GET _t 1 _tab)
    add_test(NAME l1_${_tab} COMMAND test_l1_sos ${_mode} ${_fns}/${_tab} 0)
    list(APPEND _sos_l1 l1_${_tab})
endforeach()
set(_sosrun ${SIMC_REFERENCE_DIR}/runs/eep_h_test_pctnames)
add_test(NAME l1_simcrun_ranlux_sos_eep_h_test_pctnames
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ranlux_fullrun.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                 ${_sosrun} eep_h_test_pctnames ${CMAKE_CURRENT_BINARY_DIR}/ranlux_sos_eep_h_test_pctnames)
set_tests_properties(sos_selftest PROPERTIES LABELS "l1;unit")
set_tests_properties(${_sos_l1} PROPERTIES LABELS l1 SKIP_RETURN_CODE 77)
set_tests_properties(l1_simcrun_ranlux_sos_eep_h_test_pctnames PROPERTIES LABELS "l1;replay" SKIP_RETURN_CODE 77
                     TIMEOUT 600)
