# tests/l1/hrs.cmake — task P14.2: ISpectrometer "hrsr"/"hrsl" (mc_hrsr.f / mc_hrsl.f, *_hut.f, *_recon.f).
# Included (OPTIONAL) by tests/l1/CMakeLists.txt.
#   hrs_selftest : REAL(16) quirk values, side constants, RNG use, both modes, counters, plug-ins (no table)
#   l1_<table>   : the P14.F1 tables reference/functions/hrs{r,l}/ (bitwise, tolerance 0; fortran_compat)
#   l1_simcrun_ranlux_hrs_LeRp_rad_cap : RANLUX + compat full run == the Fortran run (text outputs and tree)
add_executable(test_l1_hrs test_l1_hrs.cpp)
target_link_libraries(test_l1_hrs PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})
target_compile_options(test_l1_hrs PRIVATE -O2)
add_test(NAME hrs_selftest COMMAND test_l1_hrs selftest)
set_tests_properties(hrs_selftest PROPERTIES LABELS "l1;unit")
foreach(_s hrsr hrsl)
    foreach(_t "mc;mc_${_s}" "mc;mc_${_s}_smear" "mc;mc_${_s}_decay" "transp;transp_${_s}" "chain;transp_${_s}_chain"
            "recon;recon_${_s}" "hut;hut_${_s}")
        list(GET _t 0 _mode)
        list(GET _t 1 _tab)
        add_test(NAME l1_${_tab} COMMAND test_l1_hrs ${_mode} ${SIMC_REFERENCE_DIR}/functions/${_s}/${_tab} 0)
        set_tests_properties(l1_${_tab} PROPERTIES LABELS l1 SKIP_RETURN_CODE 77)
    endforeach()
endforeach()
# full simc++ run of the P14.F1 capture run (HRS-L e-, HRS-R p, H(e,e'p) rad, 1000 events) with RANLUX + --fortran-compat:
# .hist/.geni/.gen text and the h10 tree identical to the Fortran outputs; corrected-mode .geni slit lines = summary
add_test(NAME l1_simcrun_ranlux_hrs_LeRp_rad_cap
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ranlux_fullrun.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                 ${SIMC_REFERENCE_DIR}/runs/p14f1_hrs_LeRp_rad_cap p14f1_hrs_LeRp_rad_cap
                 ${CMAKE_CURRENT_BINARY_DIR}/ranlux_hrs_LeRp_rad_cap)
set_tests_properties(l1_simcrun_ranlux_hrs_LeRp_rad_cap PROPERTIES LABELS "l1;replay" SKIP_RETURN_CODE 77 TIMEOUT 600)
