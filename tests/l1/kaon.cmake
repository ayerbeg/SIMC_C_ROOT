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
