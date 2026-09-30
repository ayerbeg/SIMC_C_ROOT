# tests/l1/material.cmake — task P3.1: gauss1, enerloss_new, musc, musc_ext, target_musc,
# trip_thru_target, target_init, extreme_trip_thru_target (+ real-RNG distributions).
# Included (OPTIONAL) by tests/l1/CMakeLists.txt. Tolerance |rel| <= 1e-12 (D-010); RNG consumption exact.
add_executable(test_l1_material test_l1_material.cpp)
target_link_libraries(test_l1_material PRIVATE simc_l1runner simcfw_material ROOT::MathCore)

set(_fnm ${SIMC_REFERENCE_DIR}/functions)
set(_l1mat_tests
    "gauss1_fixedrng|gauss1"
    "enerloss_new|enerloss_new"
    "enerloss_new_fixedrng|enerloss_fixedrng"
    "musc_fixedrng|musc"
    "musc_ext_fixedrng|musc_ext"
    "target_musc_fixedrng|target_musc"
    "trip_thru_target_fixedrng|trip"
    "target_init|target_init"
    "extreme_trip_thru_target|extreme")
foreach(entry ${_l1mat_tests})
    string(REPLACE "|" ";" _p "${entry}")
    list(GET _p 0 _table)
    list(GET _p 1 _mode)
    add_test(NAME l1_${_table} COMMAND test_l1_material ${_mode} ${_fnm}/${_table} 1e-12)
    set_tests_properties(l1_${_table} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()

# Real-RNG distributions (Mt64Random, 1e6 samples): quantiles, binned KS and chi2 vs the Fortran histograms.
foreach(_id 1 2 3 4 5 6 7 8 11 12 13 21 22 23)
    add_test(NAME l1_distributions_${_id} COMMAND test_l1_material distributions ${_fnm}/distributions ${_id})
    set_tests_properties(l1_distributions_${_id} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()

add_test(NAME material_selftest COMMAND test_l1_material selftest)
