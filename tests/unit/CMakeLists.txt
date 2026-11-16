# tests/unit — framework unit tests (CTest). Reference data outside the repo is passed as
# arguments; a test returns 77 (SKIP) when that data is missing.

function(simc_unit_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE simcfw)
endfunction()

simc_unit_test(test_parameters)
simc_unit_test(test_inp_reader)
simc_unit_test(test_inp_reference)
simc_unit_test(test_stage_headers)

set(_csv ${SIMC_REFERENCE_DIR}/functions/_inputs/inp_variables.csv)
set(_runs ${SIMC_REFERENCE_DIR}/runs)

add_test(NAME params_completeness COMMAND test_parameters completeness ${_csv})
add_test(NAME params_json_roundtrip COMMAND test_parameters json)
add_test(NAME inp_reader_ctp_semantics COMMAND test_inp_reader)
add_test(NAME inp_convert_all_infiles COMMAND test_inp_reference all_infiles ${SIMC_FORTRAN_DIR}/infiles)
add_test(NAME inp_vs_hist_eep_h_test_pctnames COMMAND test_inp_reference hist
         ${_runs}/eep_h_test_pctnames/eep_h_test_pctnames.inp ${_runs}/eep_h_test_pctnames/eep_h_test_pctnames.hist)
add_test(NAME inp_vs_hist_hms_shms_validation COMMAND test_inp_reference hist
         ${_runs}/hms_shms_validation/test_hms_shms_validation.inp ${_runs}/hms_shms_validation/test_hms_shms_validation.hist)
add_test(NAME inp_not_registered_eep_h_test COMMAND test_inp_reference notreg
         ${_runs}/eep_h_test/test_eep_h.inp ${_runs}/eep_h_test/test_eep_h.out 46)
add_test(NAME inp_not_registered_hms_shms_validation COMMAND test_inp_reference notreg
         ${_runs}/hms_shms_validation/test_hms_shms_validation.inp ${_runs}/hms_shms_validation/test_hms_shms_validation.out 20)
add_test(NAME inp_not_registered_eep_h_test_pctnames COMMAND test_inp_reference notreg
         ${_runs}/eep_h_test_pctnames/eep_h_test_pctnames.inp ${_runs}/eep_h_test_pctnames/eep_h_test_pctnames.out 0)
add_test(NAME stage_headers_and_registry COMMAND test_stage_headers)
# (skipped, exit 77, when the reference input is absent: tests/skip_if_missing.sh)
add_test(NAME inp2json_cli COMMAND ${CMAKE_SOURCE_DIR}/tests/skip_if_missing.sh ${_runs}/eep_h_test_pctnames/eep_h_test_pctnames.inp --
         $<TARGET_FILE:inp2json> -q -o ${CMAKE_BINARY_DIR}/inp2json_eep_h_test_pctnames.json
         --table=${CMAKE_BINARY_DIR}/inp2json_eep_h_test_pctnames.table.txt
         ${_runs}/eep_h_test_pctnames/eep_h_test_pctnames.inp)
set_tests_properties(inp2json_cli PROPERTIES SKIP_RETURN_CODE 77)

set_tests_properties(params_completeness inp_convert_all_infiles
    inp_vs_hist_eep_h_test_pctnames inp_vs_hist_hms_shms_validation
    inp_not_registered_eep_h_test inp_not_registered_hms_shms_validation
    inp_not_registered_eep_h_test_pctnames
    PROPERTIES SKIP_RETURN_CODE 77)

# ---- P1.4 output sinks
add_executable(test_output_sinks test_output_sinks.cpp)
target_link_libraries(test_output_sinks PRIVATE simcfw_root)
add_executable(test_hist2json test_hist2json.cpp)
target_link_libraries(test_hist2json PRIVATE simcfw)

add_test(NAME output_sinks COMMAND test_output_sinks
         ${_runs}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.root ${CMAKE_SOURCE_DIR}/docs/manual/markdown/reference/summary.schema.json)
set(_hist_runs eep_h_elastic_norad_A eep_h_elastic_norad_B eep_h_elastic_rad_A eep_h_elastic_rad_B
               acc_phsp_hms_shms_A acc_phsp_hms_shms_B)
foreach(r ${_hist_runs})
    add_test(NAME hist2json_${r} COMMAND test_hist2json ${_runs}/${r}/${r}.hist ${CMAKE_SOURCE_DIR}/docs/manual/markdown/reference/summary.schema.json)
endforeach()
add_test(NAME hist2json_eep_h_test_pctnames COMMAND test_hist2json
         ${_runs}/eep_h_test_pctnames/eep_h_test_pctnames.hist ${CMAKE_SOURCE_DIR}/docs/manual/markdown/reference/summary.schema.json)
add_test(NAME hist2json_hms_shms_validation COMMAND test_hist2json
         ${_runs}/hms_shms_validation/test_hms_shms_validation.hist ${CMAKE_SOURCE_DIR}/docs/manual/markdown/reference/summary.schema.json)
set_tests_properties(output_sinks hist2json_eep_h_test_pctnames hist2json_hms_shms_validation PROPERTIES SKIP_RETURN_CODE 77)
foreach(r ${_hist_runs})
    set_tests_properties(hist2json_${r} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()
include(${CMAKE_CURRENT_SOURCE_DIR}/setup.cmake OPTIONAL)   # P2.1 setup vs .hist
# P10.4 (D-044): squares are x*x as in gfortran, never std::pow(x, 2) (gcc 16 calls libm pow for it)
add_test(NAME lint_no_pow_square COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/lint_no_pow_square.sh ${CMAKE_SOURCE_DIR})

# CTest labels (P7.4): every test of this directory is "unit"; the inp2json command line also "cli".
get_property(_unit_tests DIRECTORY PROPERTY TESTS)
set_property(TEST ${_unit_tests} APPEND PROPERTY LABELS unit)
set_property(TEST inp2json_cli APPEND PROPERTY LABELS cli)
