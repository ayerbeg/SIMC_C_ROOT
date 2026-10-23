# tests/l1/sidis.cmake — task P12.1: the SIDIS physics (semi_physics.f, semi_dilution.f, cteq5/Ctq5Pdf.f, fdss/fdss.f,
# F1F2IN21_v1.0.f free nucleon) and the IPhysicsModel "sidis". Included (OPTIONAL) by tests/l1/CMakeLists.txt.
# Now: data-file identity, table reading and physical sanity (labels unit). The bitwise L1 tests against the P12.F1
# tables (reference/functions/sidis/, fortran_compat, tolerance 0) are added when the tables exist.
add_executable(test_l1_sidis test_l1_sidis.cpp)
target_link_libraries(test_l1_sidis PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})

set(_sdata ${CMAKE_SOURCE_DIR}/data/sidis)
add_test(NAME l1_sidis_data_files COMMAND test_l1_sidis files ${_sdata} ${SIMC_FORTRAN_DIR})
add_test(NAME sidis_tables COMMAND test_l1_sidis tables ${_sdata})
add_test(NAME sidis_sanity COMMAND test_l1_sidis sanity ${_sdata})
set_tests_properties(l1_sidis_data_files PROPERTIES SKIP_RETURN_CODE 77 LABELS l1)
set_tests_properties(sidis_tables sidis_sanity PROPERTIES SKIP_RETURN_CODE 77 LABELS unit)

# Bitwise L1 against the P12.F1 tables (reference/functions/sidis/, fortran_compat, tolerance 0). The configuration
# inputs (sidis_state/configs.txt) are looked up in sidis_state/inp, simc_gfortran/infiles and the harness p12inp dir.
if(NOT DEFINED SIMC_HARNESS_DIR)
    get_filename_component(SIMC_HARNESS_DIR "${SIMC_REFERENCE_DIR}/../fortran_harness" ABSOLUTE)
endif()
set(_sfn ${SIMC_REFERENCE_DIR}/functions/sidis)
set(_sinp "${_sfn}/sidis_state/inp:${SIMC_FORTRAN_DIR}/infiles:${SIMC_HARNESS_DIR}/rundir_P12F1/p12inp")
set(_ssp ${CMAKE_SOURCE_DIR}/data/spectral)
add_test(NAME l1_ctq5pdf COMMAND test_l1_sidis ctq5pdf ${_sfn}/ctq5pdf ${_sdata} 0)
add_test(NAME l1_ctq5_state COMMAND test_l1_sidis ctq5_state ${_sfn}/ctq5pdf ${_sdata})
add_test(NAME l1_fdss COMMAND test_l1_sidis fdss ${_sfn}/fdss ${_sdata} ${SIMC_FORTRAN_DIR}/fdss 0)
add_test(NAME l1_f1f2in21 COMMAND test_l1_sidis f1f2in21 ${_sfn}/f1f2in21 0)
add_test(NAME l1_slacemcfit_rhad_global COMMAND test_l1_sidis slacemcfit ${_sfn}/slacemcfit_rhad_global 0)
add_test(NAME l1_peepix COMMAND test_l1_sidis peepix ${_sfn}/peepix ${_sfn}/sidis_state "${_sinp}" ${_sdata} ${_ssp} 0)
set_tests_properties(l1_ctq5pdf l1_ctq5_state l1_fdss l1_f1f2in21 l1_slacemcfit_rhad_global l1_peepix
    PROPERTIES SKIP_RETURN_CODE 77 LABELS l1)
# D-051 F1/F4/F6 (FORTRAN-BUG, corrected by default): both modes on the peepix rows (F2 in sidis_sanity)
add_test(NAME fortran_diff_D051_sidis
         COMMAND test_l1_sidis fortran_diff ${_sfn}/peepix ${_sfn}/sidis_state "${_sinp}" ${_sdata} ${_ssp})
set_tests_properties(fortran_diff_D051_sidis PROPERTIES SKIP_RETURN_CODE 77 LABELS "fortran_diff;unit")
