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
