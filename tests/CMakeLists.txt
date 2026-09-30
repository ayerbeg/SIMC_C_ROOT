# tests/CMakeLists.txt — LEGACY tests (built only with SIMC_BUILD_LEGACY=ON).
# They exercise the pre-translation classes (simc_core/simc_io) and are unverified against the
# Fortran (D-004). New tests go in tests/unit and tests/l1.
#
# P1.7 housekeeping: every test executable built here is registered in CTest. Removed from the
# build: test_integrated_generation (hung: nonsense MeV/GeV momenta, never terminated),
# test_event_generator_diagnostic (a print-only diagnostic, no checks), test_sos_transport /
# test_hrsr_transport / test_hrsl_transport (FATAL: SOS/HRS matrix files not in data/; SOS/HRS are
# out of scope until P14, sources kept for then). The other unbuilt *.cpp here are older one-off
# diagnostics that were already outside the build.

add_executable(simc_tests
    test_main.cpp
    test_event.cpp
    test_config.cpp
    test_random.cpp
    test_cross_section.cpp
    test_energy_loss.cpp
    test_multiple_scattering.cpp
    test_spectrometer_optics.cpp
)
target_link_libraries(simc_tests PRIVATE simc_core simc_io)
add_test(NAME legacy_AllTests COMMAND simc_tests)

add_executable(test_event_generator test_event_generator.cpp)
target_link_libraries(test_event_generator PRIVATE simc_core simc_io)
add_test(NAME legacy_EventGeneratorTests COMMAND test_event_generator)

set_tests_properties(legacy_AllTests legacy_EventGeneratorTests PROPERTIES TIMEOUT 120)
