# tests/CMakeLists.txt
# Simple standalone test system without Google Test dependency

message(STATUS "Building standalone tests (no Google Test required)")

# Simple test executable - ALL EXISTING TESTS
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

# Link libraries
target_link_libraries(simc_tests
    PRIVATE
        simc_core
        simc_io
)

# Add tests to CTest (will run all tests in one executable)
add_test(NAME AllTests COMMAND simc_tests)

# =============================================================================
# NEW: Event Generator test executable (separate)
# =============================================================================
add_executable(test_event_generator
    test_event_generator.cpp
)

target_link_libraries(test_event_generator
    PRIVATE
        simc_core
        simc_io
)

add_executable(test_event_generator_diagnostic test_event_generator_diagnostic.cpp)
target_link_libraries(test_event_generator_diagnostic simc_core)

# Add as separate test
add_test(NAME EventGeneratorTests COMMAND test_event_generator)

message(STATUS "Tests will run without external test framework")
message(STATUS "Run with: cd build && ./tests/simc_tests")
message(STATUS "Run event generator tests: cd build/tests && ./test_event_generator")
