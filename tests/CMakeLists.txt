# tests/CMakeLists.txt
# Simple standalone test system without Google Test dependency

message(STATUS "Building standalone tests (no Google Test required)")

# Simple test executable
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

message(STATUS "Tests will run without external test framework")
message(STATUS "Run with: cd build && ./simc_tests")
