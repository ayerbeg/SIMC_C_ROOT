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

add_executable(test_integrated_generation 
    test_integrated_generation.cpp
)

target_link_libraries(test_integrated_generation
    simc_core
    simc_io
    ${ROOT_LIBRARIES}
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

add_executable(test_sos_transport test_sos_transport.cpp)
target_link_libraries(test_sos_transport
    simc_core
    ${ROOT_LIBRARIES}
)

add_executable(test_hrsr_transport test_hrsr_transport.cpp)
target_link_libraries(test_hrsr_transport
    simc_core
    ${ROOT_LIBRARIES}
)

add_executable(test_hrsl_transport test_hrsl_transport.cpp)
target_link_libraries(test_hrsl_transport
    simc_core
    ${ROOT_LIBRARIES}
)




# Add as separate test
add_test(NAME EventGeneratorTests COMMAND test_event_generator)

message(STATUS "Tests will run without external test framework")
message(STATUS "Run with: cd build && ./tests/simc_tests")
message(STATUS "Run event generator tests: cd build/tests && ./test_event_generator")
