# tests/l1 — L1 (function-level) tests against the Fortran harness tables in
# ${SIMC_REFERENCE_DIR}/functions/<routine>/ (reference/functions/README.md). Tests return 77
# (CTest SKIP) when a table is missing.

add_library(simcfw_l1runner STATIC L1Runner.cpp)
target_include_directories(simcfw_l1runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(simcfw_l1runner PRIVATE -Wall -Wextra -pedantic)

add_executable(test_l1_runner test_l1_runner.cpp)
target_link_libraries(test_l1_runner PRIVATE simcfw_l1runner)
add_test(NAME l1_runner_selftest COMMAND test_l1_runner ${CMAKE_CURRENT_BINARY_DIR})

# physics_proton.f: sigMott, fofa_best_fit, sigep, deForest (D-010: |rel| <= 1e-12)
add_executable(test_l1_physics_proton test_l1_physics_proton.cpp)
target_link_libraries(test_l1_physics_proton PRIVATE simcfw_l1runner simcfw_physics simcfw_material simcfw_random)
set(_fn ${SIMC_REFERENCE_DIR}/functions)
foreach(r sigmott fofa_best_fit sigep deforest)
    add_test(NAME l1_${r} COMMAND test_l1_physics_proton ${r} ${_fn}/${r} 1e-12)
    set_tests_properties(l1_${r} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()
add_test(NAME physics_proton_plugins COMMAND test_l1_physics_proton plugins)

# Per-module L1 tests (separate files so parallel tasks don't collide; director, P2)
include(${CMAKE_CURRENT_SOURCE_DIR}/material.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/setup.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/optics.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/hms.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/shms.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/radiation.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/reactions.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/driver.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/eep.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/colldecay.cmake OPTIONAL)

# CTest labels (P7.4): every test of this directory is "l1"; the bitwise replays of the Fortran captures are also
# "replay", the simc++ command-line comparisons "cli" (and the long distribution tests keep "slow").
get_property(_l1_tests DIRECTORY PROPERTY TESTS)
set_property(TEST ${_l1_tests} APPEND PROPERTY LABELS l1)
foreach(_t ${_l1_tests})
    if(_t MATCHES "_replay_")
        set_property(TEST ${_t} APPEND PROPERTY LABELS replay)
    elseif(_t MATCHES "^cli_")
        set_property(TEST ${_t} APPEND PROPERTY LABELS cli)
    endif()
endforeach()
include(${CMAKE_CURRENT_SOURCE_DIR}/pion.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/kaon.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/delta.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/sidis.cmake OPTIONAL)
include(${CMAKE_CURRENT_SOURCE_DIR}/rho.cmake OPTIONAL)
