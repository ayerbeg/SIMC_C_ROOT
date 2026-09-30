# tests/l1 — L1 (function-level) tests against the Fortran harness tables in
# ${SIMC_REFERENCE_DIR}/functions/<routine>/ (reference/functions/README.md). Tests return 77
# (CTest SKIP) when a table is missing.

add_library(simc_l1runner STATIC L1Runner.cpp)
target_include_directories(simc_l1runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(simc_l1runner PRIVATE -Wall -Wextra -pedantic)

add_executable(test_l1_runner test_l1_runner.cpp)
target_link_libraries(test_l1_runner PRIVATE simc_l1runner)
add_test(NAME l1_runner_selftest COMMAND test_l1_runner ${CMAKE_CURRENT_BINARY_DIR})

# physics_proton.f: sigMott, fofa_best_fit, sigep, deForest (D-010: |rel| <= 1e-12)
add_executable(test_l1_physics_proton test_l1_physics_proton.cpp)
target_link_libraries(test_l1_physics_proton PRIVATE simc_l1runner simcfw_physics)
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
