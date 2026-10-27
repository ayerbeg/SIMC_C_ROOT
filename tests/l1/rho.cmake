# tests/l1/rho.cmake — task P13.1: the exclusive rho0 kernels (rho_physics.f:peerho, generate_rho.f, rho_decay.f) and the
# IPhysicsModel "rho". Included (OPTIONAL) by tests/l1/CMakeLists.txt.
# Now: kinematic identities, RNG order/counts, failure paths, the D-031 switches (unit), and the rho set-up/plug-in
# checks. The bitwise L1 tests against the P13.F1 tables (reference/functions/rho/, fortran_compat, tolerance 0) are added
# when the tables exist.
add_executable(test_l1_rho test_l1_rho.cpp)
target_link_libraries(test_l1_rho PRIVATE simcfw_l1runner ${SIMC_DRIVER_LIBS})
target_compile_options(test_l1_rho PRIVATE -Wall -Wextra -pedantic)

add_test(NAME rho_sanity COMMAND test_l1_rho sanity)
add_test(NAME rho_setup COMMAND test_l1_rho setup ${SIMC_REFERENCE_DIR}/runs/pion_h_k1_rad_dec1_A/pion_h_k1_rad_dec1_A.inp)
set_tests_properties(rho_sanity rho_setup PROPERTIES SKIP_RETURN_CODE 77 LABELS unit)
# P13.1 R1/R2 (FORTRAN-BUG, corrected by default): both modes
add_test(NAME fortran_diff_P131_R1_R2_rho COMMAND test_l1_rho fortran_diff)
set_tests_properties(fortran_diff_P131_R1_R2_rho PROPERTIES LABELS "fortran_diff;unit")
