# simcfw_rho (task P13.1): exclusive rho0 electroproduction of the Fortran SIMC: rho_physics.f:peerho (the IPhysicsModel
# plug-in "rho"), generate_rho.f:generate_rho and rho_decay.f:rho_decay (kernels called by the rho reaction, P13.2).
# OBJECT library so the SIMC_REGISTER registration is always linked. ROOT-free. The REAL(16) literals 1.d+06 and 1.d3 of
# rho_physics.f (-fdefault-real-8) are reproduced with __float128 divisions.
add_library(simcfw_rho OBJECT
    RhoPhysics.cpp
    RhoModel.cpp
)
target_link_libraries(simcfw_rho PUBLIC simcfw simcfw_kinematics)
target_compile_options(simcfw_rho PRIVATE -Wall -Wextra -pedantic)
