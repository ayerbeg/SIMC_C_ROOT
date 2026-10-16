# simcfw_delta (task P11.1): the doing_delta "cross section" physics_delta.f:peedelta (D-048 F3/F4) and the IPhysicsModel
# plug-in "delta". It reuses sig_blok (simcfw_pion) and transform_to_cm (reactions::transformToCm, simcfw_reactions).
# OBJECT library so the SIMC_REGISTER registration is always linked. ROOT-free. The sig_blok argument quirk of the Fortran
# (binary128 temporaries passed to REAL*8 dummies, D-048 F3) is reproduced with __float128 in fortran_compat mode.
add_library(simcfw_delta OBJECT
    DeltaModel.cpp
)
target_link_libraries(simcfw_delta PUBLIC simcfw simcfw_kinematics simcfw_reactions simcfw_pion)
target_compile_options(simcfw_delta PRIVATE -Wall -Wextra -pedantic)
