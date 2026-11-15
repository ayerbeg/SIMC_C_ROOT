# simcfw_pion (task P9.1): exclusive pion cross sections (physics_pion.f: peepi, sig_param_2021/exclfit, sigmaid with the
# MAID tables, and the uncalled sig_blok, sig_param04, sig_param_3000, sig_multipole), pizero_decay.f, and the
# IPhysicsModel plug-in "pion". transform_to_cm is the reactions::transformToCm port (simcfw_reactions, P4.1).
# OBJECT library so the SIMC_REGISTER registration is always linked. ROOT-free.
# Consumers link simcfw_pion plus ${SIMC_DRIVER_LIBS} (or simcfw_reactions and its dependencies).
add_library(simcfw_pion OBJECT
    PionCrossSections.cpp
    SigMultipole.cpp
    PizeroDecay.cpp
    PionModel.cpp
)
target_link_libraries(simcfw_pion PUBLIC simcfw simcfw_kinematics simcfw_reactions)
# Default directory of the MAID tables and MULTIPOL.file (data/pion: byte copies of the Fortran files): DataDir("pion")
# (include/simc/core/DataPath.h; P15.4, D-069).
target_compile_options(simcfw_pion PRIVATE -Wall -Wextra -pedantic)
