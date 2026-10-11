# simcfw_reactions (task P4.1, P8.2): IBeam "standard", IReaction "hydrogen_elastic"/"phase_space" and
# "deuterium_eep"/"nucleus_eep" (event.f generate, complete_ev, complete_recon_ev, complete_main), limits_update,
# transform_to_cm (jacobians.f). The spectral-function plug-ins ("theory", "benhar") are in simcfw_physics.
# OBJECT library so the SIMC_REGISTER registrations are always linked. ROOT-free.
add_library(simcfw_reactions OBJECT
    StandardBeam.cpp
    EepReactionBase.cpp
    HydrogenElastic.cpp
    NuclearEep.cpp
    EventChain.cpp
)
target_link_libraries(simcfw_reactions PUBLIC simcfw simcfw_kinematics simcfw_material)
target_compile_options(simcfw_reactions PRIVATE -Wall -Wextra -pedantic)

# simcfw_reactions_pion (task P9.2): IReaction "pion_production" (event.f doing_pion branches). Separate OBJECT library
# because it uses the pion kernels (pizero_decay, simcfw_pion) and the Benhar table of generate_em (simcfw_physics).
add_library(simcfw_reactions_pion OBJECT pion/PionProduction.cpp)
target_link_libraries(simcfw_reactions_pion PUBLIC simcfw simcfw_kinematics simcfw_material simcfw_reactions
    simcfw_physics simcfw_pion)
target_compile_definitions(simcfw_reactions_pion PRIVATE SIMC_PION_DATA_DIR="${CMAKE_SOURCE_DIR}/data/pion")
target_compile_options(simcfw_reactions_pion PRIVATE -Wall -Wextra -pedantic)

# simcfw_reactions_kaon (task P10.2): IReaction "kaon_production" (event.f doing_kaon branches = the pion code with the
# kaon flags, plus the complete_main kaon branch).
add_library(simcfw_reactions_kaon OBJECT kaon/KaonProduction.cpp)
target_link_libraries(simcfw_reactions_kaon PUBLIC simcfw simcfw_reactions simcfw_reactions_pion)
target_compile_options(simcfw_reactions_kaon PRIVATE -Wall -Wextra -pedantic)
