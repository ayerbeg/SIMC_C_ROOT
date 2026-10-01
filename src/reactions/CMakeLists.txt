# simcfw_reactions (task P4.1): IBeam "standard", IReaction "hydrogen_elastic"/"phase_space" (event.f generate,
# complete_ev, complete_recon_ev), complete_main, limits_update, transform_to_cm (jacobians.f).
# OBJECT library so the SIMC_REGISTER registrations are always linked. ROOT-free.
add_library(simcfw_reactions OBJECT
    StandardBeam.cpp
    HydrogenElastic.cpp
    EventChain.cpp
)
target_link_libraries(simcfw_reactions PUBLIC simcfw simcfw_kinematics simcfw_material)
target_compile_options(simcfw_reactions PRIVATE -Wall -Wextra -pedantic)
