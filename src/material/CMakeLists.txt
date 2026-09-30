# simcfw_material (task P3.1): the SIMC material model
# (enerloss_new, musc, musc_ext) and the Hall C cryo target (trip_thru_target, extreme_trip_thru_target,
# target_musc, target_init). OBJECT library so the SIMC_REGISTER registrations are always linked
# ("simc_material" IMaterialModel, "hallc_cryo" ITarget).
# The random sources (src/core/random) are in their own target simcfw_random (top-level CMakeLists, P2.2).
add_library(simcfw_material OBJECT
    FortranFormat.cpp
    SimcMaterial.cpp
    HallCCryoTarget.cpp
)
target_link_libraries(simcfw_material PUBLIC simcfw simcfw_random)
target_compile_options(simcfw_material PRIVATE -Wall -Wextra -pedantic)
