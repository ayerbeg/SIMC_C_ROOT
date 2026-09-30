# simcfw_material (task P3.1): IRandom implementations (src/core/random), the SIMC material model
# (enerloss_new, musc, musc_ext) and the Hall C cryo target (trip_thru_target, extreme_trip_thru_target,
# target_musc, target_init). OBJECT library so the SIMC_REGISTER registrations are always linked
# ("mt64" IRandom, "simc_material" IMaterialModel, "hallc_cryo" ITarget).
# The random sources live in src/core/random but are built here because this task may not edit the
# top-level CMakeLists.txt; the director may move them to their own target later.
add_library(simcfw_material OBJECT
    ${CMAKE_SOURCE_DIR}/src/core/random/Gauss1.cpp
    ${CMAKE_SOURCE_DIR}/src/core/random/Mt64Random.cpp
    ${CMAKE_SOURCE_DIR}/src/core/random/ScriptedRandom.cpp
    FortranFormat.cpp
    SimcMaterial.cpp
    HallCCryoTarget.cpp
)
target_link_libraries(simcfw_material PUBLIC simcfw)
target_compile_options(simcfw_material PRIVATE -Wall -Wextra -pedantic)
