# simcfw_kinematics (P2.1/P2.2): deterministic kinematics utilities (physics_angles, spectrometer_angles,
# loren). The single loren implementation, also used by src/optics.
add_library(simcfw_kinematics STATIC KinematicsUtils.cpp)
target_link_libraries(simcfw_kinematics PUBLIC simcfw)
target_compile_options(simcfw_kinematics PRIVATE -Wall -Wextra -pedantic)

# simcfw_setup (P2.1/P2.2): run setup from the input parameters (dbase.f / init.f equivalents), with
# material::HallCCryoTarget as the default target_init / extreme_trip_thru_target provider.
add_library(simcfw_setup STATIC
    SetupBuilder.cpp
    HallCTargetLosses.cpp
)
target_link_libraries(simcfw_setup PUBLIC simcfw simcfw_kinematics simcfw_material)
target_compile_options(simcfw_setup PRIVATE -Wall -Wextra -pedantic)
