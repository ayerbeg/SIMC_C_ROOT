# simcfw_setup (task P2.1): run setup from the input parameters (dbase.f / init.f equivalents)
# and the deterministic kinematics utilities (physics_angles, spectrometer_angles, loren).
add_library(simcfw_setup STATIC
    SetupBuilder.cpp
    KinematicsUtils.cpp
)
target_link_libraries(simcfw_setup PUBLIC simcfw)
target_compile_options(simcfw_setup PRIVATE -Wall -Wextra -pedantic)
