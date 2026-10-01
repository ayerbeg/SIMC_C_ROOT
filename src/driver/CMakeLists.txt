# simcfw_driver (tasks P7.1/P7.2): the per-event driver.
#   MonteCarlo.cpp : simc.f:montecarlo — both arms through target multiple scattering, ISpectrometer "hms"/"shms"
#                    (from the registry) and the track back to the target (P7.1).
# Consumers must also link simcfw_hms, simcfw_shms (OBJECT libraries: their registrations), simcfw_optics,
# simcfw_material and simcfw_setup.
add_library(simcfw_driver STATIC
    MonteCarlo.cpp
)
target_link_libraries(simcfw_driver PUBLIC simcfw simcfw_kinematics simcfw_material)
target_compile_options(simcfw_driver PRIVATE -Wall -Wextra -pedantic)
