# simcfw_radiation (task P5.1): IRadiator "simc_rad" — radc.f, brem.f, init.f radc_init/radc_init_ev/basicrad_init_ev.
# OBJECT library so the SIMC_REGISTER registration is always linked. ROOT-free.
add_library(simcfw_radiation OBJECT
    Brem.cpp
    Radc.cpp
    RadcInit.cpp
    SimcRadiator.cpp
)
target_link_libraries(simcfw_radiation PUBLIC simcfw simcfw_material)
target_compile_options(simcfw_radiation PRIVATE -Wall -Wextra -pedantic)
