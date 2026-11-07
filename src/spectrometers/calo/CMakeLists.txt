# simcfw_calo (task P14.3): ISpectrometer "calo" (calorimeter arms 7/8, calo/mc_calo.f). The pi0 two-photon logic of
# simc.f is driver::PizeroCalorimeter (src/driver). OBJECT library so the registration is always linked.
add_library(simcfw_calo OBJECT
    CaloSpectrometer.cpp   # ISpectrometer "calo"
    McCalo.cpp             # kernel: mc_calo.f, mc_calo_recon.f
)
target_link_libraries(simcfw_calo PUBLIC simcfw simcfw_optics)
target_compile_options(simcfw_calo PRIVATE -Wall -Wextra -pedantic)
