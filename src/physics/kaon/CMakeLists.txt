# simcfw_kaon (task P10.1): kaon electroproduction cross sections (physics_kaon.f: peeK, sig_factorized, eekeek,
# eekeeks, phspwght), CERNLIB fint (cern/fint.f), the /sigkaon/ tables (dbase.f:624-679) and the IPhysicsModel plug-in
# "kaon". transform_to_cm is the reactions::transformToCm port (simcfw_reactions, P4.1).
# OBJECT library so the SIMC_REGISTER registration is always linked. ROOT-free. The Saghai amplitudes are evaluated in
# binary128 (__float128, glibc hypotf128 from libm), as the Fortran object (D-043 F4).
add_library(simcfw_kaon OBJECT
    KaonCrossSections.cpp
    KaonModel.cpp
)
target_link_libraries(simcfw_kaon PUBLIC simcfw simcfw_kinematics simcfw_reactions m)
# Default directory of the kaon model tables (data/kaon: byte copies of the Fortran files): DataDir("kaon") (P15.4).
target_compile_options(simcfw_kaon PRIVATE -Wall -Wextra -pedantic)
