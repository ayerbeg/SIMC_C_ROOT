# simcfw_sidis (task P12.1): the semi-inclusive (e,e'h)X cross section of the Fortran SIMC (semi_physics.f: peepiX,
# slacemcfit, Rhad_global; semi_dilution.f, dead code), the CTEQ5 PDFs (cteq5/Ctq5Pdf.f), the fDSS fragmentation
# functions (fdss/fdss.f), the free-nucleon part of the Christy 2021 inclusive fit (F1F2IN21_v1.0.f: F1F2IN21 A = 1,
# SF, rescsp, rescsn, RESMODP, RESMODN) and the IPhysicsModel plug-in "sidis".
# OBJECT library so the SIMC_REGISTER registration is always linked. ROOT-free. fdss.f's REAL(16) parts are evaluated in
# binary128 (__float128) as the Fortran object.
add_library(simcfw_sidis OBJECT
    Cteq5Pdf.cpp
    Fdss.cpp
    F1F2In21.cpp
    SidisCrossSection.cpp
    SemiDilution.cpp
    SidisModel.cpp
)
target_link_libraries(simcfw_sidis PUBLIC simcfw m)
# Default directory of the SIDIS model tables (data/sidis: byte copies of the Fortran files)
target_compile_definitions(simcfw_sidis PRIVATE SIMC_SIDIS_DATA_DIR="${CMAKE_SOURCE_DIR}/data/sidis")
target_compile_options(simcfw_sidis PRIVATE -Wall -Wextra -pedantic)
