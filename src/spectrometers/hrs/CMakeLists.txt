# simcfw_hrs (task P14.2): ISpectrometer "hrsr" (arm 3) and "hrsl" (arm 4) — mc_hrsr.f / mc_hrsl.f, *_hut.f, recon via the
# shared COSY engine. One implementation with a side parameter. OBJECT library so the registrations are always linked.
add_library(simcfw_hrs OBJECT
    HrsSpectrometer.cpp   # ISpectrometer "hrsr", "hrsl"
    McHrs.cpp             # kernel: mc_hrs*.f
    McHrsHut.cpp          # kernel: mc_hrs*_hut.f
)
target_link_libraries(simcfw_hrs PUBLIC simcfw simcfw_optics)
target_compile_options(simcfw_hrs PRIVATE -Wall -Wextra -pedantic)
target_compile_definitions(simcfw_hrs PRIVATE
    SIMC_HRS_DEFAULT_MATRIX_DIR="${CMAKE_SOURCE_DIR}/data/matrices/simc_fortran")
