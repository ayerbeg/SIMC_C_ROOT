# simcfw_hms (task P6.2): ISpectrometer "hms" — mc_hms, mc_hms_coll, mc_hms_hut (CERNLIB lfit is in simcfw_optics), recon via the
# shared COSY engine (simcfw_optics). OBJECT library so the SIMC_REGISTER registration is always linked.
# Consumers must also link simcfw_optics and simcfw_material (object files of OBJECT libraries do not propagate).
add_library(simcfw_hms OBJECT
    HmsSpectrometer.cpp   # ISpectrometer "hms" plug-in
    McHms.cpp             # kernel: mc_hms.f, hit_dipole
    McHmsColl.cpp
    McHmsHut.cpp
)
target_link_libraries(simcfw_hms PUBLIC simcfw simcfw_optics)
target_compile_options(simcfw_hms PRIVATE -Wall -Wextra -pedantic)
# Default location of hms/forward_cosy.dat and hms/recon_cosy.dat (byte-identical copies of the Fortran files):
# DataDir("matrices/simc_fortran") (include/simc/core/DataPath.h; P15.4, D-069).
