# simcfw_sos (task P14.1): ISpectrometer "sos" — mc_sos, mc_sos_hut (CERNLIB lfit is in simcfw_optics), recon via the
# shared COSY engine (simcfw_optics). OBJECT library so the SIMC_REGISTER registration is always linked.
# Consumers must also link simcfw_optics and simcfw_material (object files of OBJECT libraries do not propagate).
add_library(simcfw_sos OBJECT
    SosSpectrometer.cpp   # ISpectrometer "sos" plug-in
    McSos.cpp             # kernel: mc_sos.f
    McSosHut.cpp          # kernel: mc_sos_hut.f
)
target_link_libraries(simcfw_sos PUBLIC simcfw simcfw_optics)
target_compile_options(simcfw_sos PRIVATE -Wall -Wextra -pedantic)
# Default location of sos/forward_cosy.dat and sos/recon_cosy.dat (byte-identical copies of the Fortran files):
# DataDir("matrices/simc_fortran") (include/simc/core/DataPath.h; P15.4, D-069).
