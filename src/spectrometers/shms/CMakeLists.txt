# simcfw_shms (task P6.3): ISpectrometer "shms" — mc_shms, mc_shms_coll, mc_shms_hut (hut.inc), recon via the shared
# COSY engine; CERNLIB lfit is in simcfw_optics. mc_shms_coll calls pion_coll_absorb, which lives in simcfw_hms
# (hms/pion_coll_absorb.f is shared by both spectrometers in SIMC).
# OBJECT library so the SIMC_REGISTER registration is always linked. Consumers must also link simcfw_hms,
# simcfw_optics and simcfw_material (object files of OBJECT libraries do not propagate).
add_library(simcfw_shms OBJECT
    ShmsSpectrometer.cpp  # ISpectrometer "shms" plug-in
    McShms.cpp            # kernel: mc_shms.f
    McShmsColl.cpp
    McShmsHut.cpp
)
target_link_libraries(simcfw_shms PUBLIC simcfw simcfw_optics)
target_compile_options(simcfw_shms PRIVATE -Wall -Wextra -pedantic)
# Default matrix directory: DataDir("matrices/simc_fortran") (include/simc/core/DataPath.h; P15.4, D-069).
