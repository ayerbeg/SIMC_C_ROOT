# simcfw_optics (task P6.1): the shared COSY optics engine — transp / transp_init (shared/transp.f),
# project (shared/project.f), rotate_haxis / rotate_vaxis, mc_hms_recon / mc_shms_recon arithmetic — plus the
# aperture-free HMS/SHMS forward chains. ROOT-free.
add_library(simcfw_optics STATIC
    CosyEngine.cpp
    ForwardChains.cpp
    Lfit.cpp          # CERNLIB lfit (cern/lfit.f), used by mc_hms_hut and mc_shms_hut (moved here in P6.3)
    TargetField.cpp   # polarized-target field tracking, trg_track.f (P15.2b)
)
target_link_libraries(simcfw_optics PUBLIC simcfw simcfw_kinematics)
target_compile_options(simcfw_optics PRIVATE -Wall -Wextra -pedantic)
