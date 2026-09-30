# simcfw_optics (task P6.1): the shared COSY optics engine — transp / transp_init (shared/transp.f),
# project (shared/project.f), rotate_haxis / rotate_vaxis, mc_hms_recon / mc_shms_recon arithmetic — plus the
# aperture-free HMS/SHMS forward chains. ROOT-free.
add_library(simcfw_optics STATIC
    CosyEngine.cpp
    ForwardChains.cpp
)
target_link_libraries(simcfw_optics PUBLIC simcfw)
target_compile_options(simcfw_optics PRIVATE -Wall -Wextra -pedantic)
