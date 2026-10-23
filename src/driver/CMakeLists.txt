# simcfw_driver (tasks P7.1/P7.2): the SIMC main program.
#   MonteCarlo.cpp    : simc.f:montecarlo — both arms through target multiple scattering, ISpectrometer "hms"/"shms"
#                       (from the registry) and the track back to the target (P7.1).
#   Central.cpp       : simc.f:calculate_central (P7.2)
#   RunReport.cpp     : simc.f normalisation (l.360-432), report (.hist) and .geni text, summary sections (P7.2)
#   FortranWriter.cpp : gfortran formatted-output emulation for the .hist/.geni text
#   Pipeline.cpp      : the JSON plug-in choice per stage (docs/DRIVER.md)
#   SimcRun.cpp       : the main loop (simc.f l.1-623)
# Consumers must also link simcfw_hms, simcfw_shms (OBJECT libraries: their registrations), simcfw_optics,
# simcfw_material and simcfw_setup — or simply ${SIMC_DRIVER_LIBS} (everything a full run needs).
add_library(simcfw_driver STATIC
    MonteCarlo.cpp
    Central.cpp
    FortranWriter.cpp
    Pipeline.cpp
    RunReport.cpp
    SimcRun.cpp
)
target_link_libraries(simcfw_driver PUBLIC simcfw simcfw_kinematics simcfw_material simcfw_setup)
# Default directory of the spectral-function files (data/spectral: copies of the Fortran theory/Benhar files, P8.1)
target_compile_definitions(simcfw_driver PRIVATE SIMC_SPECTRAL_DATA_DIR="${CMAKE_SOURCE_DIR}/data/spectral")
target_compile_options(simcfw_driver PRIVATE -Wall -Wextra -pedantic)

# The "fortran_hist" sink (SIMC_REGISTER): OBJECT library so the registration is linked.
add_library(simcfw_driver_sinks OBJECT ReportSink.cpp)
target_link_libraries(simcfw_driver_sinks PUBLIC simcfw_driver)
target_compile_options(simcfw_driver_sinks PRIVATE -Wall -Wextra -pedantic)

# Everything a complete run links (OBJECT libraries must be linked directly by the executable).
set(SIMC_DRIVER_LIBS simcfw_driver simcfw_driver_sinks simcfw_reactions simcfw_reactions_pion simcfw_reactions_kaon simcfw_reactions_delta simcfw_reactions_sidis simcfw_radiation
    simcfw_physics simcfw_pion simcfw_kaon simcfw_delta simcfw_sidis simcfw_root
    simcfw_hms simcfw_shms simcfw_optics simcfw_setup simcfw_material simcfw_kinematics simcfw_random simcfw
    CACHE INTERNAL "libraries of a complete SIMC++ run")
# The command line simc++ is defined in the top-level CMakeLists.txt (it also links the example plug-ins).
