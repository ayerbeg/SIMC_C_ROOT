# simcfw_driver (tasks P7.1/P7.2): the SIMC main program.
#   MonteCarlo.cpp    : simc.f:montecarlo — both arms through target multiple scattering, ISpectrometer "hms"/"shms"/"sos"
#                       (from the registry) and the track back to the target (P7.1).
#   PizeroCalorimeter.cpp : simc.f:1489-1556 — pi0 -> 2 gamma through a calorimeter arm (P14.3)
#   Central.cpp       : simc.f:calculate_central (P7.2)
#   RunReport*.cpp    : simc.f normalisation (l.360-432), report (.hist), .geni text, summary sections (P7.2; split P15.3)
#   FortranWriter.cpp : gfortran formatted-output emulation for the .hist/.geni text
#   Pipeline.cpp      : the JSON plug-in choice per stage (docs/manual/markdown/reference/DRIVER.md)
#   SimcRun.cpp       : the main loop (simc.f l.1-623)
#   SimcRunParallel.cpp : the parallel event loop (--threads N; not Fortran; P15.4, docs/manual/markdown/reference/PERFORMANCE.md)
# Consumers must also link simcfw_hms, simcfw_shms, simcfw_sos (OBJECT libraries: their registrations), simcfw_optics,
# simcfw_material and simcfw_setup — or simply ${SIMC_DRIVER_LIBS} (everything a full run needs).
add_library(simcfw_driver STATIC
    MonteCarlo.cpp
    PizeroCalorimeter.cpp
    Central.cpp
    CommandLine.cpp         # the simc++ command line (CommandLineMain, P15.3)
    FortranWriter.cpp
    Pipeline.cpp
    RunReport.cpp           # normalisation, summary sections (P15.3 split:)
    RunReportHist.cpp       #   .hist (simc.f:report)
    RunReportGeni.cpp       #   .geni
    RunReportRecords.cpp    #   reaction lines and FORMAT records of the report
    SimcRun.cpp
    SimcRunParallel.cpp     # the parallel event loop (--threads N, P15.4)
    TargetFieldAngles.cpp   # simc.f:120-155 arm angles of the polarized-target field (P15.2b)
)
find_package(Threads REQUIRED)
target_link_libraries(simcfw_driver PUBLIC simcfw simcfw_kinematics simcfw_material simcfw_setup Threads::Threads)
# Default directories of the spectral-function files (data/spectral: copies of the Fortran theory/Benhar files, P8.1) and
# of trg_field_map.dat (data/field, P15.2b): DataDir("spectral"), DataDir("field") (include/simc/core/DataPath.h, P15.4).
target_compile_options(simcfw_driver PRIVATE -Wall -Wextra -pedantic)

# The "fortran_hist" sink (SIMC_REGISTER): OBJECT library so the registration is linked.
add_library(simcfw_driver_sinks OBJECT ReportSink.cpp)
target_link_libraries(simcfw_driver_sinks PUBLIC simcfw_driver)
target_compile_options(simcfw_driver_sinks PRIVATE -Wall -Wextra -pedantic)

# Everything a complete run links (OBJECT libraries must be linked directly by the executable).
set(SIMC_DRIVER_LIBS simcfw_driver simcfw_driver_sinks simcfw_reactions simcfw_reactions_pion simcfw_reactions_kaon simcfw_reactions_delta simcfw_reactions_sidis simcfw_reactions_rho simcfw_radiation
    simcfw_physics simcfw_pion simcfw_kaon simcfw_delta simcfw_sidis simcfw_rho simcfw_root
    simcfw_hms simcfw_shms simcfw_sos simcfw_hrs simcfw_calo simcfw_optics simcfw_setup simcfw_material simcfw_kinematics simcfw_random simcfw
    CACHE INTERNAL "libraries of a complete SIMC++ run")
# The command line simc++ is defined in the top-level CMakeLists.txt (it also links the example plug-ins).
