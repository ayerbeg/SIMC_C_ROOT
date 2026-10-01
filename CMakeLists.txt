cmake_minimum_required(VERSION 3.15)
project(SIMC_CPP VERSION 2.0 LANGUAGES CXX)

# SIMC++ (translation branch): docs/ARCHITECTURE.md, docs/DRIVER.md.
# The pre-translation ("legacy") code and its SIMC_BUILD_LEGACY option were removed in P7.2; they are in git history.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

enable_testing()

find_package(ROOT 6.20 REQUIRED COMPONENTS Core Tree Hist RIO)

# ============================================================================
# Find nlohmann/json
# ============================================================================

find_package(nlohmann_json 3.9.0 QUIET)
if(NOT nlohmann_json_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.2
    )
    FetchContent_MakeAvailable(json)
endif()

# ============================================================================
# New framework (translation branch, docs/ARCHITECTURE.md). ROOT-free core.
# ============================================================================

set(SIMC_REFERENCE_DIR "/home/ayerbe/NewDisk/reference" CACHE PATH
    "Fortran reference data (reference/runs, reference/functions); tests skip if absent")
set(SIMC_FORTRAN_DIR "/home/ayerbe/NewDisk/simc_gfortran" CACHE PATH
    "Fortran SIMC source tree (read-only; infiles/ used by tests)")

add_library(simcfw STATIC
    src/core/SimcParameters.cpp
    src/core/RunSetup.cpp
    src/input/InpReader.cpp
    src/input/EffectiveValues.cpp
    src/output/NtupleLayout.cpp
    src/output/SummaryJson.cpp
)
target_include_directories(simcfw PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(simcfw PUBLIC nlohmann_json::nlohmann_json)
target_compile_options(simcfw PRIVATE -Wall -Wextra -pedantic)

add_executable(inp2json tools/inp2json.cpp)
target_link_libraries(inp2json PRIVATE simcfw)
set_target_properties(inp2json PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(hist2json tools/hist2json.cpp)
target_link_libraries(hist2json PRIVATE simcfw)
set_target_properties(hist2json PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ROOT-based output sinks. OBJECT library so that the SIMC_REGISTER static
# registrations are always linked in (see core/Registry.h).
add_library(simcfw_root OBJECT
    src/output/RootTreeSink.cpp
    src/output/RunSummarySink.cpp
    src/output/GenHistSink.cpp
)
target_link_libraries(simcfw_root PUBLIC simcfw ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)
target_compile_options(simcfw_root PRIVATE -Wall -Wextra -pedantic)

# Physics-model plug-ins (ROOT-free). OBJECT library so the SIMC_REGISTER registrations are linked.
add_library(simcfw_physics OBJECT
    src/physics/proton/ProtonCrossSections.cpp
    src/physics/proton/ProtonModels.cpp
)
target_link_libraries(simcfw_physics PUBLIC simcfw)
target_compile_options(simcfw_physics PRIVATE -Wall -Wextra -pedantic)

# IRandom implementations and gauss1 (src/core/random; P3.1, own target since P2.2). STATIC: Mt64Random.cpp has a
# SIMC_REGISTER(IRandom, "mt64") that is only linked when its object file is pulled in (SimcRun references
# Mt64Random directly for this reason, P7.2).
add_library(simcfw_random STATIC
    src/core/random/Gauss1.cpp
    src/core/random/Mt64Random.cpp
    src/core/random/ScriptedRandom.cpp
)
target_link_libraries(simcfw_random PUBLIC simcfw)
target_compile_options(simcfw_random PRIVATE -Wall -Wextra -pedantic)

# Phase 2/3 modules: each has its own CMakeLists so parallel tasks don't edit this file (director, P2).
#   src/setup    : RunSetup from init.f / calculate_central (target simcfw_setup)
#   src/optics   : shared COSY transport engine (target simcfw_optics)
#   src/material : enerloss_new, musc, musc_ext, target_musc, gauss1, trip_thru_target (target simcfw_material)
add_subdirectory(src/material)
add_subdirectory(src/setup)
add_subdirectory(src/optics)      # shared COSY engine: transp, project, rotate_*axis, locforunt (P6.1)
add_subdirectory(src/spectrometers/hms)   # ISpectrometer "hms" (P6.2)
add_subdirectory(src/spectrometers/shms)  # ISpectrometer "shms" (P6.3)
add_subdirectory(src/radiation)           # IRadiator "simc_rad": radc.f, brem.f, radc_init_ev (P5.1)
add_subdirectory(src/reactions)           # IBeam "standard", IReaction "hydrogen_elastic"/"phase_space": event.f (P4.1)
add_subdirectory(src/driver)              # simc.f: montecarlo (P7.1), main loop + normalisation (P7.2)

add_subdirectory(tests/unit)
add_subdirectory(tests/l1)
add_subdirectory(tests/regression)

# ============================================================================
# Installation
# ============================================================================
install(TARGETS simcpp inp2json hist2json RUNTIME DESTINATION bin)
install(DIRECTORY include/simc DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY data/ DESTINATION share/simc/data)
