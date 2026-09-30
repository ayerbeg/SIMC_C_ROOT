cmake_minimum_required(VERSION 3.15)
project(SIMC_CPP VERSION 2.0 LANGUAGES CXX)

# C++17 standard required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

# ============================================================================
# Data Directory Configuration
# ============================================================================

# Option 1: Set data directory at configure time
set(SIMC_DATA_DIR "${CMAKE_SOURCE_DIR}/data" CACHE PATH 
    "Path to SIMC data files (COSY matrices, config, tables)")

# Option 2: Install data files and set install path
set(SIMC_INSTALL_DATADIR "${CMAKE_INSTALL_PREFIX}/share/simc/data")

# Create a configuration header with data path
configure_file(
    "${CMAKE_SOURCE_DIR}/include/simc/SimcConfig.h.in"
    "${CMAKE_BINARY_DIR}/include/simc/SimcConfig.h"
    @ONLY
)

# Add the build directory to include path for generated config
include_directories(${CMAKE_BINARY_DIR}/include)

# ============================================================================
# Copy data files to build directory (for development)
# ============================================================================

# This makes development easier - data files available in build/
add_custom_target(copy_data_files ALL
    COMMENT "Copying data files to build directory"
)

# Copy all data files
add_custom_command(TARGET copy_data_files POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data
        ${CMAKE_BINARY_DIR}/data
    COMMENT "Copying data files to build directory"
)

# Copy COSY matrices
#add_custom_command(TARGET copy_data_files POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_directory
#        ${CMAKE_SOURCE_DIR}/data/matrices
#        ${CMAKE_BINARY_DIR}/data/matrices
#    COMMENT "Copying COSY matrix files"
#)

# Copy configuration files
#add_custom_command(TARGET copy_data_files POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_directory
#        ${CMAKE_SOURCE_DIR}/data/config
#        ${CMAKE_BINARY_DIR}/data/config
#    COMMENT "Copying configuration files"
#)

# Copy physics tables
#add_custom_command(TARGET copy_data_files POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_directory
#        ${CMAKE_SOURCE_DIR}/data/tables
#        ${CMAKE_BINARY_DIR}/data/tables
#    COMMENT "Copying physics tables"
#)

# ============================================================================
# Find ROOT
# ============================================================================

# SIMC_BUILD_LEGACY: the pre-translation code (simc, simc_core, simc_io, old tests).
# It keeps building for now; nothing in the new framework may depend on it (P1.3).
# P1.7: still ON by default because tools/run_simc.sh runs the legacy `simc` binary; the new
# framework and all tests/unit + tests/l1 build and pass with -DSIMC_BUILD_LEGACY=OFF.
option(SIMC_BUILD_LEGACY "Build the legacy (pre-translation) targets" ON)

enable_testing()

find_package(ROOT 6.20 REQUIRED COMPONENTS Core Tree Hist RIO)
if(SIMC_BUILD_LEGACY)
if(ROOT_FOUND)
    message(STATUS "Found ROOT ${ROOT_VERSION} at ${ROOT_USE_FILE}")
    include(${ROOT_USE_FILE})
endif()
endif()

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
# Include directories
# ============================================================================

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

# Phase 2/3 modules: each has its own CMakeLists so parallel tasks don't edit this file (director, P2).
#   src/setup    : RunSetup from init.f / calculate_central (target simcfw_setup)
#   src/optics   : shared COSY transport engine (target simcfw_optics)
#   src/material : enerloss_new, musc, musc_ext, target_musc, gauss1, trip_thru_target (target simcfw_material)
add_subdirectory(src/material)
add_subdirectory(src/setup)
add_subdirectory(src/optics)      # shared COSY engine: transp, project, rotate_*axis, locforunt (P6.1)

add_subdirectory(tests/unit)
add_subdirectory(tests/l1)

if(SIMC_BUILD_LEGACY)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

# ============================================================================
# Add subdirectories (legacy)
# ============================================================================

add_subdirectory(src)

# ============================================================================
# Phase 5d: Link all spectrometer sources
# ============================================================================

message(STATUS "Linking Phase 5d spectrometer sources...")

target_sources(simc PRIVATE
    # Phase 5d: All 5 spectrometers
    src/spectrometers/SHMS.cpp
    src/spectrometers/HMS.cpp
    src/spectrometers/SOS.cpp
    src/spectrometers/HRSr.cpp
    src/spectrometers/HRSl.cpp
)

set_target_properties(simc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

# ============================================================================
# Tests (optional)
# ============================================================================

option(BUILD_TESTS "Build legacy unit tests" ON)
if(BUILD_TESTS)
    if(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
        add_subdirectory(tests)
    else()
        message(STATUS "tests/CMakeLists.txt not found - skipping tests")
    endif()
endif()
endif()  # SIMC_BUILD_LEGACY

# ============================================================================
# Installation
# ============================================================================

# Install headers
install(DIRECTORY include/simc 
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")

# Install data files
install(DIRECTORY data/
        DESTINATION ${SIMC_INSTALL_DATADIR}
        FILES_MATCHING 
        PATTERN "*.dat"
        PATTERN "*.json"
        PATTERN "*.txt")

# Install executables (will be added by src/CMakeLists.txt)
# install(TARGETS simc DESTINATION bin)

# ============================================================================
# Create uninstall target
# ============================================================================

#if(NOT TARGET uninstall)
#    configure_file(
#        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
#        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
#        IMMEDIATE @ONLY)

#    add_custom_target(uninstall
#        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
#endif()

# ============================================================================
# Documentation
# ============================================================================

option(BUILD_DOCUMENTATION "Build documentation with Doxygen" OFF)
if(BUILD_DOCUMENTATION)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        configure_file(${CMAKE_SOURCE_DIR}/docs/Doxyfile.in 
                      ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
        add_custom_target(doc ALL
            ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
        )
    endif()
endif()

# ============================================================================
# Print summary
# ============================================================================

message(STATUS "")
message(STATUS "========================================")
message(STATUS "SIMC C++/ROOT Configuration Summary")
message(STATUS "========================================")
message(STATUS "Build type:       ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler:     ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ standard:     C++${CMAKE_CXX_STANDARD}")
message(STATUS "C++ flags:        ${CMAKE_CXX_FLAGS}")
message(STATUS "ROOT version:     ${ROOT_VERSION}")
message(STATUS "ROOT libraries:   ${ROOT_LIBRARIES}")
message(STATUS "Install prefix:   ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Data directory:   ${SIMC_DATA_DIR}")
message(STATUS "Build tests:      ${BUILD_TESTS}")
message(STATUS "========================================")
message(STATUS "")
message(STATUS "Data files will be copied to: ${CMAKE_BINARY_DIR}/data/")
message(STATUS "After 'make install', data will be in: ${SIMC_INSTALL_DATADIR}")
message(STATUS "")

