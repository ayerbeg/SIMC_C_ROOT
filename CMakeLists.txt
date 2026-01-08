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

# Copy COSY matrices
add_custom_command(TARGET copy_data_files POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data/matrices
        ${CMAKE_BINARY_DIR}/data/matrices
    COMMENT "Copying COSY matrix files"
)

# Copy configuration files
add_custom_command(TARGET copy_data_files POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data/config
        ${CMAKE_BINARY_DIR}/data/config
    COMMENT "Copying configuration files"
)

# Copy physics tables
add_custom_command(TARGET copy_data_files POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/data/tables
        ${CMAKE_BINARY_DIR}/data/tables
    COMMENT "Copying physics tables"
)

# ============================================================================
# Find ROOT
# ============================================================================

find_package(ROOT 6.20 REQUIRED COMPONENTS Core Tree Hist RIO)
if(ROOT_FOUND)
    message(STATUS "Found ROOT ${ROOT_VERSION} at ${ROOT_USE_FILE}")
    include(${ROOT_USE_FILE})
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

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

# ============================================================================
# Add subdirectories
# ============================================================================

add_subdirectory(src)

# ============================================================================
# Tests (optional)
# ============================================================================

option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    if(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
        enable_testing()
        add_subdirectory(tests)
    else()
        message(STATUS "tests/CMakeLists.txt not found - skipping tests")
    endif()
endif()

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
