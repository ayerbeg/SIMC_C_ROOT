# src/CMakeLists.txt
# PHASE 5c.2 - SHMS with symlink-based data file handling
# VERSION: 2.0 - Compatible with CMake 3.10+

# Collect source files
set(CORE_SOURCES
    core/SimcEvent.cpp
    core/ConfigManager.cpp
    core/RandomGenerator.cpp
)

set(IO_SOURCES
    io/OutputManager.cpp
)

# Physics modules - collect all .cpp files in physics/
file(GLOB PHYSICS_SOURCES "physics/*.cpp")

message(STATUS "Physics modules found:")
foreach(src ${PHYSICS_SOURCES})
    get_filename_component(filename ${src} NAME)
    message(STATUS "  - ${filename}")
endforeach()

# ============================================================================
# PHASE 5c.1-5c.2: Transport and Spectrometer modules
# ============================================================================

# Transport modules
file(GLOB TRANSPORT_SOURCES "transport/*.cpp")

if(TRANSPORT_SOURCES)
    message(STATUS "Transport modules found:")
    foreach(src ${TRANSPORT_SOURCES})
        get_filename_component(filename ${src} NAME)
        message(STATUS "  - ${filename}")
    endforeach()
endif()

# PHASE 5c.2: Spectrometer modules
file(GLOB SPECTROMETER_SOURCES "spectrometers/*.cpp")

message(STATUS "Spectrometer modules found:")
foreach(src ${SPECTROMETER_SOURCES})
    get_filename_component(filename ${src} NAME)
    message(STATUS "  - ${filename}")
endforeach()

# Generate ROOT dictionary for SimcEvent
set(DICT_HEADERS 
    ${CMAKE_SOURCE_DIR}/include/simc/core/SimcEvent.h
)

ROOT_GENERATE_DICTIONARY(G__SimcEvent 
    ${DICT_HEADERS}
    LINKDEF ${CMAKE_SOURCE_DIR}/include/simc/core/SimcEventLinkDef.h
)

# Core library
add_library(simc_core SHARED
    ${CORE_SOURCES}
    ${PHYSICS_SOURCES}
    ${TRANSPORT_SOURCES}
    ${SPECTROMETER_SOURCES}
    G__SimcEvent.cxx
)

target_link_libraries(simc_core
    PUBLIC
        ${ROOT_LIBRARIES}
        nlohmann_json::nlohmann_json
)

target_include_directories(simc_core
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${ROOT_INCLUDE_DIRS}
)

# I/O library
add_library(simc_io SHARED
    ${IO_SOURCES}
)

target_link_libraries(simc_io
    PUBLIC
        simc_core
        ${ROOT_LIBRARIES}
)

# Main executable
add_executable(simc main.cpp)

target_link_libraries(simc
    PRIVATE
        simc_core
        simc_io
)

# ============================================================================
# PHASE 5c.2: Handle spectrometer matrix data files with symlinks
# ============================================================================

# SHMS data directory
set(SHMS_DATA_DIR "${CMAKE_SOURCE_DIR}/data/shms")

# Symlink names (Fortran-style short names)
set(SHMS_FORWARD_LINK "shms_forward.dat")
set(SHMS_RECON_LINK "shms_recon.dat")

# Check if data directory exists
if(EXISTS ${SHMS_DATA_DIR})
    message(STATUS "Found SHMS data directory: ${SHMS_DATA_DIR}")
    
    # Look for symlinks in data/shms/
    set(forward_symlink "${SHMS_DATA_DIR}/${SHMS_FORWARD_LINK}")
    set(recon_symlink "${SHMS_DATA_DIR}/${SHMS_RECON_LINK}")
    
    # Check if symlinks exist
    if(EXISTS ${forward_symlink})
        # Get the target of the symlink
        get_filename_component(forward_target ${forward_symlink} REALPATH)
        get_filename_component(forward_name ${forward_target} NAME)
        message(STATUS "  Forward matrix: ${SHMS_FORWARD_LINK} → ${forward_name}")
        
        # Create symlink in build directory pointing to the actual file
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${forward_symlink}
                ${CMAKE_BINARY_DIR}/${SHMS_FORWARD_LINK}
        )
    else()
        message(WARNING "  Forward matrix symlink not found: ${forward_symlink}")
        message(WARNING "  Please create symlink in ${SHMS_DATA_DIR}/")
        message(WARNING "  Example:")
        message(WARNING "    cd ${SHMS_DATA_DIR}")
        message(WARNING "    ln -s shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_forward.dat ${SHMS_FORWARD_LINK}")
    endif()
    
    if(EXISTS ${recon_symlink})
        # Get the target of the symlink
        get_filename_component(recon_target ${recon_symlink} REALPATH)
        get_filename_component(recon_name ${recon_target} NAME)
        message(STATUS "  Recon matrix: ${SHMS_RECON_LINK} → ${recon_name}")
        
        # Create symlink in build directory pointing to the actual file
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${recon_symlink}
                ${CMAKE_BINARY_DIR}/${SHMS_RECON_LINK}
        )
    else()
        message(WARNING "  Recon matrix symlink not found: ${recon_symlink}")
        message(WARNING "  Please create symlink in ${SHMS_DATA_DIR}/")
        message(WARNING "  Example:")
        message(WARNING "    cd ${SHMS_DATA_DIR}")
        message(WARNING "    ln -s shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_recon_fit.dat ${SHMS_RECON_LINK}")
    endif()
    
else()
    message(WARNING "SHMS data directory not found: ${SHMS_DATA_DIR}")
    message(WARNING "Please create directory structure:")
    message(WARNING "  mkdir -p ${SHMS_DATA_DIR}")
    message(WARNING "  cd ${SHMS_DATA_DIR}")
    message(WARNING "  # Download matrix files, then create symlinks:")
    message(WARNING "  ln -s shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_forward.dat ${SHMS_FORWARD_LINK}")
    message(WARNING "  ln -s shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_recon_fit.dat ${SHMS_RECON_LINK}")
endif()

# ============================================================================
# Optional: Add custom target to help set up symlinks
# ============================================================================

# Create a helper script to set up data files
file(WRITE "${CMAKE_BINARY_DIR}/setup_shms_data.sh"
"#!/bin/bash
# Helper script to set up SHMS matrix data files

SHMS_DATA_DIR=\"${SHMS_DATA_DIR}\"

echo \"Setting up SHMS matrix data files...\"
echo \"\"

# Create directory if it doesn't exist
if [ ! -d \"\$SHMS_DATA_DIR\" ]; then
    echo \"Creating directory: \$SHMS_DATA_DIR\"
    mkdir -p \"\$SHMS_DATA_DIR\"
fi

cd \"\$SHMS_DATA_DIR\"

# Download files if they don't exist
FORWARD_FILE=\"shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_forward.dat\"
RECON_FILE=\"shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_recon_fit.dat\"

if [ ! -f \"\$FORWARD_FILE\" ]; then
    echo \"Downloading forward matrix...\"
    wget https://raw.githubusercontent.com/ayerbeg/simc_gfortran/master/shms/\$FORWARD_FILE
else
    echo \"Forward matrix already exists: \$FORWARD_FILE\"
fi

if [ ! -f \"\$RECON_FILE\" ]; then
    echo \"Downloading reconstruction matrix...\"
    wget https://raw.githubusercontent.com/ayerbeg/simc_gfortran/master/shms/\$RECON_FILE
else
    echo \"Recon matrix already exists: \$RECON_FILE\"
fi

# Create symlinks
echo \"\"
echo \"Creating symlinks...\"

if [ -f \"${SHMS_FORWARD_LINK}\" ] || [ -L \"${SHMS_FORWARD_LINK}\" ]; then
    rm \"${SHMS_FORWARD_LINK}\"
fi
ln -s \"\$FORWARD_FILE\" \"${SHMS_FORWARD_LINK}\"
echo \"  Created: ${SHMS_FORWARD_LINK} -> \$FORWARD_FILE\"

if [ -f \"${SHMS_RECON_LINK}\" ] || [ -L \"${SHMS_RECON_LINK}\" ]; then
    rm \"${SHMS_RECON_LINK}\"
fi
ln -s \"\$RECON_FILE\" \"${SHMS_RECON_LINK}\"
echo \"  Created: ${SHMS_RECON_LINK} -> \$RECON_FILE\"

echo \"\"
echo \"Done! Now re-run cmake to pick up the symlinks.\"
")

# Make the script executable
execute_process(
    COMMAND chmod +x "${CMAKE_BINARY_DIR}/setup_shms_data.sh"
)

message(STATUS "")
message(STATUS "SHMS Data Setup:")
message(STATUS "  If matrix files are missing, run: ${CMAKE_BINARY_DIR}/setup_shms_data.sh")
message(STATUS "")

# Installation
install(TARGETS simc simc_core simc_io
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# ============================================================================
# Install matrix files (CMake 3.10+ compatible)
# Resolves symlinks and installs the actual files with short names
# ============================================================================

if(EXISTS "${SHMS_DATA_DIR}/${SHMS_FORWARD_LINK}")
    # Resolve the symlink to get the actual file path
    get_filename_component(forward_real "${SHMS_DATA_DIR}/${SHMS_FORWARD_LINK}" REALPATH)
    if(EXISTS ${forward_real})
        # Install the actual file, but rename it to the short symlink name
        install(FILES ${forward_real}
            DESTINATION bin
            RENAME ${SHMS_FORWARD_LINK}
        )
    endif()
endif()

if(EXISTS "${SHMS_DATA_DIR}/${SHMS_RECON_LINK}")
    # Resolve the symlink to get the actual file path
    get_filename_component(recon_real "${SHMS_DATA_DIR}/${SHMS_RECON_LINK}" REALPATH)
    if(EXISTS ${recon_real})
        # Install the actual file, but rename it to the short symlink name
        install(FILES ${recon_real}
            DESTINATION bin
            RENAME ${SHMS_RECON_LINK}
        )
    endif()
endif()
