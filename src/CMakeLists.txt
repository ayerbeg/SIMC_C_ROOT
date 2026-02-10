# src/CMakeLists.txt
# MODIFIED FOR PHASE 5c.2 - Added SHMS spectrometer transport

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

# Debug: Print what physics files were found
message(STATUS "Physics modules found:")
foreach(src ${PHYSICS_SOURCES})
    get_filename_component(filename ${src} NAME)
    message(STATUS "  - ${filename}")
endforeach()

# ============================================================================
# PHASE 5c.1-5c.2: Transport and Spectrometer modules
# ============================================================================

# Transport modules (if you have them in transport/ directory)
file(GLOB TRANSPORT_SOURCES "transport/*.cpp")

# Debug: Print what transport files were found
if(TRANSPORT_SOURCES)
    message(STATUS "Transport modules found:")
    foreach(src ${TRANSPORT_SOURCES})
        get_filename_component(filename ${src} NAME)
        message(STATUS "  - ${filename}")
    endforeach()
endif()

# PHASE 5c.2: Spectrometer modules - collect all .cpp files in spectrometers/
file(GLOB SPECTROMETER_SOURCES "spectrometers/*.cpp")

# Debug: Print what spectrometer files were found
message(STATUS "Spectrometer modules found:")
foreach(src ${SPECTROMETER_SOURCES})
    get_filename_component(filename ${src} NAME)
    message(STATUS "  - ${filename}")
endforeach()

# Generate ROOT dictionary for SimcEvent
set(DICT_HEADERS 
    ${CMAKE_SOURCE_DIR}/include/simc/SimcEvent.h
)

ROOT_GENERATE_DICTIONARY(G__SimcEvent 
    ${DICT_HEADERS}
    LINKDEF ${CMAKE_SOURCE_DIR}/include/simc/SimcEventLinkDef.h
)

# Core library
add_library(simc_core SHARED
    ${CORE_SOURCES}
    ${PHYSICS_SOURCES}
    ${TRANSPORT_SOURCES}      # PHASE 5c.1: Added transport sources
    ${SPECTROMETER_SOURCES}   # PHASE 5c.2: Added spectrometer sources
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
# PHASE 5c.2: Copy SHMS matrix data files to build directory
# ============================================================================

# Define matrix files to copy
set(SHMS_MATRIX_FILES
    shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_forward.dat
    shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_recon_fit.dat
)

# Copy each matrix file to build directory
foreach(matrix_file ${SHMS_MATRIX_FILES})
    if(EXISTS ${CMAKE_SOURCE_DIR}/${matrix_file})
        configure_file(
            ${CMAKE_SOURCE_DIR}/${matrix_file}
            ${CMAKE_BINARY_DIR}/${matrix_file}
            COPYONLY
        )
        message(STATUS "Copying SHMS matrix file: ${matrix_file}")
    else()
        message(WARNING "SHMS matrix file not found: ${matrix_file}")
        message(WARNING "  Please download from: https://raw.githubusercontent.com/ayerbeg/simc_gfortran/master/shms/${matrix_file}")
    endif()
endforeach()

# Installation
install(TARGETS simc simc_core simc_io
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Install matrix files
install(FILES ${SHMS_MATRIX_FILES}
    DESTINATION bin
    OPTIONAL
)
