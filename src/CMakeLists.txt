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




