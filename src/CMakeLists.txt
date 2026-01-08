# src/CMakeLists.txt

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

# Installation
install(TARGETS simc simc_core simc_io
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
