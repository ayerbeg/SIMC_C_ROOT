# examples/out_of_tree_plugin — a plug-in built OUTSIDE the SIMC++ source tree against the installed package (P15.3).
#   cmake --install <simc build> --prefix <prefix>
#   cmake -S examples/out_of_tree_plugin -B build -DCMAKE_PREFIX_PATH="<prefix>;<ROOT install>"
#   cmake --build build && build/my_simc <file.inp> --pipeline toy_window.json -o out
# Mechanism (docs/EXTENDING.md): plug-ins register themselves at static-initialisation time (SIMC_REGISTER); they are
# LINKED into an executable, not loaded at run time. my_simc = the simc++ command line (simc::driver::CommandLineMain)
# + this plug-in; SIMC::simc brings the built-in stages.
cmake_minimum_required(VERSION 3.16)
project(simc_out_of_tree_plugin CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SIMC 2.0 REQUIRED)

add_executable(my_simc main.cpp ToyWindowSpectrometer.cpp)
target_link_libraries(my_simc PRIVATE SIMC::simc)
target_compile_options(my_simc PRIVATE -Wall -Wextra -pedantic)
