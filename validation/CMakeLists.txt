# SIMC++ validation tools (standalone CMake project; P1.5).
#   cmake -S validation -B validation/build -DCMAKE_PREFIX_PATH=/home/ayerbe/Builds/ROOT6_install
#   cmake --build validation/build -j && (cd validation/build && ctest --output-on-failure)
# Kept independent of the top-level CMakeLists.txt on purpose (the core is built
# by another task); it can later be pulled in with add_subdirectory(validation).
cmake_minimum_required(VERSION 3.16)
project(simc_validation CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

list(APPEND CMAKE_PREFIX_PATH /home/ayerbe/Builds/ROOT6_install)
find_package(ROOT 6.20 REQUIRED COMPONENTS Core RIO Tree Hist Gpad Graf MathCore)

find_package(nlohmann_json 3.9 QUIET)
if(NOT nlohmann_json_FOUND)
  find_path(NLOHMANN_INCLUDE nlohmann/json.hpp PATHS /usr/local/include /usr/include REQUIRED)
  add_library(nlohmann_json_hdr INTERFACE)
  target_include_directories(nlohmann_json_hdr INTERFACE ${NLOHMANN_INCLUDE})
  add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json_hdr)
endif()

add_library(simcval STATIC compare/RunSummary.cpp compare/Stats.cpp)
target_include_directories(simcval PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/compare)
target_link_libraries(simcval PUBLIC ROOT::Core ROOT::MathCore nlohmann_json::nlohmann_json)
target_compile_options(simcval PRIVATE -Wall -Wextra)

add_executable(simc_compare compare/simc_compare.cpp)
target_link_libraries(simc_compare PRIVATE simcval ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
                      ROOT::Gpad ROOT::Graf ROOT::MathCore)
target_compile_options(simc_compare PRIVATE -Wall -Wextra)

# ---------------------------------------------------------------- self-tests
enable_testing()
set(SIMC_REFERENCE_DIR "/home/ayerbe/NewDisk/reference" CACHE PATH "Fortran reference data")
set(R ${SIMC_REFERENCE_DIR}/runs)
set(T ${CMAKE_CURRENT_BINARY_DIR}/selftest)
set(CFG ${CMAKE_CURRENT_SOURCE_DIR}/config/compare_default.json)

if(EXISTS ${R}/eep_h_elastic_norad_A/eep_h_elastic_norad_A.hist)
  # 1. A vs A (same file) -> all p = 1, pulls = 0
  foreach(tag eep_h_elastic_norad eep_h_elastic_rad acc_phsp_hms_shms)
    add_test(NAME compare_AvsA_${tag}
      COMMAND simc_compare --a ${R}/${tag}_A/${tag}_A.hist --b ${R}/${tag}_A/${tag}_A.hist
              --out ${T}/AvsA_${tag} --config ${CFG} --no-pdf --expect IDENTICAL)
  endforeach()
  # 2. Fortran A vs B (independent seeds) -> PASS
  foreach(tag eep_h_elastic_norad eep_h_elastic_rad acc_phsp_hms_shms)
    add_test(NAME compare_FvsF_${tag}
      COMMAND simc_compare --a ${R}/${tag}_A/${tag}_A.hist --b ${R}/${tag}_B/${tag}_B.hist
              --label-a ${tag}_A --label-b ${tag}_B
              --out ${T}/FvsF_${tag} --config ${CFG} --expect PASS)
  endforeach()
  # 3a. different physics -> FAIL
  add_test(NAME compare_sensitivity_rad_vs_norad
    COMMAND simc_compare --a ${R}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist
            --b ${R}/eep_h_elastic_norad_A/eep_h_elastic_norad_A.hist
            --label-a rad_A --label-b norad_A --out ${T}/rad_vs_norad --config ${CFG} --expect FAIL)
  # 3b. HMS delta shifted by +0.1 % -> must flag hsdelta
  add_test(NAME compare_sensitivity_hsdelta_shift
    COMMAND simc_compare --a ${R}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist
            --b ${R}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist --perturb hsdelta:0.1
            --label-a rad_A --label-b "rad_A(hsdelta+0.1)" --out ${T}/hsdelta_shift --config ${CFG}
            --expect FAIL --expect-flagged hsdelta)
  # parser: .hist -> summary JSON round trip
  add_test(NAME compare_dump_summary
    COMMAND simc_compare --dump-summary ${R}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist
            ${T}/eep_h_elastic_rad_A.summary.json)
  add_test(NAME compare_json_vs_hist_identical
    COMMAND simc_compare --a ${R}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist
            --b ${CMAKE_CURRENT_SOURCE_DIR}/examples/eep_h_elastic_rad_A.summary.json
            --out ${T}/json_vs_hist --config ${CFG} --no-pdf --expect IDENTICAL)
  set_tests_properties(compare_json_vs_hist_identical PROPERTIES DEPENDS compare_dump_summary)
else()
  message(STATUS "reference runs not found under ${R}: compare self-tests not registered")
endif()
