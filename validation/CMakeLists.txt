# SIMC++ validation tools (standalone CMake project; P1.5).
#   cmake -S validation -B validation/build -DCMAKE_PREFIX_PATH=$ROOTSYS
#   cmake --build validation/build -j && (cd validation/build && ctest --output-on-failure)
# Kept independent of the top-level CMakeLists.txt on purpose (the core is built
# by another task); it can later be pulled in with add_subdirectory(validation).
cmake_minimum_required(VERSION 3.16)
project(simc_validation CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

if(DEFINED ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
endif()
find_package(ROOT 6.40 REQUIRED COMPONENTS Core RIO Tree Hist Gpad Graf MathCore)

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

add_executable(simc_compare compare/simc_compare.cpp compare/CompareCommon.cpp compare/CompareRunLevel.cpp
               compare/CompareTreeLevel.cpp compare/CompareReport.cpp)   # (P15.3 split)
target_link_libraries(simc_compare PRIVATE simcval ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
                      ROOT::Gpad ROOT::Graf ROOT::MathCore)
target_compile_options(simc_compare PRIVATE -Wall -Wextra)

# P8.2: event-by-event difference of two same-seed SIMC++ runs (fortran_compat vs corrected)
add_executable(simc_mode_diff compare/mode_diff.cpp)
target_link_libraries(simc_mode_diff PRIVATE ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_options(simc_mode_diff PRIVATE -Wall -Wextra)

# ---------------------------------------------------------------- self-tests
enable_testing()
set(SIMC_REFERENCE_DIR "$ENV{SIMC_REFERENCE}" CACHE PATH "Fortran reference data (empty = absent: self-tests not registered)")
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
  # the example summaries name their tree file as @SIMC_REFERENCE_DIR@/runs/...: configured into the build tree (D-073)
  foreach(ex eep_h_elastic_rad_A acc_phsp_hms_shms_A)
    configure_file(examples/${ex}.summary.json ${CMAKE_CURRENT_BINARY_DIR}/examples/${ex}.summary.json @ONLY)
  endforeach()
  add_test(NAME compare_json_vs_hist_identical
    COMMAND simc_compare --a ${R}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist
            --b ${CMAKE_CURRENT_BINARY_DIR}/examples/eep_h_elastic_rad_A.summary.json
            --out ${T}/json_vs_hist --config ${CFG} --no-pdf --expect IDENTICAL)
  set_tests_properties(compare_json_vs_hist_identical PROPERTIES DEPENDS compare_dump_summary)
  # P7.4: exactly one tree empty -> deterministic FAIL (the run-level numbers are identical here)
  add_test(NAME compare_one_tree_empty
    COMMAND simc_compare --a ${R}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist
            --b ${R}/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist --empty-tree-b
            --out ${T}/one_tree_empty --config ${CFG} --no-pdf --expect FAIL)
  # P12.3 (D-053): Fortran A vs B with the wtcontr sigma from the pass_cuts subset (SIDIS, pion) -> PASS, and the
  # subset reproduces each run's wtcontr (no "pass_cuts subset of the tree gives wtcontr" warning)
  foreach(tag sidis_h_k1_pip_rad pion_h_k1_rad_dec1)
    if(EXISTS ${R}/${tag}_B/${tag}_B.hist)
      add_test(NAME compare_FvsF_${tag}
        COMMAND simc_compare --a ${R}/${tag}_A/${tag}_A.hist --b ${R}/${tag}_B/${tag}_B.hist
                --label-a ${tag}_A --label-b ${tag}_B --out ${T}/FvsF_${tag} --config ${CFG} --no-pdf --expect PASS)
    endif()
  endforeach()
  foreach(tag eep_h_elastic_rad sidis_h_k1_pip_rad pion_h_k1_rad_dec1)
    if(EXISTS ${R}/${tag}_B/${tag}_B.hist)
      add_test(NAME compare_wtcontr_subset_${tag}
        COMMAND ${CMAKE_COMMAND} -DF=${T}/FvsF_${tag}/summary.txt
                -DWANT=A+pass_cuts+subset,+compat+window\;B+pass_cuts+subset,+compat+window
                -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_wtcontr_subset.cmake)
      set_tests_properties(compare_wtcontr_subset_${tag} PROPERTIES DEPENDS compare_FvsF_${tag})
    endif()
  endforeach()
else()
  message(STATUS "reference runs not found under ${R}: compare self-tests not registered")
endif()

set(SIMC_BASELINE_DIR "$ENV{SIMC_BASELINES}" CACHE PATH "SIMC++ baseline outputs (empty = absent)")
# P7.4: both trees empty (Ncontribute = 0 in both runs, K2 acceptance) -> per-variable tests skipped, PASS
# P10.4: regenerated on glibc 2.43 / gcc 16.2.1 (the original f9cea51_K2K4_s4242 is kept)
set(_K2C ${SIMC_BASELINE_DIR}/1aad1f0_glibc243_K2K4_s4242)
if(EXISTS ${R}/acc_phsp_K2_A/acc_phsp_K2_A.hist AND EXISTS ${_K2C}/cpp_acc_phsp_K2_A.summary.json)
  add_test(NAME compare_both_trees_empty
    COMMAND simc_compare --a ${R}/acc_phsp_K2_A --b ${_K2C}/cpp_acc_phsp_K2_A.summary.json
            --label-a F_K2_A --label-b cpp_K2_A --out ${T}/both_trees_empty --config ${CFG} --no-pdf --expect PASS)
endif()
# P12.3: a C++ compat run compared with the wrong mode (--mode-b corrected) is detected: its pass_cuts subset no longer
# reproduces its wtcontr (warning); with --mode-b compat it does. (The verdict is a deterministic FAIL: Ngen 20000 vs 100000.)
set(_SB ${SIMC_BASELINE_DIR}/1aad1f0_glibc243_sidis_h_k1_pip_rad_n20k_s4242)
if(EXISTS ${R}/sidis_h_k1_pip_rad_A/sidis_h_k1_pip_rad_A.hist AND EXISTS ${_SB}/dir_sidis_h_k1_pip_rad_s4242.summary.json)
  foreach(m compat corrected)
    add_test(NAME compare_wtcontr_mode_${m}
      COMMAND simc_compare --a ${R}/sidis_h_k1_pip_rad_A/sidis_h_k1_pip_rad_A.hist
              --b ${_SB}/dir_sidis_h_k1_pip_rad_s4242.summary.json --mode-b ${m}
              --out ${T}/wtcontr_mode_${m} --config ${CFG} --no-pdf --expect FAIL)   # Ngen 20000 vs 100000: det FAIL
  endforeach()
  add_test(NAME compare_wtcontr_mode_compat_consistent
    COMMAND ${CMAKE_COMMAND} -DF=${T}/wtcontr_mode_compat/summary.txt -DWANT=B+pass_cuts+subset,+compat+window
            -DFORBID=pass_cuts+subset+of+the+tree+gives -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_wtcontr_subset.cmake)
  add_test(NAME compare_wtcontr_mode_corrected_detected
    COMMAND ${CMAKE_COMMAND} -DF=${T}/wtcontr_mode_corrected/summary.txt -DWANT=B:+pass_cuts+subset+of+the+tree+gives+wtcontr
            -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_wtcontr_subset.cmake)
  set_tests_properties(compare_wtcontr_mode_compat_consistent PROPERTIES DEPENDS compare_wtcontr_mode_compat)
  set_tests_properties(compare_wtcontr_mode_corrected_detected PROPERTIES DEPENDS compare_wtcontr_mode_corrected)
endif()
# P15.1 (D-054): without --mode-b the pass_cuts mode comes from the summary's run.fortran_compat (the compat baseline: compat
# window); an explicit --mode-b that contradicts it is applied with a warning
if(EXISTS ${R}/sidis_h_k1_pip_rad_A/sidis_h_k1_pip_rad_A.hist AND EXISTS ${_SB}/dir_sidis_h_k1_pip_rad_s4242.summary.json)
  add_test(NAME compare_wtcontr_mode_auto
    COMMAND simc_compare --a ${R}/sidis_h_k1_pip_rad_A/sidis_h_k1_pip_rad_A.hist
            --b ${_SB}/dir_sidis_h_k1_pip_rad_s4242.summary.json
            --out ${T}/wtcontr_mode_auto --config ${CFG} --no-pdf --expect FAIL)   # Ngen 20000 vs 100000: det FAIL
  add_test(NAME compare_wtcontr_mode_auto_compat
    COMMAND ${CMAKE_COMMAND} -DF=${T}/wtcontr_mode_auto/summary.txt -DWANT=B+pass_cuts+subset,+compat+window
            -DFORBID=pass_cuts+subset+of+the+tree+gives;overrides+the+run+summary
            -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_wtcontr_subset.cmake)
  set_tests_properties(compare_wtcontr_mode_auto_compat PROPERTIES DEPENDS compare_wtcontr_mode_auto)
  add_test(NAME compare_wtcontr_mode_conflict_warned
    COMMAND ${CMAKE_COMMAND} -DF=${T}/wtcontr_mode_corrected/summary.txt
            -DWANT=--mode-b+corrected+overrides+the+run+summary+(compat)
            -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_wtcontr_subset.cmake)
  set_tests_properties(compare_wtcontr_mode_conflict_warned PROPERTIES DEPENDS compare_wtcontr_mode_corrected)
endif()
# P7.4: a Fortran tree with a NaN event (K2 norad, D-028 F2): reported, excluded, identical to itself
if(EXISTS ${R}/eep_h_K2_norad_A/eep_h_K2_norad_A.hist)
  add_test(NAME compare_nan_event_AvsA
    COMMAND simc_compare --a ${R}/eep_h_K2_norad_A/eep_h_K2_norad_A.hist --b ${R}/eep_h_K2_norad_A/eep_h_K2_norad_A.hist
            --out ${T}/nan_AvsA --config ${CFG} --no-pdf --expect IDENTICAL)
  add_test(NAME compare_nan_event_reported
    COMMAND grep -q "non-finite tree events (excluded from the per-variable tests): A=1 \\[.*\\] B=1 \\["
            ${T}/nan_AvsA/summary.txt)
  set_tests_properties(compare_nan_event_reported PROPERTIES DEPENDS compare_nan_event_AvsA)
endif()
