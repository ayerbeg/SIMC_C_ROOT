# tests/regression — before/after output identity of complete runs (task P7.4).
# A refactor must not change a single number: the outputs of `simc++ <inp> --seed 4242` must equal the files that
# the reference commit produced (SIMC_BASELINE_DIR/<commit>_<tag>_s4242/, see the README.txt there).
# P7.5 (D-031): regression_identity_<tag> runs with --fortran-compat against the f9cea51 baselines (exact Fortran
# behaviour, the outputs before the corrections existed).
# Label "regression"; the tests are skipped (77) when the baseline is absent.
add_executable(simc_root_identity root_identity.cpp)
target_link_libraries(simc_root_identity PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)
set_target_properties(simc_root_identity PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(SIMC_BASELINE_DIR "/home/ayerbe/NewDisk/SIMC++_Project/baselines" CACHE PATH
    "Reference outputs of earlier commits for the output-identity tests")
set(SIMC_BASELINE_COMMIT "f9cea51" CACHE STRING "Commit whose outputs are the identity baseline")

foreach(_p "rad;eep_h_elastic_rad_A" "norad;eep_h_elastic_norad_A" "phsp;acc_phsp_hms_shms_A")
    list(GET _p 0 _tag)
    list(GET _p 1 _inp)
    add_test(NAME regression_identity_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${SIMC_FORTRAN_DIR}/infiles/${_inp}.inp 4242 dir_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_BASELINE_COMMIT}_${_tag}_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/identity_${_tag} --fortran-compat)
    set_tests_properties(regression_identity_${_tag} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 1200
                         LABELS regression)
endforeach()

# P7.5 (D-031): the DEFAULT (corrected) mode against its own baselines, produced by commit 08812e6
# (SIMC_BASELINE_DIR/<commit>_corrected_<tag>_s4242/, README.txt there lists the differences to the f9cea51 ones).
set(SIMC_CORRECTED_BASELINE_COMMIT "08812e6" CACHE STRING "Commit whose corrected-mode outputs are the identity baseline")
foreach(_p "rad;eep_h_elastic_rad_A" "norad;eep_h_elastic_norad_A" "phsp;acc_phsp_hms_shms_A")
    list(GET _p 0 _tag)
    list(GET _p 1 _inp)
    add_test(NAME regression_corrected_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${SIMC_FORTRAN_DIR}/infiles/${_inp}.inp 4242 dir_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_CORRECTED_BASELINE_COMMIT}_corrected_${_tag}_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/corrected_${_tag})
    set_tests_properties(regression_corrected_${_tag} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 1200
                         LABELS regression)
endforeach()

# P8.2: D(e,e'p) and Fe-56(e,e'p) (rad, seed 4242, 20000 events) against the baselines of commit 508633a, in both modes
# (SIMC_BASELINE_DIR/508633a_[corrected_]eep_<tag>_n20k_s4242/). Inputs: the P8.F1 L2 inputs in reference/runs.
set(SIMC_EEP_BASELINE_COMMIT "508633a" CACHE STRING "Commit whose D/A(e,e'p) outputs are the identity baseline")
foreach(_tag d_rad fe56_rad)
    set(_inp ${SIMC_REFERENCE_DIR}/runs/eep_${_tag}_A/eep_${_tag}_A.inp)
    add_test(NAME regression_identity_eep_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_eep_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_EEP_BASELINE_COMMIT}_eep_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/identity_eep_${_tag} --fortran-compat --ngen 20000)
    add_test(NAME regression_corrected_eep_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_eep_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_EEP_BASELINE_COMMIT}_corrected_eep_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/corrected_eep_${_tag} --ngen 20000)
    set_tests_properties(regression_identity_eep_${_tag} regression_corrected_eep_${_tag} PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 1200 LABELS regression)
endforeach()
