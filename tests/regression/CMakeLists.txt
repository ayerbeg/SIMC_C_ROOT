# tests/regression — before/after output identity of complete runs (task P7.4).
# A refactor must not change a single number: the outputs of `simc++ <inp> --seed 4242` must equal the files that
# the reference commit produced (SIMC_BASELINE_DIR/<commit>_<tag>_s4242/, see the README.txt there).
# P7.5 (D-031): regression_identity_<tag> runs with --fortran-compat against the f9cea51 baselines (exact Fortran
# behaviour, the outputs before the corrections existed).
# Label "regression"; the tests are skipped (77) when the baseline is absent.
# P10.4 (D-044): the baselines are bitwise valid for one toolchain only. All five families were regenerated with
# glibc 2.43 / gcc 16.2.1 by commit 4408493 (validation/scripts/regen_baselines.sh) into NEW folders
# <commit>_glibc243_<tag>; the *_COMMIT variables below therefore hold "4408493_glibc243". The original folders
# (f9cea51_*, 08812e6_corrected_*, 508633a_*, 3b27d5c_*, 1275770_*) are kept; their README.txt explain the history.
# P15.1 (D-054): every family was regenerated ONCE by commit 55cf7a0 (run.fortran_compat and run.rng added to the summary
# JSON/text); validation/scripts/compare_baseline_sets.sh with p151_baseline_map.txt proves that the 27 new folders
# 55cf7a0_glibc243_* differ from the previous ones (4408493, 51e86c2, 31b4232, e6feaf0, 17310b4) only by those keys.
# The *_COMMIT variables below hold "55cf7a0_glibc243"; the earlier folders are kept.
# On a failure, output_identity.sh prints the baseline toolchain (README.txt "toolchain:") next to this machine's.
# P10.5: the baselines were made with ROOT 6.30.06; under ROOT 6.40.04 (C++20) every regression test passes unchanged
# (contents identical; only ROOT metadata in the .root bytes differ, which simc_root_identity ignores), so they stay.
add_executable(simc_root_identity root_identity.cpp)
target_link_libraries(simc_root_identity PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)
set_target_properties(simc_root_identity PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(SIMC_BASELINE_DIR "$ENV{SIMC_BASELINES}" CACHE PATH
    "Reference outputs of earlier commits for the output-identity tests")
if(NOT SIMC_BASELINE_DIR)   # empty: a non-existent placeholder, the identity tests skip (77)
    set(SIMC_BASELINE_DIR "${CMAKE_BINARY_DIR}/not_configured_SIMC_BASELINE_DIR")
endif()
set(SIMC_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose outputs are the identity baseline")

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
set(SIMC_CORRECTED_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose corrected-mode outputs are the identity baseline")
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
set(SIMC_EEP_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose D/A(e,e'p) outputs are the identity baseline")
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

# P9.2: exclusive pion production (pi+ H with decay in flight, pi- D; K1, rad, seed 4242, 20000 events) against the
# baselines of commit 3b27d5c, in both modes (SIMC_BASELINE_DIR/3b27d5c_[corrected_]pion_<tag>_n20k_s4242/). Inputs: the
# P9.F1 L2 inputs in reference/runs.
set(SIMC_PION_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose pion outputs are the identity baseline")
foreach(_tag h_k1_rad_dec1 d_k1_pim_rad_dec1)
    set(_inp ${SIMC_REFERENCE_DIR}/runs/pion_${_tag}_A/pion_${_tag}_A.inp)
    add_test(NAME regression_identity_pion_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_pion_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_PION_BASELINE_COMMIT}_pion_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/identity_pion_${_tag} --fortran-compat --ngen 20000)
    add_test(NAME regression_corrected_pion_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_pion_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_PION_BASELINE_COMMIT}_corrected_pion_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/corrected_pion_${_tag} --ngen 20000)
    set_tests_properties(regression_identity_pion_${_tag} regression_corrected_pion_${_tag} PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 1200 LABELS regression)
endforeach()

# P10.2: kaon production (K+ Lambda H with decay in flight, K+ Sigma- D; K1, rad, seed 4242, 20000 events) against the
# baselines of commit 1275770, in both modes (SIMC_BASELINE_DIR/1275770_[corrected_]kaon_<tag>_n20k_s4242/). Inputs: the
# P10.F1 L2 inputs in reference/runs.
set(SIMC_KAON_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose kaon outputs are the identity baseline")
foreach(_tag h_k1_rad_dec1 d_k1_sigm_rad_dec1)
    set(_inp ${SIMC_REFERENCE_DIR}/runs/kaon_${_tag}_A/kaon_${_tag}_A.inp)
    add_test(NAME regression_identity_kaon_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_kaon_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_KAON_BASELINE_COMMIT}_kaon_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/identity_kaon_${_tag} --fortran-compat --ngen 20000)
    add_test(NAME regression_corrected_kaon_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_kaon_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_KAON_BASELINE_COMMIT}_corrected_kaon_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/corrected_kaon_${_tag} --ngen 20000)
    set_tests_properties(regression_identity_kaon_${_tag} regression_corrected_kaon_${_tag} PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 1200 LABELS regression)
endforeach()

# P11.1: Delta region H(e,e'p)pi (D1: SHMS e-, HMS p; rad, seed 4242, 20000 events) against the baselines of commit 074cced
# (glibc 2.43 / gcc 16.2.1 / ROOT 6.40.04), in both modes (SIMC_BASELINE_DIR/51e86c2_glibc243_[corrected_]delta_h_d1_rad_n20k_s4242/ since P14.2: the 074cced compat
# baseline differs only in the sign of 9 sigcm NaNs, now the Fortran sign (sig_blok operand order, commit 51e86c2);
# validation/scripts/regen_baselines.sh --only delta). Input: the P11.F1 L2 input in reference/runs.
set(SIMC_DELTA_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose Delta outputs are the identity baseline")
foreach(_tag h_d1_rad)
    set(_inp ${SIMC_REFERENCE_DIR}/runs/delta_${_tag}_A/delta_${_tag}_A.inp)
    add_test(NAME regression_identity_delta_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_delta_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_DELTA_BASELINE_COMMIT}_delta_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/identity_delta_${_tag} --fortran-compat --ngen 20000)
    add_test(NAME regression_corrected_delta_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_delta_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_DELTA_BASELINE_COMMIT}_corrected_delta_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/corrected_delta_${_tag} --ngen 20000)
    set_tests_properties(regression_identity_delta_${_tag} regression_corrected_delta_${_tag} PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 1200 LABELS regression)
endforeach()

# P12.2: SIDIS H(e,e'pi+)X (K1: HMS e-, SHMS pi+; rad, seed 4242, 20000 events) against the baselines of commit 31b4232
# (glibc 2.43 / gcc 16.2.1 / ROOT 6.40.04), in both modes (SIMC_BASELINE_DIR/31b4232_glibc243_[corrected_]sidis_h_k1_pip_rad_n20k_s4242/,
# validation/scripts/regen_baselines.sh --only sidis). Input: the P12.F1 L2 input in reference/runs.
set(SIMC_SIDIS_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose SIDIS outputs are the identity baseline")
foreach(_tag h_k1_pip_rad)
    set(_inp ${SIMC_REFERENCE_DIR}/runs/sidis_${_tag}_A/sidis_${_tag}_A.inp)
    add_test(NAME regression_identity_sidis_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_sidis_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_SIDIS_BASELINE_COMMIT}_sidis_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/identity_sidis_${_tag} --fortran-compat --ngen 20000)
    add_test(NAME regression_corrected_sidis_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_sidis_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_SIDIS_BASELINE_COMMIT}_corrected_sidis_${_tag}_n20k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/corrected_sidis_${_tag} --ngen 20000)
    set_tests_properties(regression_identity_sidis_${_tag} regression_corrected_sidis_${_tag} PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 1200 LABELS regression)
endforeach()

# P13.2: exclusive rho H(e,e'pi+)[rho] (K1: HMS e-, SHMS pi; rad, seed 4242, 1000 events: ~2.9e6 tries) against the
# baselines of commit e6feaf0 (glibc 2.43 / gcc 16.2.1 / ROOT 6.40.04), in both modes
# (SIMC_BASELINE_DIR/e6feaf0_glibc243_[corrected_]rho_h_k1_rad_n1k_s4242/, validation/scripts/regen_baselines.sh --only rho).
set(SIMC_RHO_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose rho outputs are the identity baseline")
foreach(_tag h_k1_rad)
    set(_inp ${SIMC_REFERENCE_DIR}/runs/rho_${_tag}_A/rho_${_tag}_A.inp)
    add_test(NAME regression_identity_rho_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_rho_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_RHO_BASELINE_COMMIT}_rho_${_tag}_n1k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/identity_rho_${_tag} --fortran-compat --ngen 1000)
    add_test(NAME regression_corrected_rho_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_rho_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_RHO_BASELINE_COMMIT}_corrected_rho_${_tag}_n1k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/corrected_rho_${_tag} --ngen 1000)
    set_tests_properties(regression_identity_rho_${_tag} regression_corrected_rho_${_tag} PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 1200 LABELS regression)
endforeach()

# P14.3: exclusive H(e,e'pi0)p with the calorimeter (HMS e-, NPS-like arm 8, pizero_ngamma 2; seed 4242, 2000 events) against
# the baselines of commit 17310b4 (glibc 2.43 / gcc 16.2.1 / ROOT 6.40.04), in both modes
# (SIMC_BASELINE_DIR/17310b4_glibc243_[corrected_]calo_excl_h_ng2_n2k_s4242/, validation/scripts/regen_baselines.sh --only calo).
set(SIMC_CALO_BASELINE_COMMIT "55cf7a0_glibc243" CACHE STRING "Commit whose calorimeter outputs are the identity baseline")
foreach(_tag excl_h_ng2)
    set(_inp ${SIMC_REFERENCE_DIR}/runs/p14f2_${_tag}_A/p14f2_${_tag}_A.inp)
    add_test(NAME regression_identity_calo_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_calo_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_CALO_BASELINE_COMMIT}_calo_${_tag}_n2k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/identity_calo_${_tag} --fortran-compat --ngen 2000)
    add_test(NAME regression_corrected_calo_${_tag}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/output_identity.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                     ${_inp} 4242 dir_calo_${_tag}_s4242
                     ${SIMC_BASELINE_DIR}/${SIMC_CALO_BASELINE_COMMIT}_corrected_calo_${_tag}_n2k_s4242
                     ${CMAKE_CURRENT_BINARY_DIR}/corrected_calo_${_tag} --ngen 2000)
    set_tests_properties(regression_identity_calo_${_tag} regression_corrected_calo_${_tag} PROPERTIES
                         SKIP_RETURN_CODE 77 TIMEOUT 1200 LABELS regression)
endforeach()

# P15.2 (D-064/D-065): label "regression_shipped" — full-run identity of every runnable shipped input of simc_gfortran/infiles
# (and the P15.F1 spect_mode / pi0-in-HMS/SHMS runs) against the Fortran archives in reference/runs: simc++ with the RANLUX
# port (seeded from the .inp, as the Fortran) and --fortran-compat; .hist except the time lines, .geni, .gen text-identical,
# h10 tree bitwise (validation/scripts/p15_shipped_regression.sh, one folder per test). Every run here takes < 10 min in
# simc++ (longest p15f1_pi0_shms_A ~2 min), so all are full runs, no short captures. Folders that cover a shipped input:
#   shipped_<name>[_pct]     the P15.F1 sweep (as shipped, or the convert_inputfile.sh renamed copy of a legacy input)
#   nps_* (5), hms_shms_validation (test_hms_shms_validation), eep_h_test_pctnames (test_eep_h renamed).
# Not repeated here to keep check_all short (the reaction families are covered by the shipped_*_pct SOS runs above); their
# full-run identity is on disk (simc_runs/20261114_142803_P15.2_shipped_identity_other): p14f1_{sidis,delta,kaon,pion}_sos_A
# (mc_neg_h_z85_zbin, test_delta_h, test_kaon_h, test_pion_h renamed, 1e5 events, 71-498 s). rhotest (p14f1_rho_sos_A) and
# p14f1_eep_sos_A are the slow l1_simcrun_ranlux_* tests (kaon.cmake, hrs.cmake).
# P15.2b: shipped_poltar_exclusivepion_example and shipped_inputfile_protoype_pct (using_tgt_field, trg_track.f) are full
# 50000-event runs (~615 s each in simc++: the only ones above 10 min; the 1000-event capture p15f2_poltar_cap of the same
# input and the 5 P15.F2 probes p15f2_poltar_* are the short ones). Labels "regression_shipped;slow"; skipped (77) without
# reference.
set(SIMC_SHIPPED_FOLDERS
    p15f1_eep_sm1_A p15f1_eep_smm1_A p15f1_eep_smm2_A p15f1_phsp_sm1_A p15f1_pion_sm1_A
    p15f1_pi0_shms_A p15f1_pi0_hms_A p15f1_pi0_shms_cap
    shipped_eep_hydrogen_q8 shipped_fpi12_test shipped_mc_neg_h_z37_ptbin shipped_test_eep_fe shipped_test_eep_fe_bh
    shipped_poltar_exclusivepion_example shipped_inputfile_protoype_pct shipped_mc_neg_d_z85_zbin_pct
    shipped_test_delta_pct shipped_test_delta2_pct shipped_test_delta_h2_pct shipped_test_eep_d_pct shipped_test_eep_fe_pct
    shipped_test_eep_fe_bh_pct shipped_test_kaon_c12_pct shipped_test_kaon_d_pct shipped_test_kaon_he3_pct
    shipped_test_kaon_he4_pct shipped_test_pion_d_pct shipped_test_pion_he_pct
    nps_excl_pi0_test nps_excl_pi0_ld2_proton_test nps_excl_pi0_ld2_neutron_test nps_delta_pi0_test nps_sidis_pi0_test
    hms_shms_validation eep_h_test_pctnames
    p15f2_poltar_cap p15f2_poltar_pol0 p15f2_poltar_bphi90 p15f2_poltar_bangle5 p15f2_poltar_pi0 p15f2_poltar_bphi0)
foreach(_f ${SIMC_SHIPPED_FOLDERS})
    add_test(NAME regression_shipped_${_f}
             COMMAND ${CMAKE_SOURCE_DIR}/validation/scripts/p15_shipped_regression.sh $<TARGET_FILE:simcpp>
                     $<TARGET_FILE:simc_root_identity> ${CMAKE_CURRENT_BINARY_DIR}/shipped_${_f} ${_f})
    set_tests_properties(regression_shipped_${_f} PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 1800
                         ENVIRONMENT "REF=${SIMC_REFERENCE_DIR}/runs;JOBS=1" LABELS "regression_shipped;slow")
endforeach()

# P15.4: the parallel event loop (--threads N; docs/PERFORMANCE.md). Reproducibility: two runs with the same (seed, N) are
# identical (H(e,e'p) rad: HMS+SHMS; the polarized-target pion example: target field + pion + SHMS); --threads 1 is the
# sequential loop bit for bit; --threads > 1 is refused with --fortran-compat and with the rng "ranlux".
add_test(NAME threads_reproducibility_eep
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/threads_test.sh repro $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                 ${SIMC_FORTRAN_DIR}/infiles/eep_h_elastic_rad_A.inp ${CMAKE_CURRENT_BINARY_DIR}/threads_repro_eep 3 3001)
add_test(NAME threads_reproducibility_poltar
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/threads_test.sh repro $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                 ${SIMC_REFERENCE_DIR}/runs/shipped_poltar_exclusivepion_example/poltar_exclusivepion_example.inp
                 ${CMAKE_CURRENT_BINARY_DIR}/threads_repro_poltar 4 200)
add_test(NAME threads_one_is_sequential
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/threads_test.sh one $<TARGET_FILE:simcpp> $<TARGET_FILE:simc_root_identity>
                 ${SIMC_FORTRAN_DIR}/infiles/eep_h_elastic_rad_A.inp ${CMAKE_CURRENT_BINARY_DIR}/threads_one 2000)
add_test(NAME threads_refusal
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/threads_test.sh refuse $<TARGET_FILE:simcpp>
                 ${SIMC_FORTRAN_DIR}/infiles/eep_h_elastic_rad_A.inp ${CMAKE_CURRENT_BINARY_DIR}/threads_refusal)
set_tests_properties(threads_reproducibility_eep threads_reproducibility_poltar threads_one_is_sequential threads_refusal
                     PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900 LABELS "regression;threads")
