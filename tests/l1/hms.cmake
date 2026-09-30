# tests/l1/hms.cmake — task P6.2: ISpectrometer "hms" (mc_hms, mc_hms_coll, mc_hms_hut, lfit) against the Fortran
# tables of P6.F1. Included (OPTIONAL) by tests/l1/CMakeLists.txt.
#   l1_mc_hms, l1_mc_hms_smear : every drv_mc.f column; |rel| <= 1e-12 for real*8, exact for the real*4 focal plane,
#                                stop code / RNG consumption exact, counters summed over the table
#   l1_mc_p12_check_hms        : hs* predictions of chk_p12.f (eep_h_elastic_norad_A events)
#   l1_mc_dist_hms_{nosmear,smear} : 1e6 Mt64 rays vs the RANLUX mc_distributions (pulls, chi2, KS, quantiles)
add_executable(test_l1_hms test_l1_hms.cpp)
target_link_libraries(test_l1_hms PRIVATE simc_l1runner simcfw_hms simcfw_optics simcfw_material)
target_compile_options(test_l1_hms PRIVATE -O2)

set(_fnh ${SIMC_REFERENCE_DIR}/functions)
add_test(NAME l1_mc_hms COMMAND test_l1_hms mc ${_fnh}/mc_hms 1e-12)
add_test(NAME l1_mc_hms_smear COMMAND test_l1_hms mc ${_fnh}/mc_hms_smear 1e-12)
add_test(NAME l1_mc_p12_check_hms COMMAND test_l1_hms p12 ${_fnh}/mc_p12_check 1e-12)
add_test(NAME l1_mc_dist_hms_nosmear COMMAND test_l1_hms dist ${_fnh}/mc_distributions hms_nosmear 1000000 20260930)
add_test(NAME l1_mc_dist_hms_smear COMMAND test_l1_hms dist ${_fnh}/mc_distributions hms_smear 1000000 20260931)
add_test(NAME hms_selftest COMMAND test_l1_hms selftest)
set_tests_properties(l1_mc_hms l1_mc_hms_smear l1_mc_p12_check_hms l1_mc_dist_hms_nosmear l1_mc_dist_hms_smear
    PROPERTIES SKIP_RETURN_CODE 77)
set_tests_properties(l1_mc_dist_hms_nosmear l1_mc_dist_hms_smear PROPERTIES TIMEOUT 900 LABELS "slow")
