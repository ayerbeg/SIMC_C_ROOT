# tests/l1/shms.cmake — task P6.3: ISpectrometer "shms" (mc_shms, mc_shms_coll, mc_shms_hut) against the Fortran
# tables of P6.F1. Included (OPTIONAL) by tests/l1/CMakeLists.txt.
#   l1_mc_shms, l1_mc_shms_smear : every drv_mc.f column; |rel| <= 1e-12 for real*8, exact for the real*4 focal
#                                  plane, stop code / stop_id / RNG consumption exact, counters summed over the table
#   l1_mc_p12_check_shms         : ss* predictions of chk_p12.f (eep_h_elastic_norad_A events)
#   l1_mc_dist_shms_{nosmear,smear} : 1e6 Mt64 rays vs the RANLUX mc_distributions (pulls, chi2, KS, quantiles)
add_executable(test_l1_shms test_l1_shms.cpp)
# simcfw_hms: hms::pionCollAbsorb (hms/pion_coll_absorb.f, called by mc_shms_coll)
target_link_libraries(test_l1_shms PRIVATE simcfw_l1runner simcfw_shms simcfw_hms simcfw_optics simcfw_material)
target_compile_options(test_l1_shms PRIVATE -O2)

set(_fns ${SIMC_REFERENCE_DIR}/functions)
add_test(NAME l1_mc_shms COMMAND test_l1_shms mc ${_fns}/mc_shms 1e-12)
add_test(NAME l1_mc_shms_smear COMMAND test_l1_shms mc ${_fns}/mc_shms_smear 1e-12)
add_test(NAME l1_mc_p12_check_shms COMMAND test_l1_shms p12 ${_fns}/mc_p12_check 1e-12)
add_test(NAME l1_mc_dist_shms_nosmear COMMAND test_l1_shms dist ${_fns}/mc_distributions shms_nosmear 1000000 20260932)
add_test(NAME l1_mc_dist_shms_smear COMMAND test_l1_shms dist ${_fns}/mc_distributions shms_smear 1000000 20260933)
add_test(NAME shms_selftest COMMAND test_l1_shms selftest)
set_tests_properties(l1_mc_shms l1_mc_shms_smear l1_mc_p12_check_shms l1_mc_dist_shms_nosmear l1_mc_dist_shms_smear
    PROPERTIES SKIP_RETURN_CODE 77)
set_tests_properties(l1_mc_dist_shms_nosmear l1_mc_dist_shms_smear PROPERTIES TIMEOUT 900 LABELS "slow")
