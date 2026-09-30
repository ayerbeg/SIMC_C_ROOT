# tests/l1/optics.cmake — task P6.1: COSY optics engine (transp, chains, recon) against the Fortran tables.
# Included (OPTIONAL) by tests/l1/CMakeLists.txt. Tolerance |rel| <= 1e-12 (all real*8 paths).
# The COSY files are the repository copies in data/matrices/simc_fortran/{hms,shms}/ (byte-identical to the Fortran ones,
# checked by l1_optics_cosy_files).
add_executable(test_l1_optics test_l1_optics.cpp)
target_link_libraries(test_l1_optics PRIVATE simc_l1runner simcfw_optics simcfw_material)

set(_fno ${SIMC_REFERENCE_DIR}/functions)
set(_mat ${CMAKE_SOURCE_DIR}/data/matrices/simc_fortran)
foreach(_t transp_hms transp_shms)
    add_test(NAME l1_${_t} COMMAND test_l1_optics transp ${_mat} ${_fno}/${_t} 1e-12)
endforeach()
foreach(_t transp_hms_chain transp_shms_chain)
    add_test(NAME l1_${_t} COMMAND test_l1_optics chain ${_mat} ${_fno}/${_t} 1e-12)
endforeach()
foreach(_t recon_hms recon_shms)
    add_test(NAME l1_${_t} COMMAND test_l1_optics recon ${_mat} ${_fno}/${_t} 1e-12)
endforeach()
add_test(NAME l1_optics_cosy_files COMMAND test_l1_optics files ${_mat} ${SIMC_FORTRAN_DIR})
add_test(NAME optics_selftest COMMAND test_l1_optics selftest ${_mat})
set_tests_properties(l1_transp_hms l1_transp_shms l1_transp_hms_chain l1_transp_shms_chain l1_recon_hms l1_recon_shms
    l1_optics_cosy_files PROPERTIES SKIP_RETURN_CODE 77)
