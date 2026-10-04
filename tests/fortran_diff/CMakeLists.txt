# tests/fortran_diff — task P7.5 (D-031, D-032): the documented differences with the Fortran (docs/FORTRAN_DIFFERENCES.md).
# For every corrected clear Fortran bug (class a) one test shows BOTH behaviours: fortran_compat = true gives the
# Fortran result, the default (corrected) gives the corrected one. Label "fortran_diff" (+ "unit"/"cli").
add_executable(test_fortran_differences test_fortran_differences.cpp)
target_link_libraries(test_fortran_differences PRIVATE ${SIMC_DRIVER_LIBS})
target_compile_definitions(test_fortran_differences PRIVATE SIMC_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# One CTest per item id (the ids of docs/FORTRAN_DIFFERENCES.md); the test program takes the id as argument.
set(SIMC_FORTRAN_DIFF_ITEMS
    D032_spectrometer_id
)
foreach(_item ${SIMC_FORTRAN_DIFF_ITEMS})
    add_test(NAME fortran_diff_${_item}
             COMMAND test_fortran_differences ${_item} ${SIMC_REFERENCE_DIR} ${SIMC_FORTRAN_DIR})
    set_tests_properties(fortran_diff_${_item} PROPERTIES SKIP_RETURN_CODE 77 LABELS "fortran_diff;unit")
endforeach()

add_test(NAME cli_arm_code_6_rejected
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cli_arm_code_rejected.sh $<TARGET_FILE:simcpp> $<TARGET_FILE:inp2json>
                 ${CMAKE_CURRENT_BINARY_DIR}/cli_arm_code_rejected)
set_tests_properties(cli_arm_code_6_rejected PROPERTIES LABELS "fortran_diff;cli")
