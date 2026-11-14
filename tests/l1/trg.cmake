# tests/l1/trg.cmake — task P15.2b: polarized-target field tracking (optics::TargetField, trg_track.f) against the
# P15.F2 tables reference/functions/trg_track/<table>/ (bitwise, tolerance 0). The map is the repository copy
# data/field/trg_field_map.dat (l1_trg_field_map_file checks it equals the Fortran file); track_to_tgt uses the COSY
# recon files of data/matrices/simc_fortran.
add_executable(test_l1_trg test_l1_trg.cpp)
target_link_libraries(test_l1_trg PRIVATE simcfw_l1runner simcfw_optics simcfw_material)
set(_trg ${SIMC_REFERENCE_DIR}/functions/trg_track)
set(_trgmap ${CMAKE_SOURCE_DIR}/data/field)
set(_trgmat ${CMAKE_SOURCE_DIR}/data/matrices/simc_fortran)
set(_trg_tests)
foreach(_t trginit trginit_map trgfield trgderiv trgrk4 trgtracktoplane track_from_tgt track_to_tgt)
    add_test(NAME l1_trg_${_t} COMMAND test_l1_trg ${_t} ${_trg}/${_t} ${_trgmap} ${_trgmat})
    list(APPEND _trg_tests l1_trg_${_t})
endforeach()
add_test(NAME l1_trg_field_map_file COMMAND test_l1_trg files ${_trgmap} ${SIMC_FORTRAN_DIR})
list(APPEND _trg_tests l1_trg_field_map_file)
set_tests_properties(${_trg_tests} PROPERTIES SKIP_RETURN_CODE 77 LABELS l1 TIMEOUT 600)
# the capture replay: all 704778 trgTrackToPlane, 203423 track_from_tgt and the 9174 tracking track_to_tgt calls of the
# 1000-event capture of the shipped poltar example (bitwise); needs xz
add_test(NAME l1_trg_replay_capture COMMAND test_l1_trg replay ${SIMC_REFERENCE_DIR}/runs/p15f2_poltar_cap/capture
         ${_trgmap} ${_trgmat})
set_tests_properties(l1_trg_replay_capture PROPERTIES SKIP_RETURN_CODE 77 LABELS "l1;replay" TIMEOUT 900)
