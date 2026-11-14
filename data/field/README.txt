Polarized-target field map (task P15.2b). Byte-identical copy of simc_gfortran/trg_field_map.dat, the map that
trg_track.f:trgInit (l.296-336) reads when the input has using_tgt_field = 1 and tgt_field_file = 'trg_field_map.dat'
(the shipped poltar_exclusivepion_example.inp). md5 in MD5SUMS; the CTest l1_trg_field_map_file checks that the copy
still equals the Fortran file.

Format (read list-directed, 2601 lines): r outer loop, z inner loop, 51 x 51 nodes on the hard-coded grid
z = r = 0, 2, ..., 100 cm; 7 numbers per line, of which only column 3 (Bz [T]) and column 4 (Br [T]) are used.

simc++ resolves tgt_field_file as the Fortran working directory would: an absolute path as is; otherwise the file in
--data-dir DIR (if given and present), then the current directory, then this directory (compiled in as
SIMC_FIELD_DATA_DIR). tgt_field_file = '0' is the zero field and ' ' the uniform 5 T test field (trgInit).
