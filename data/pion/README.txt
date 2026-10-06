Pion cross-section data (task P9.1). Byte-identical copies of the files that the Fortran SIMC opens in its working
directory (simc_gfortran/), md5 sums in MD5SUMS. The CTest l1_pion_data_files checks that the copies still equal the
Fortran files (when the Fortran tree is present) and that their md5 sums are the ones listed here.

  maidpi0p.dat   physics_pion.f:sigmaid (l.577-736)   MAID2007 table, ipi = 1: pi0 p
  maidpi0n.dat   "                                    ipi = 2: pi0 n
  maidpipn.dat   "                                    ipi = 3: pi+ n
  maidpimp.dat   "                                    ipi = 4: pi- p
  MULTIPOL.file  physics_pion.f:sig_multipole (l.197-398), multipole amplitudes (dead code in the Fortran: its call
                 in peepi is commented out; ported for the L1 table only)

Formats (read exactly as the Fortran does):
  maid*.dat      25 (Q2) x 46 (W) x 23 (theta) records, format (f11.6,18f8.4) (physics_pion.f:620-635). The records
                 are 147 characters long, so the 18th field is blank (= 0, Fortran PAD mode). The lines after record
                 26450 (a MAID text trailer) are not read.
  MULTIPOL.file  8 Q2 blocks: format (/11X,f5.1/) (-Q**2 in fm**-2), then 10 x [(11X,f5.0) KLAB, (6(/9X,7f8.3))].
                 The nu_table is overwritten by every Q2 block (the last block's nodes are used); the third block
                 lists 360 MeV instead of 380 MeV (as in the Fortran file). Lines 665-673 are not read.

simc++ looks for these files in the "pion" physics model's JSON key "data_dir" (default: this directory, compiled in
as SIMC_PION_DATA_DIR):
  "pipeline": {"physics": {"type": "pion", "data_dir": "/path/to/dir"}}
maidpi0pborn.dat (also in simc_gfortran) is not read by any Fortran routine and is not copied.
