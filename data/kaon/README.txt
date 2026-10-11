Kaon cross-section model data (task P10.1). Byte-identical copies of the files that the Fortran SIMC opens in its
working directory (simc_gfortran/) for doing_kaon (dbase.f:624-679), md5 sums in MD5SUMS. The CTest l1_kaon_data_files
checks that the copies still equal the Fortran files and that their md5 sums are the ones listed here.

  weightc.dat        weightc(20,50), list-directed, one value per record   (physics_kaon.f:phspwght, e- HMS + K SOS; dead code)
  weightd.dat        weightd(8,40,30), list-directed                        (phspwght, e- SOS + K HMS; dead code)
  saghai_proton.dat  zrff1..6/ziff1..6(10,11,19): per (s, Q2) a list-directed header "s Q2", then per angle 3 records
                     (6e12.4): a 5-value kinematics line (not used) and the two amplitude lines. 14 s-blocks in the file,
                     the first 10 are read (s = 2.6 ... 4.7, 5.1, 5.4 GeV^2; D-043 F3: eekeek's grid says 5.0, 5.3).
  saghai_sigma0.dat  zsrff1..6/zsiff1..6(20,10,19): per record a header "s Q2", a 5-value kinematics line, ff1..ff3,
                     ff4..ff6. Line 1 is a truncated fragment, so dbase.f reads every record ONE LINE OFF (D-043 F2):
                     fortran_compat reproduces that; the default mode skips the fragment and reads each record from its
                     own lines. The file content also starts with 18 stray records (s = 2.898, Q2 = 0.52/0.75) in the
                     slots of s = 2.851, Q2 = 0 (kept as the file orders them; FORTRAN-DIFF-CANDIDATE), and line 72 has a
                     truncated field '-00E-02' (read as -0, Fortran Ew.d rules).
These tables feed only the Saghai model (ntup%sigcm1, ntuple column "saghai"); the event weight uses the factorized
model (sig_factorized), which needs no data file.

simc++ looks for these files in the "kaon" physics model's JSON key "data_dir" (default: this directory, compiled in
as SIMC_KAON_DATA_DIR):  "pipeline": {"physics": {"type": "kaon", "data_dir": "/path/to/dir"}}
