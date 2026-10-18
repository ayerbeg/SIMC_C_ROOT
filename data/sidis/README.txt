SIDIS model data (task P12.1). Byte-identical copies of the files that the Fortran SIMC opens (relative to its working
directory simc_gfortran/) for doing_semi; md5 sums in MD5SUMS. The CTest l1_sidis_data_files checks that the copies still
equal the Fortran files and that their md5 sums are the ones listed here.

  cteq5m.tbl   simc_gfortran/cteq5/cteq5m.tbl   CTEQ5M parton distributions (Ctq5Pdf.f:SetCtq5(1) + ReadTbl): two
               header lines, "Dr Fl Al Amass(1:6)", "NX NT NfMx" = 90 13 5, QINI QMAX QL(0:13), XMIN XV(0:90), then
               Npts = 91*14*8 = 10192 list-directed values UPD(1:Npts) (x f(x,Q) on the (x, Q) grid for the flavours
               -NfMx..2 of Ctq5Pdf's offset JFL = Ip + NfMx). Read by peepiX (semi_physics.f:300) and semi_dilution
               (dead code) with iset = 1; the other cteq5*.tbl sets of simc_gfortran/cteq5/ are never read by SIMC.
  KANLO.GRID   simc_gfortran/fdss/KANLO.GRID    fDSS NLO kaon fragmentation functions (fdss.f:fDSS, IH = 2, IO = 1),
               34*24 = 816 records of FORMAT (9(1PE10.3)) (the last record has no newline): PARTON(1:9, IQ, IX) =
               the total u, d, s, c, b, g and the valence u, d, s combinations (z times the FF, fdss.f:113-121, 141-158)
               at the first 34 z and the 24 Q2 grid nodes of fdss.f:57-63. Read by peepiX for doing_semika only (K+ and K-); the pion, proton and
               hadron grids of simc_gfortran/fdss/ are never read by SIMC (the pion FFs are P. Bosted's fit).

The f1f2tables/*cdbonn* deuteron smearing tables that peepiX's first call of SQESUB reads have no effect on the SIDIS
cross section (F1F2IN21 is called with A = 1) and are not needed.

simc++ looks for these files in the "sidis" physics model's JSON key "data_dir" (default: this directory, compiled in
as SIMC_SIDIS_DATA_DIR):  "pipeline": {"physics": {"type": "sidis", "data_dir": "/path/to/dir"}}
