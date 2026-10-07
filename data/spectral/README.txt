Spectral-function data for A(e,e'p) (task P8.1). Byte-identical copies of the files that the Fortran SIMC opens in its
working directory (simc_gfortran/), md5 sums in MD5SUMS. The CTest l1_spectral_files checks that the copies still equal
the Fortran files.

  h2.theory      init.f:theory_init  D(e,e'p): 1-D momentum distribution (doing_deuterium)          A = 2
  c12.theory     init.f:theory_init  A(e,e'p) shell model: 2 shells (1s, 1p), E_Fermi 15.96 MeV     A = 12
  fe56.theory    init.f:theory_init  5 shells, E_Fermi 10.18 MeV                                     A = 56
  au197.theory   init.f:theory_init  10 shells, E_Fermi 5.78 MeV                                     A = 197
  benharsf_3mod.dat  sf_lookup.f:sf_lookup_init  Benhar SF (use_benhar_sf), He-3 "mod" version   A = 3
  benharsf_4.dat     "                                                               A = 4
  benharsf_12.dat    "                                                               A = 12 (and the Fortran default)
  benharsf_56.dat    "                                                               A = 56
  benharsf_197.dat   "                                                               A = 197

The Fortran chooses the file from nint(targ%A) (init.f:839-851, dbase.f:598-616) and silently falls back to the C-12
file for every other A (He-3, He-4, Ca-40, ...; D-033 F1). SIMC++ takes the same default, prints a warning when it
falls back, and lets the user name the file explicitly:

  "pipeline": {"reaction": {"type": "nucleus_eep",
                            "spectral_function": {"type": "theory", "file": "c12.theory"}}}
  "pipeline": {"reaction": {"type": "nucleus_eep",
                            "spectral_function": {"type": "benhar", "file": "benharsf_4.dat"}}}

simc++ looks for these files in --data-dir DIR (default: this directory, compiled in as SIMC_SPECTRAL_DATA_DIR).
File formats: theory files = init.f:855-886 (nrhoPm absorption E_Fermi / per shell nprot Em Emsig bs_norm /
(Pm, rho) points, a shell ends when Pm decreases); Benhar files = sf_lookup.f:26-30 (numPm numEm / numPm*numEm lines
Pm Em SFp SFn dPm dEm).

Momentum distributions of the struck nucleon for quasifree pion/kaon/SIDIS production on D and A > 2 (task P9.2;
dbase.f:563-587: pval [MeV/c], mprob read list-directed, mprob normalised to mprob(nump)):
  deut.dat       doing_deutpi / doing_deutkaon / doing_deutsemi                        A = 2
  he3.dat        doing_hepi / doing_hekaon                                             A = 3
  he4.dat        "                                                                     A = 4
  c12.dat        "                                                                     A = 12 (and the Fortran default)
For A > 2 quasifree production the Benhar table of the same A (benharsf_*.dat) is also read: VERTEXedge%Em%max =
Emval(numEm) (init.f:356) and the missing energy of each event (sf_lookup.f:generate_em).
