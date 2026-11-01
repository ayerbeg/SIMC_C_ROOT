Byte-identical copies of the COSY files the Fortran SIMC reads (task P6.1), under the relative paths
that transp_init / mc_hms_recon / mc_shms_recon / mc_sos_recon open (so this directory plays the role of the SIMC
directory for simc::optics::CosyEngine::initForward / initRecon):

  hms/forward_cosy.dat   <- simc_gfortran/hms/forward_cosy.dat                       (12 classes)
  hms/recon_cosy.dat     <- simc_gfortran/hms/recon_cosy.dat                         (461 terms)
  shms/shms_forward.dat  <- simc_gfortran/shms/shms_forward.dat -> shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_hb80_forward.dat (32 classes)
  shms/shms_recon.dat    <- simc_gfortran/shms/shms_recon.dat   -> shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_hb80_recon_fit.dat (461 terms)
  sos/forward_cosy.dat   <- simc_gfortran/sos/forward_cosy.dat  (10 classes; transp.f:302, task P14.1)
  sos/recon_cosy.dat     <- simc_gfortran/sos/recon_cosy.dat    (mc_sos_recon.f:73, task P14.1)
  hrs/hrs_forward_cosy.dat <- simc_gfortran/hrsl/hrs_forward_cosy.dat = hrsr/hrs_forward_cosy.dat (byte-identical, 12 classes;
                              transp.f:304/306, task P14.2: ONE copy for both HRS arms)
  hrs/hrs_recon_cosy.dat   <- simc_gfortran/hrsl/hrs_recon_cosy.dat = hrsr/hrs_recon_cosy.dat (mc_hrs*_recon.f:73)

md5 sums: MD5SUMS (identical to those listed in reference/functions/transp_hms/README.txt).
The CTest l1_optics_cosy_files checks that these copies still equal the Fortran files.

Why not ../hms/recon_cosy.dat: the older copy there (used by the legacy code) was hand-edited — blanks were
inserted between fixed-width fields such as "0.273541506E-01-0.671779120E-02" — so it differs from the
Fortran file (md5 aac2e08a... vs 639364a8...). The Fortran reads these files with fixed-column formats
(format(1x,5g14.7,1x,6i1) forward, format(1x,4g16.9,1x,5i1) recon); CosyEngine does the same and needs
the original layout.
