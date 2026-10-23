# SIMC++ output formats (P1.4)

The C++ writes outputs that can be compared directly with the Fortran SIMC outputs:

| Fortran output | C++ sink (registry name) | C++ file |
|---|---|---|
| `worksim/<run>.bin` → `make_root_tree` → `<run>.root`, tree `h10` | `RootTreeSink` (`root_tree`) | `<output_dir>/<run>.root`, tree `h10` |
| `outfiles/<run>.hist` (simc.f:report) + `outfiles/<run>.geni` (stop counters) | `RunSummarySink` (`run_summary`) | `<run>.summary.json` (+ `<run>.summary.txt`) |
| `outfiles/<run>.gen` (diagnostic histograms) | `GenHistSink` (`gen_hist`) | `<run>.gen.root` (TH1D) |

`tools/hist2json` converts a Fortran `.hist` + `.geni` into the same `summary.json` layout. The comparison tool can
therefore read both codes in one format.

All sinks implement `IOutputSink` (include/simc/stages/IOutputSink.h). The driver calls `beginRun(setup)` once,
`fill(main, vertex, orig, recon, success)` once per **tried** event, and `endRun(summary)` once.
Each sink accepts an optional JSON config `{"file": "..."}`.

## 1. Tree `h10` (RootTreeSink)

- Source of truth: NtupleInit.f (tags and order), results_write.f:results_ntu_write (values),
  util/root_tree/make_root_tree.f + froot.c (tree `h10`, one `Float_t` leaf per tag, branch title = name).
- Only successful events are written (results_write.f:34).
- Each value is computed in double precision and stored as `Float_t`, like make_root_tree's `real*4 ntup_out`.
- Implemented layouts: the (e,e'p) family (`doing_hyd_elast`, `doing_deuterium`, `doing_heavy`), 46 columns. This includes
  `doing_phsp`, which keeps `doing_hyd_elast`, so the acceptance runs also have these 46 columns. Exclusive pion production
  (`doing_pion`, P9.2): 53 columns, 65 with `doing_pizero` (§1b). Kaon production (`doing_kaon`, P10.2): 55 columns (§1c).
  Delta region (`doing_delta`, P11.1): the 53 pion columns (§1d). SIDIS (`doing_semi`, P12.2): 56 columns (§1e). The rho
  layout (NtupleInit.f l.191-295) will come with that reaction; `beginRun` throws for it (and for `using_tgt_field` and π0 SIDIS). The layout is chosen from the run
  flags (`NtupleColumnsFor(setup)`, `FillNtuple`).
- **Arm assignment** (results_write.f:64 and 89): `hs*` = the spectrometer on the beam **right** (HMS, HRSR) and `ss*` = the
  one on the **left** (SOS, HRSL, SHMS).
  - `electron_arm` = 1, 3, 7 (beam-right spectrometers, `core/SpectrometerId.h`): `hs*` holds the electron, `ss*` the hadron.
  - `electron_arm` = 2, 4, 5, 8: `hs*` holds the hadron and `ss*` the electron. (The Fortran also lists 6, the retired SHMS
    LSA tune; SIMC++ rejects 6 at input, D-032.)
  - Any other value: columns 1-24 keep the previous event's values. This matches the static Fortran array, which the Fortran
    reports as "results_write not yet set up".
  - In the table, R and L are the right-side and left-side particle: `e` or `p` depending on the case above.
- Units are those of the Fortran ntuple. The `!mr` comments in results_write.f are wrong: the stored values are slopes
  (rad), as the reference tree shows (hsyptar is about ±0.03).

| # | branch | C++ source | units | Fortran |
|---|---|---|---|---|
| 1 | `hsdelta` | `recon.R.delta` | % | results_write.f:65/91 |
| 2 | `hsyptar` | `recon.R.yptar` | slope (rad) | results_write.f:66/92 |
| 3 | `hsxptar` | `recon.R.xptar` | slope (rad) | results_write.f:67/93 |
| 4 | `hsytar` | `recon.R.z` | cm | results_write.f:68/94 |
| 5 | `hsxfp` | `main.FP.R.x` | cm | results_write.f:69/95 |
| 6 | `hsxpfp` | `main.FP.R.dx` | slope (rad) | results_write.f:70/96 |
| 7 | `hsyfp` | `main.FP.R.y` | cm | results_write.f:71/97 |
| 8 | `hsypfp` | `main.FP.R.dy` | slope (rad) | results_write.f:72/98 |
| 9 | `hsdeltai` | `orig.R.delta` | % | results_write.f:73/102 |
| 10 | `hsyptari` | `orig.R.yptar` | slope (rad) | results_write.f:74/103 |
| 11 | `hsxptari` | `orig.R.xptar` | slope (rad) | results_write.f:75/104 |
| 12 | `hsytari` | `main.target.z * sin(theta_R)` | cm | results_write.f:76/105 |
| 13 | `ssdelta` | `recon.L.delta` | % | results_write.f:77/106 |
| 14 | `ssyptar` | `recon.L.yptar` | slope (rad) | results_write.f:78/107 |
| 15 | `ssxptar` | `recon.L.xptar` | slope (rad) | results_write.f:79/108 |
| 16 | `ssytar` | `recon.L.z` | cm | results_write.f:80/109 |
| 17 | `ssxfp` | `main.FP.L.x` | cm | results_write.f:81/110 |
| 18 | `ssxpfp` | `main.FP.L.dx` | slope (rad) | results_write.f:82/111 |
| 19 | `ssyfp` | `main.FP.L.y` | cm | results_write.f:83/112 |
| 20 | `ssypfp` | `main.FP.L.dy` | slope (rad) | results_write.f:84/113 |
| 21 | `ssdeltai` | `orig.L.delta` | % | results_write.f:85/114 |
| 22 | `ssyptari` | `orig.L.yptar` | slope (rad) | results_write.f:86/115 |
| 23 | `ssxptari` | `orig.L.xptar` | slope (rad) | results_write.f:87/116 |
| 24 | `ssytari` | `-main.target.z * sin(theta_L)` | cm | results_write.f:88/117 |
| 25 | `q` | `recon.q/1000` | GeV/c | results_write.f:121 |
| 26 | `nu` | `recon.nu/1000` | GeV | results_write.f:122 |
| 27 | `Q2` | `recon.Q2/1e6` | (GeV/c)^2 | results_write.f:123 |
| 28 | `W` | `recon.W/1000` | GeV | results_write.f:124 |
| 29 | `epsilon` | `recon.epsilon` |  | results_write.f:125 |
| 30 | `Em` | `recon.Em/1000` | GeV | results_write.f:126 |
| 31 | `Pm` | `recon.Pm/1000` | GeV/c | results_write.f:127 |
| 32 | `thetapq` | `recon.theta_pq` | rad | results_write.f:128 |
| 33 | `phipq` | `recon.phi_pq` | rad | results_write.f:129 |
| 34 | `corrsing` | `(recon.e.P - Mp*Ebeam/(2*Ebeam*sin^2(recon.e.theta/2)+Mp))/1000` | GeV | results_write.f:55-56,247 |
| 35 | `Pmx` | `-recon.Pmx/1000` | GeV/c | results_write.f:61,248 |
| 36 | `Pmy` | `(recon.Pmz*uq.y - recon.Pmy*uq.z)/sqrt(uq.y^2+uq.z^2)/1000` | GeV/c | results_write.f:59-60,249 |
| 37 | `Pmz` | `-(recon.Pmy*uq.y + recon.Pmz*uq.z)/sqrt(uq.y^2+uq.z^2)/1000` | GeV/c | results_write.f:57-58,250 |
| 38 | `PmPar` | `recon.PmPar/1000` | GeV/c | results_write.f:251 |
| 39 | `PmPer` | `recon.PmPer/1000` | GeV/c | results_write.f:252 |
| 40 | `PmOop` | `recon.PmOop/1000` | GeV/c | results_write.f:253 |
| 41 | `fry` | `-main.target.rastery` | cm | results_write.f:254 |
| 42 | `radphot` | `main.ntup.radphot/1000` | GeV | results_write.f:255 |
| 43 | `sigcc` | `main.sigcc` | ub/sr (H) or ub*MeV^2/sr^2 | results_write.f:256 |
| 44 | `Weight` | `main.weight` |  | results_write.f:257 |
| 45 | `theta_e` | `recon.e.theta` | rad | results_write.f:258 |
| 46 | `theta_p` | `recon.p.theta` | rad | results_write.f:259 |

Signs:
- `hsytari = +z·sinθ(right arm)` and `ssytari = −z·sinθ(left arm)`.
- `fry = −main.target.rastery`.
- `Pmx/Pmy/Pmz` are the "Heepcheck" components (results_write.f:57-61).
- `corrsing` uses the nominal `Ebeam` and Mp = 938.27231.

### 1b. Pion layout (doing_pion; NtupleInit.f:101-190, results_write.f:131-186)

Columns 1-33 are those of the table above (same arm assignment). Then:

| # | branch | C++ source | units | Fortran |
|---|---|---|---|---|
| 34 | `missmass` | `main.ntup.mm/1000` (complete_recon_ev: sign(mm²)·sqrt\|Em²−Pm²\| with Em = ν + Mtar_struck − E_π) | GeV | results_write.f:134 |
| 35 | `mmnuc` | `main.ntup.mmA/1000` (missing mass on the nucleus, M instead of Mtar_struck) | GeV | results_write.f:135 |
| 36 | `phad` | `recon.p.P/1000` | GeV/c | results_write.f:136 |
| 37 | `t` | `main.ntup.t/1e6` (recon; SIMC sign: −t_Mandelstam) | GeV² | results_write.f:137 |
| 38 | `pmpar` | `recon.PmPar/1000` | GeV/c | results_write.f:138 |
| 39 | `pmper` | `recon.PmPer/1000` | GeV/c | results_write.f:139 |
| 40 | `pmoop` | `recon.PmOop/1000` | GeV/c | results_write.f:140 |
| 41 | `fry` | `-main.target.rastery` | cm | results_write.f:141 |
| 42 | `radphot` | `main.ntup.radphot/1000` | GeV | results_write.f:142 |
| 43 | `pfermi` | `pfer/1000·\|d\|/d`, d = p̂fermi·q̂(vertex) (d = 0 → 1e-20) | GeV/c | results_write.f:143-145 |
| 44 | `siglab` | `main.sigcc` (d⁵σ incl. Delta factor and Coulomb factor) | µb/MeV/sr² | results_write.f:146 |
| 45 | `sigcm` | `main.ntup.sigcm` (peepi) | µb/MeV²/rad | results_write.f:147 |
| 46 | `Weight` | `main.weight` |  | results_write.f:148 |
| 47 | `decdist` | `main.ntup.decdist` (0 without a decay) | cm | results_write.f:149 |
| 48 | `Mhadron` | `sqrt(main.ntup.Mh2_final)` (the muon mass after a π decay) | MeV | results_write.f:150 |
| 49 | `pdotqhat` | `pfer/1000·d` | GeV/c | results_write.f:151 |
| 50 | `Q2i` | `vertex.Q2/1e6` | GeV² | results_write.f:152 |
| 51 | `Wi` | `main.w/1000` (vertex, Mtar_struck) | GeV | results_write.f:153 |
| 52 | `ti` | `main.t/1e6` (vertex) | GeV² | results_write.f:154 |
| 53 | `phipqi` | `main.phi_pq` (vertex) | rad | results_write.f:155 |
| 54-65 | `xcal_gamma1`, `ycal_gamma1`, `Egamma1`, `Pgamma1x/y/z`, `xcal_gamma2`, `ycal_gamma2`, `Egamma2`, `Pgamma2x/y/z` | doing_pizero only: the calorimeter positions (written only by mc_calo, not ported: 0) and the lab photons of pizero_decay (MeV) | cm, MeV | results_write.f:172-185 |

`pfer`, `pferx..z` are COMMON /pfermi/ (`MainEvent::pfermi`) of the event. The validation captures prove the 53 columns
bit for bit (`l1_simcrun_replay_pion_*`).

### 1c. Kaon layout (doing_kaon; NtupleInit.f:101-166, results_write.f:131-166)

Columns 1-53 are those of the pion layout (§1b), with the kaon meanings below, then two kaon columns:

| # | branch | C++ source | units | Fortran |
|---|---|---|---|---|
| 44 | `siglab` | `main.sigcc` (peeK, incl. the Coulomb factor; without the survival probability) | µb/MeV/sr² | results_write.f:146 |
| 45 | `sigcm` | `main.ntup.sigcm` (= sig_factorized) | model units | results_write.f:147 |
| 46 | `Weight` | `main.weight` (× survival probability when doing_decay = 0) |  | results_write.f:148 |
| 47 | `decdist` | `main.ntup.decdist`: the decay distance [cm] with decay ON; the **survival probability** with decay OFF (physics_kaon.f:164) | cm or 1 | results_write.f:149 |
| 48 | `Mhadron` | `sqrt(main.ntup.Mh2_final)` (μ or π mass after a K decay; a first-half K→ππ decay in transp keeps Mk in fortran_compat, D-006) | MeV | results_write.f:150 |
| 54 | `saghai` | `main.ntup.sigcm1` (Saghai model; 0 for every event in fortran_compat, D-043 F1) | model units | results_write.f:164 |
| 55 | `factor` | `main.ntup.sigcm2` (factorized model, = `sigcm`) | model units | results_write.f:165 |

With `using_tgt_field` the two kaon columns would move to 62/63 (not ported). The kaon captures prove the 55 columns bit for
bit (`l1_simcrun_replay_kaon_*`).

### 1d. Delta layout (doing_delta; NtupleInit.f:101-158, results_write.f:131-155)

The 53 columns of the pion layout (§1b), with these meanings for H/D/A(e,e'p)π (the proton is detected, the pion missing):

| # | branch | content for doing_delta |
|---|---|---|
| 34 | `missmass` | the missing (pion) mass, recon |
| 36 | `phad` | the proton momentum, recon |
| 37, 52 | `t`, `ti` | t of the proton vertex (Mh² = Mp²), recon / vertex |
| 44 | `siglab` | `main.sigcc` = peedelta: **phase space × photon flux × jacobian, NO cross section** (D-048 F4) |
| 45 | `sigcm` | `main.ntup.sigcm` = sig_blok (the Brauel π fit at the proton CM angles; not in the weight). fortran_compat: the Fortran value read from the low 8 bytes of binary128 arguments, meaningless, NaN/Inf in ~0.4 % of the events (D-048 F3); default: the fit with the intended arguments |
| 46 | `Weight` | `main.weight` (no cross section, F4) |
| 47, 48 | `decdist`, `Mhadron` | 0 and Mp (no decay) |

The arm prefixes follow `electron_arm` as for every layout (D1: SHMS e⁻, HMS p → `hs*` = proton). The Delta captures prove the
53 columns bit for bit (`l1_simcrun_replay_delta_*`).

### 1e. SIDIS layout (doing_semi; NtupleInit.f:192-237, results_write.f:187-215)

Columns 1-33 as the (e,e'p) layout, then:

| # | branch | source | units |
|---|---|---|---|
| 34 | `missmass` | `main.ntup.mm/1000` (complete_recon_ev, nucleon) | GeV |
| 35 | `ppi` | `recon.p.P/1000` | GeV/c |
| 36 | `t` | `main.ntup.t/1e6` (recon) | GeV² |
| 37 | `fry` | `-main.target.rastery` | cm |
| 38 | `radphot` | `main.ntup.radphot/1000` | GeV |
| 39 | `siglab` | `main.sigcc` = peepiX | µb/MeV²/sr² |
| 40 | `sigcent` | `main.sigcent`: never assigned, always 0 (D-051 F5) | |
| 41 | `Weight` | `main.weight` (× the survival probability for K± with doing_decay = 0) | |
| 42 | `decdist` | decay distance [cm] with doing_decay = 1; the survival probability (peepiX) with doing_decay = 0 | |
| 43 | `Mhadron` | `sqrt(Mh2_final)` | MeV |
| 44, 45 | `z`, `zi` | recon / vertex z = E_h/ν | |
| 46, 47 | `pt2`, `pt2i` | recon / vertex p_T² /1e6 | GeV² |
| 48, 49 | `xbj`, `xbji` | recon / vertex x_Bj (proton mass) | |
| 50 | `thqi` | `acos(vertex.uq.z)` (the vertex q angle) | rad |
| 51 | `sighad` | `main.ntup.sigcm` = peepiX σ_had | GeV⁻² |
| 52 | `jacobian` | `main.davejac` (peepiX jacobian × 1000) | |
| 53 | `centjac` | `main.johnjac`: never assigned, always 0 (D-051 F5) | |
| 54 | `pfermi` | `pfer/1000 * sign(p̂_fermi·q̂)`: fortran_compat 0/0 = **NaN** for every event with pfer = 0 (H, A > 2; D-051 F3), default 0 | GeV/c |
| 55 | `xfermi` | `main.ntup.xfermi` (x_Bj with Fermi motion; unchanged without do_fermi) | |
| 56 | `phipqi` | `main.phi_pq` | rad |

The two SIDIS captures prove the 56 columns bit for bit (`l1_simcrun_replay_semi_{h,d}`).

## 2. Run summary JSON (RunSummarySink, hist2json)

- Schema: [`docs/summary.schema.json`](summary.schema.json), JSON Schema draft-07, `schema = "simc_run_summary"`,
  `schema_version = 1`.
- The layout is the one read by the P1.5 comparison tool (validation/README.md, validation/compare/RunSummary.cpp).
- `run.producer` is `"simc++"` for the C++ and `"hist2json (Fortran .hist/.geni)"` for Fortran-derived files.
  simc_compare recognises "hist2json" and applies the `.hist` print precision.
- Units follow the `.hist` print-out.
- Acceptance: `simc_compare --a <run>.hist --b <run>.summary.json --expect IDENTICAL` passes for the 6 P1.2 runs and the
  2 P0 runs.

**Required keys** (these are what simc_compare reads):

| key | meaning | Fortran (.hist / .geni, simc.f) |
|---|---|---|
| `run.tag`, `run.producer` | run name, program | file stem |
| `run.random_seed` | seed (C++ always writes it; the Fortran prints it only when non-zero) | `Random Seed =` |
| `counters.ngen` | with the Fortran sign convention (ngen < 0 means \|ngen\| tries) | `Ngen (request)` |
| `counters.ntried`, `counters.ncontribute`, `counters.nco_no_rad_prot` | | `Ntried`, `Ncontribute`, `Nco_no_rad_prot` |
| `normalization.wtcontr` | wtcontribute·normfac / nevent | `MeV: wtcontr=` (simc.f:920) |
| `normalization.sigcc_central`, `normalization.sigcc_average` | central σ, Σσ/nevent | `CENTRAL.sigcc`, `AVERAGE.sigcc` |
| `normalization.charge_mC, targetfac, luminosity_ub, luminosity_GeV2, genvol, normfac` | | `charge`, `targetfac`, `luminosity` (ub^-1, then GeV^2), `genvol`, `normfac` |
| `flags.<name>` | every logical `name= T/F` of the FLAGS and RADIATIVE blocks, **named exactly as printed** (e.g. `using_E_arm_monteca`, `use_offshell_r`, which the Fortran truncates) | FLAGS, RADIATIVE CORRECTIONS |
| `stops.<spectrometer>.<counter>` | integer stop counters (`trials` and `successes` are required). The keys are listed below. | `.geni` |
| `tree.file, tree.name, tree.weight_branch` | the event tree. A relative `file` is relative to the JSON's folder; the C++ writes `<run>.root`, and hist2json writes the absolute path of `<run>.root` when it sits next to the `.hist`. | |

**Optional sections** (hist2json fills all of them; the C++ fills the ones marked *):

| key | content | .hist block |
|---|---|---|
| `run.reaction`*, `run.reaction_detail`, `run.begin_time`, `run.end_time`, `run.theory_file`, `run.random_state_file` | | header, MISCELLANEOUS |
| `counters.pct_no_rad_proton` | | `-> %no_rad_prot` |
| `settings`* | integer/real values of the FLAGS block: which_pion, which_kaon, pizero_ngamma, electron_arm, hadron_arm, ctau, transparency | FLAGS |
| `stop_roles`* | {spectrometer: electron/hadron/both} | from electron_arm/hadron_arm |
| `kinematics`* | Ebeam, dE_over_E, x/y_width, fr_pattern, fr1, fr2; `e`/`p`: angle_deg, momentum, offset_x/y/z, offset_xptar_mr, offset_yptar_mr | KINEMATICS |
| `central` | Q2_GeV2, q, nu, Em, Pm, W, MM; `e`/`p`: delta, xptar_mr, yptar_mr | VALUES FOR "CENTRAL" EVENT |
| `target` | A, Z, mass_amu, mass_MeV, Mrec_amu, Mrec_MeV, Mtar_struck, Mrec_struck, rho, thick_g_cm2, angle_deg, abundancy, X0, X0_cm, length, x/y/zoffset, musc_*; `losses.<name>` = [ave, lo, hi] | TARGET specs |
| `radiative` | the integer/real values (rad_flag, extrad_flag, one_tail, intcor_mode, dE_edge_test, Egamma_max) and `central` = the "Central Values" | RADIATIVE CORRECTIONS |
| `recon_summary.e/p.<delta, xptar, yptar, ytar>` | {ave_error, resolution, units} | RECON SUMMARY |
| `spedge`*.e/p.<delta, yptar, xptar> | [min, max] in % and rad. The Fortran prints these **after slop widening** (init.f:209-216); the C++ writes the converted input until the SetupBuilder exists, and says so in `notes`. | Input Spectrometer Limits |
| `limits.vertex/original/spectrometer/radiation` | rows {label, used:[min,max], found:[lo,hi], units} | the four "Limiting … values" tables |
| `slop.MC/total` | rows {label, used, min, max, units} | ACTUAL and LIMITING SLOP values |
| `notes`* | caveats | |

`RunSummary::extra` lets a later driver add or override any optional section; it is merged verbatim.

### Stop counters (per-aperture keys)

These names are shared with the P1.5 reader. `trials` and `successes` are always present.

| spectrometer | keys (in `.geni` order) | Fortran counters |
|---|---|---|
| hms | trials, slit_hor, slit_vert, slit_oct, Q1_in, Q1_mid, Q1_out, Q2_in, Q2_mid, Q2_out, Q3_in, Q3_mid, Q3_out, D1_in, D1_out, hut, dc1, dc2, scin, cal, successes (+ `coll`: not printed by the Fortran; the corrected SIMC++ .geni prints it as "Collimator stops (pi/mu)" when using_HMScoll is on, D-037 F2) | hSTOP_* (hms/mc_hms.f, struct_hms.inc) |
| shms | trials, hb_in, hb_men, hb_mex, hb_out, slit_hor, slit_vert, slit_oct, q1_in, q1_men, q1_mid, q1_mex, q1_out, q2_in, q2_men, q2_mid, q2_mex, q2_out, q3_in, q3_men, **q2_mid_dup**, q3_mex, q3_out, D1_in, D1_flr, D1_mid1…D1_mid7, D1_mex, D1_out, hut, dc1, dc2, s1, s2, s3, cal, successes (+ q3_mid, D1_men, targ_hor/vert/oct, cer1, cer2, BP_in, BP_out: not printed; coll: printed only by the corrected SIMC++ .geni ("Collimator stops (pi/mu)", using_SHMScoll, D-037 F2)) | shmsSTOP_* (shms/mc_shms.f, struct_shms.inc) |
| sos | trials, slit_hor, slit_vert, slit_oct, quad_in/mid/out, bm01_in/out, bm02_in/out, exit, hut, dc1, dc2, scin, cal, successes | sSTOP_* |
| hrsr, hrsl | trials, slit_hor, slit_vert, Q1_…, Q2_…, D1_in, D1_out, Q3_…, hut, dc1, dc2, s1, s2, cal, successes | rSTOP_*, lSTOP_* |
| calo | trials, slit_hor, slit_vert, successes | caloSTOP_* |

Fortran `.geni` label bugs (simc.f; docs/FORTRAN_DIFFERENCES.md D-011):
- **l.519-520**: the SHMS "Q3 … mid" column prints `shmsSTOP_q2_mid`. It is stored as `q2_mid_dup`; the C++ summary also
  writes `q2_mid_dup = q2_mid` so that both sides carry it (for simc_compare). The C++ `.geni` prints `q3_mid` in that column
  by default and `q2_mid` with `--fortran-compat`.
- **l.477, 489, 501**: the SOS "Slit hor/vert/corners" line and the HRS "Slit hor/vert" line print `slit_vert` first. The
  counters are stored under their physical keys (the C++ does not write SOS/HRS blocks yet, P14).
- **l.526**: SHMS "DC1, DC2, Scin, Cal" prints only dc1 and dc2. The C++ label is "DC1, DC2" by default, the Fortran one
  with `--fortran-compat`.

## 3. Diagnostic histograms (GenHistSink)

- Histograms `<set>_<arm>_<var>`, with set ∈ {geni, gen, recon}, arm ∈ {e, p}, var ∈ {delta, yptar, xptar}, plus `<set>_Em` and
  `<set>_Pm`.
- They have 50 bins (histograms.inc nHbins). The axes come from `RunSetup::hist` (init.f:limits_init l.515-570).
- They are filled with Fortran `inc()` semantics (simc.f:624-640): `ibin = nint(0.5+(val-min)/bin)`, counted only for 1 ≤ ibin ≤ 50.
  There is no under/overflow.

| set | when | values, weight | Fortran |
|---|---|---|---|
| geni | every tried event | vertex e/p delta, yptar, **−xptar**, Em, Pm; weight 1 | simc.f:254-261 |
| recon | successes | main.RECON e/p delta, yptar, xptar (**not** negated: FORTRAN-DIFF-CANDIDATE) with weight main.weight; recon Em, Pm with weight 1 | simc.f:272-279 |
| gen | successes | vertex e/p delta, yptar, −xptar, Em; weight 1 (**gen_Pm is never filled**) | simc.f:280-286 |

The Fortran `.gen` text print-out has its own quirks, recorded here for the comparison tool:
- The "Original E arm Events" block prints `H%gen` for delta but `H%geni` for yptar/xptar, both on the `H%gen` bin centres
  (simc.f .gen writer).
- The `.gen` has no `gen_Pm`.
