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
- Implemented layout: the (e,e'p) family (`doing_hyd_elast`, `doing_deuterium`, `doing_heavy`). This includes `doing_phsp`,
  which keeps `doing_hyd_elast`, so the acceptance runs also have these 46 columns. The pion/kaon/delta and semi/rho layouts
  (NtupleInit.f l.101-295) will come with those reactions. Until then `beginRun` throws for them.
- **Arm assignment** (results_write.f:64 and 89): `hs*` = the spectrometer on the beam **right** (HMS, HRSR) and `ss*` = the
  one on the **left** (SOS, HRSL, SHMS).
  - `electron_arm` = 1, 3, 7: the electron is on the right, so `hs*` holds the electron and `ss*` the hadron.
  - `electron_arm` = 2, 4, 5, 6, 8: `hs*` holds the hadron and `ss*` the electron.
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
| hms | trials, slit_hor, slit_vert, slit_oct, Q1_in, Q1_mid, Q1_out, Q2_in, Q2_mid, Q2_out, Q3_in, Q3_mid, Q3_out, D1_in, D1_out, hut, dc1, dc2, scin, cal, successes (+ `coll`, not printed) | hSTOP_* (hms/mc_hms.f, struct_hms.inc) |
| shms | trials, hb_in, hb_men, hb_mex, hb_out, slit_hor, slit_vert, slit_oct, q1_in, q1_men, q1_mid, q1_mex, q1_out, q2_in, q2_men, q2_mid, q2_mex, q2_out, q3_in, q3_men, **q2_mid_dup**, q3_mex, q3_out, D1_in, D1_flr, D1_mid1…D1_mid7, D1_mex, D1_out, hut, dc1, dc2, s1, s2, s3, cal, successes (+ q3_mid, D1_men, targ_hor/vert/oct, cer1, cer2, BP_in, BP_out, coll: not printed) | shmsSTOP_* (shms/mc_shms.f, struct_shms.inc) |
| sos | trials, slit_hor, slit_vert, slit_oct, quad_in/mid/out, bm01_in/out, bm02_in/out, exit, hut, dc1, dc2, scin, cal, successes | sSTOP_* |
| hrsr, hrsl | trials, slit_hor, slit_vert, Q1_…, Q2_…, D1_in, D1_out, Q3_…, hut, dc1, dc2, s1, s2, cal, successes | rSTOP_*, lSTOP_* |
| calo | trials, slit_hor, slit_vert, successes | caloSTOP_* |

FORTRAN-QUIRKs in the `.geni` (simc.f):
- **l.519-520**: the SHMS "Q3 … mid" column prints `shmsSTOP_q2_mid`. It is stored as `q2_mid_dup`; the C++ also writes
  `q2_mid_dup = q2_mid` so that both sides carry it. `q3_mid` is never printed by the Fortran.
- **l.477, 489, 501**: the SOS "Slit hor/vert/corners" line and the HRS "Slit hor/vert" line print `slit_vert` first. The
  counters are stored under their physical keys.
- **l.526**: SHMS "DC1, DC2, Scin, Cal" prints only dc1 and dc2.

## 3. Diagnostic histograms (GenHistSink)

- Histograms `<set>_<arm>_<var>`, with set ∈ {geni, gen, recon}, arm ∈ {e, p}, var ∈ {delta, yptar, xptar}, plus `<set>_Em` and
  `<set>_Pm`.
- They have 50 bins (histograms.inc nHbins). The axes come from `RunSetup::hist` (init.f:limits_init l.515-570).
- They are filled with Fortran `inc()` semantics (simc.f:624-640): `ibin = nint(0.5+(val-min)/bin)`, counted only for 1 ≤ ibin ≤ 50.
  There is no under/overflow.

| set | when | values, weight | Fortran |
|---|---|---|---|
| geni | every tried event | vertex e/p delta, yptar, **−xptar**, Em, Pm; weight 1 | simc.f:254-261 |
| recon | successes | main.RECON e/p delta, yptar, xptar (**not** negated: FORTRAN-QUIRK) with weight main.weight; recon Em, Pm with weight 1 | simc.f:272-279 |
| gen | successes | vertex e/p delta, yptar, −xptar, Em; weight 1 (**gen_Pm is never filled**) | simc.f:280-286 |

The Fortran `.gen` text print-out has its own quirks, recorded here for the comparison tool:
- The "Original E arm Events" block prints `H%gen` for delta but `H%geni` for yptar/xptar, both on the `H%gen` bin centres
  (simc.f .gen writer).
- The `.gen` has no `gen_Pm`.
