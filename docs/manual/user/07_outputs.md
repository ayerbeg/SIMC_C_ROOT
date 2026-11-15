# 7. Outputs

A run writes into the output directory (`-o DIR`, default `.`), with `<run>` = the input file name without its extension
(or `--run-name`):

| file | written by (sink) | content | Fortran counterpart |
|---|---|---|---|
| `<run>.root` | `root_tree` | tree `h10`, one entry per accepted event; only if `Nntu > 0` | `worksim/<run>.bin` → `make_root_tree` → `<run>.root` |
| `<run>.hist` | `fortran_hist` | the run report | `outfiles/<run>.hist` |
| `<run>.geni` | `fortran_hist` | spectrometer stop counters | `outfiles/<run>.geni` |
| `<run>.gen` | `fortran_hist` | diagnostic histograms (text) | `outfiles/<run>.gen` |
| `<run>.gen.root` | `gen_hist` | the same histograms as ROOT TH1D | — |
| `<run>.summary.json`, `.summary.txt` | `run_summary` | run summary (schema `docs/summary.schema.json`) | — (`hist2json` builds it from a Fortran `.hist`/`.geni`) |

`ngen = 0` writes nothing (simc.f:446). There is no `.bin` file. The detailed reference, with the Fortran source line of
every value, is `docs/OUTPUT_FORMAT.md`.

## 7.1 The event tree `h10`

* One entry per **accepted** event (success), as the Fortran (results_write.f:34).
* Every branch is a `Float_t` with the Fortran ntuple name and order (NtupleInit.f), computed in double precision and
  stored as float, as `make_root_tree` does. A Fortran tree and a SIMC++ tree can be compared branch by branch.
* The layout depends on the reaction (`NtupleColumnsFor`):

| run | columns | section |
|---|---|---|
| H/D/A(e,e'p) and `doing_phsp` | 46 | [7.1.3](#713-eep-layout-46-columns) |
| pion (`doing_pion`), Delta (`doing_delta`) | 53; π0: 65 | [7.1.4](#714-pion-kaon-and-delta-layouts) |
| kaon (`doing_kaon`) | 55 | [7.1.4](#714-pion-kaon-and-delta-layouts) |
| SIDIS (`doing_semi`) | 56; π0: 68 | [7.1.5](#715-sidis-layout-56-columns) |
| rho (`doing_rho`) | 59 | [7.1.6](#716-rho-layout-59-columns) |
| any of the above except (e,e'p), with `using_tgt_field` | + 8 | [7.1.7](#717-polarized-target-columns-using_tgt_field) |

### 7.1.1 Weights and normalisation

* `Weight` is the event weight computed by event.f:complete_main: the cross section of the event times the generation
  weight (radiative weight, jacobians, spectral-function or survival factors, depending on the reaction).
* The yield of one event for the simulated charge `EXPER%charge` is `Weight × normfac / nevent`. `normfac` is printed in
  the `.hist` and the summary, and `nevent` = Ncontribute when ngen > 0, Ntried when ngen < 0 (simc.f:341-347).
  `normfac = luminosity × genvol / Ntried × nevent`, and 1 for `doing_phsp` (simc.f:368-398).
* The `.hist` `wtcontr` is the sum of these yields over the accepted events that pass the `SPedge` window
  (pass_cuts): the expected number of counts.
* `sigcc` / `siglab` is the cross section of the event (units in the tables below).

### 7.1.2 Arm prefixes `hs*` and `ss*`

As in results_write.f, `hs*` is the spectrometer on the beam **right** and `ss*` the one on the **left**:

* `electron_arm` = 1, 3, 7 (HMS, HRS-R, calorimeter HMS side): `hs*` = electron, `ss*` = hadron;
* `electron_arm` = 2, 4, 5, 8 (SOS, HRS-L, SHMS, calorimeter SHMS side): `hs*` = hadron, `ss*` = electron.

Slopes (`*yptar`, `*xptar`, `*xpfp`, `*ypfp`) are in radians (dx/dz), even where the Fortran comments say mr; deltas are in %.

### 7.1.3 (e,e'p) layout (46 columns)

| # | branches | content | units |
|---|---|---|---|
| 1-4 | `hsdelta hsyptar hsxptar hsytar` | right arm, reconstructed | %, rad, rad, cm |
| 5-8 | `hsxfp hsxpfp hsyfp hsypfp` | right arm focal plane | cm, rad, cm, rad |
| 9-12 | `hsdeltai hsyptari hsxptari hsytari` | right arm, generated (at the vertex, before the spectrometer) | %, rad, rad, cm |
| 13-24 | `ssdelta ... ssytari` | the same for the left arm | |
| 25-29 | `q nu Q2 W epsilon` | reconstructed kinematics | GeV/c, GeV, GeV², GeV, — |
| 30-31 | `Em Pm` | missing energy and momentum | GeV, GeV/c |
| 32-33 | `thetapq phipq` | angles of the proton w.r.t. q | rad |
| 34 | `corrsing` | e' momentum minus the elastic value | GeV |
| 35-37 | `Pmx Pmy Pmz` | missing-momentum components ("Heepcheck" frame) | GeV/c |
| 38-40 | `PmPar PmPer PmOop` | missing momentum parallel, perpendicular, out of plane | GeV/c |
| 41 | `fry` | −(raster vertical position) | cm |
| 42 | `radphot` | radiated photon energy | GeV |
| 43 | `sigcc` | cross section (µb/sr for H, µb/MeV/sr² for D/A) | |
| 44 | `Weight` | event weight | |
| 45-46 | `theta_e theta_p` | reconstructed scattering angles | rad |

### 7.1.4 Pion, kaon and Delta layouts

Columns 1-33 as above, then:

| # | branch | content | units |
|---|---|---|---|
| 34 | `missmass` | missing mass (on the struck nucleon) | GeV |
| 35 | `mmnuc` | missing mass on the nucleus | GeV |
| 36 | `phad` | hadron momentum | GeV/c |
| 37 | `t` | −t (SIMC sign), reconstructed | GeV² |
| 38-40 | `pmpar pmper pmoop` | missing-momentum components | GeV/c |
| 41 | `fry` | −(raster vertical position) | cm |
| 42 | `radphot` | radiated photon energy | GeV |
| 43 | `pfermi` | Fermi momentum along q (signed) | GeV/c |
| 44 | `siglab` | lab cross section | µb/MeV/sr² |
| 45 | `sigcm` | CM cross section | model units |
| 46 | `Weight` | event weight | |
| 47 | `decdist` | decay distance (cm); for kaons with decay OFF the survival probability | |
| 48 | `Mhadron` | mass of the detected particle (μ or π after a decay) | MeV |
| 49 | `pdotqhat` | Fermi momentum · q̂ | GeV/c |
| 50-53 | `Q2i Wi ti phipqi` | vertex (generated) Q², W, t, φpq | GeV², GeV, GeV², rad |
| 54-65 | `xcal_gamma1 ycal_gamma1 Egamma1 Pgamma1x Pgamma1y Pgamma1z` and the same for photon 2 | π0 runs only: photon positions on the calorimeter (cm; −1e10 = missed, 0 with HMS/SHMS) and lab photon energy and momentum (MeV) | |

Kaon runs add 54 `saghai` (Saghai-model cross section; 0 in every event with `--fortran-compat`, D-043 F1) and 55 `factor`
(factorized model, = `sigcm`). Delta runs use the 53 pion columns: `missmass` is the missing pion mass, `phad` the proton
momentum, `siglab` and `Weight` contain **no cross section** ([5.5](05_reactions.md#55-eepπ-in-the-delta-region-doing_delta)),
and `decdist`, `Mhadron` are 0 and Mp.

### 7.1.5 SIDIS layout (56 columns)

Columns 1-33 as (e,e'p), then: 34 `missmass` (GeV), 35 `ppi` (hadron momentum, GeV/c), 36 `t` (GeV²), 37 `fry`, 38
`radphot`, 39 `siglab` (µb/MeV²/sr²), 40 `sigcent` (always 0), 41 `Weight`, 42 `decdist` (cm, or the survival probability
with decay OFF), 43 `Mhadron` (MeV), 44-45 `z zi`, 46-47 `pt2 pt2i` (GeV²), 48-49 `xbj xbji`, 50 `thqi` (rad), 51 `sighad`
(GeV⁻²), 52 `jacobian`, 53 `centjac` (always 0), 54 `pfermi` (GeV/c; NaN for H and A ≥ 3 with `--fortran-compat`, 0 in
the default mode, D-051 F3), 55 `xfermi`, 56 `phipqi` (rad). π0 SIDIS appends the 12 photon columns (57-68). The `i`
suffix means vertex (generated) values.

### 7.1.6 Rho layout (59 columns)

The 56 SIDIS columns with the rho meanings, then 57 `Mrho` (Breit-Wigner ρ mass of the event, MeV), 58 `Thrho` (decay
angle in the ρ rest frame w.r.t. the lab z axis, rad) and 59 `mmnuc` (GeV). The hadron-arm columns describe the detected
decay pion; `t` is the pion t; `siglab` (µb/MeV/sr) is the ρ cross section; `sigcent`, `zi`, `pt2i`, `xfermi` are 0.

### 7.1.7 Polarized-target columns (`using_tgt_field`)

Eight columns follow `phipqi` (pion/kaon/Delta: 54-61; SIDIS/rho: 57-64), in radians: `th_tarq` (angle between q and the
polarization, NaN when `targ_pol = 0`), `phitarq`, `beta`, `phis` (Sivers angle φpq − φtarg), `phic` (Collins angle φpq
+ φtarg), and the vertex values `betai`, `phisi`, `phici`. The kaon columns `saghai`/`factor` and the rho columns
`Mrho`/`Thrho`/`mmnuc` move after them. The shipped polarized-target pion example has 61 columns.

### 7.1.8 Values that are not measurements

* With an arm switched off (`spect_mode` ≠ 0), that arm's focal-plane columns and `*ytar` are 0.0 in every event
  ([6.4](06_spectrometers.md#64-spect_mode-switching-the-spectrometer-monte-carlo-off)).
* With `--fortran-compat`, some columns carry the Fortran's NaN or meaningless values: `pfermi` (SIDIS/rho, NaN), Delta
  `sigcm`, kaon `saghai` (0). [Chapter 12](12_faq.md#nan-or-constant-columns-in-the-tree) lists them.

## 7.2 `.hist`

The `.hist` is the Fortran report, in the same layout (a small gfortran-format emulator writes it): the run title and
time stamps, the FLAGS block, the kinematics, the target, the central event and its cross section (`CENTRAL.sigcc`), the
radiative settings, the event counters (`Ngen (request)`, `Ntried`, `Ncontribute`, `Nco_no_rad_prot`), the
normalisation (`wtcontr`, `AVERAGE.sigcc`, `charge`, `targetfac`, luminosity, `genvol`, `normfac`), the resolutions
(RECON SUMMARY), the input spectrometer limits (printed **after** the slop widening, as in the Fortran) and the limit
tables. In the default mode a few lines differ from the Fortran on purpose: some integers that the Fortran writes as T/F,
angle offsets in mr, and the CENTRAL.sigcc of pion/kaon/Delta runs ([chapter 10](10_differences_from_fortran.md)).
With RANLUX and `--fortran-compat` the `.hist` equals the Fortran's except lines 2-3 (BEGIN/END time).

## 7.3 `.geni`, `.gen` and `.gen.root`

* `.geni`: for each spectrometer used, the number of trials and where they were stopped (slits, magnet entrances,
  middles and exits, hut, chambers, scintillators, calorimeter), and the successes. Counter keys per spectrometer are
  listed in `docs/OUTPUT_FORMAT.md` §2. The default mode fixes Fortran label mistakes (SHMS Q3 mid printed q2_mid; SOS/HRS
  slit hor/vert swapped; SHMS "DC1, DC2, Scin, Cal" with two values) and adds a collimator line and HRS "Successes" lines.
* `.gen`: the Fortran diagnostic histograms (50 bins) of delta, yptar, xptar of both arms, Em and Pm: `geni` (all tries,
  vertex), `gen` (accepted, vertex) and `recon` (accepted, reconstructed, weighted). `.gen.root` has the same histograms
  as TH1D named `<set>_<arm>_<var>` (for example `gen_e_delta`, `recon_Em`), filled with the Fortran binning (no
  under/overflow). As in the Fortran, `gen` and `geni` use −xptar and `recon` uses +xptar, and `gen_Pm` is never filled.

## 7.4 Run summary (`.summary.json`, `.summary.txt`)

A JSON document (`"schema": "simc_run_summary"`, version 1, `docs/summary.schema.json`) with the run identification, the
counters, the normalisation, all FLAGS, the stop counters, the kinematics, the central event, the limits and the tree
location. `.summary.txt` is the same in readable form. `simc_compare` reads it ([chapter 9](09_comparing_with_fortran.md));
`hist2json` converts a Fortran `.hist`/`.geni` into the same format. A threaded run adds `"run": {"threads": N}`.

## 7.5 Exit status

| status | meaning |
|---|---|
| 0 | success |
| 1 | usage error (unknown option, missing input file, bad `--threads` value) |
| 2 | input or configuration error: unreadable file, rejected value (arm code 6, π0 from the neutron on H, rho on A ≥ 4, decay of a non-π/K hadron, `--threads` with `--fortran-compat` or RANLUX, ...) |
| 3 | run error: a Fortran STOP reproduced by SIMC++ (the message says `Fortran STOP ...`), a missing data file, an I/O error |

The Fortran exits with status 0 even after a STOP; SIMC++ returns 3 so that scripts can detect it (decision D-005).
