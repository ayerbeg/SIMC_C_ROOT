# Differences between SIMC++ and the Fortran SIMC (user guide)

The complete, line-by-line list is `docs/FORTRAN_DIFFERENCES.md`. Each item there gives the Fortran file and line, the
original behaviour, the C++ behaviour, the measured effect and the test that covers it. This guide summarises it for
users: **what changes in your numbers, and when**. Everything in 10.1 is undone by `--fortran-compat`.

## 1. Corrections of the default mode (class a)

### Changes that can move yields, weights or acceptances

| id | when it matters | what changes in the default mode |
|---|---|---|
| D-006 pass_cuts | every run where the reconstructed hadron delta can lie between `SPedge%p%delta%max` and `SPedge%e%delta%max` | the Fortran tests the hadron delta against the **electron** delta maximum. The default mode uses the hadron limit. This changes `wtcontr`, the resolutions, and with `hard_cuts` the accepted events: −0.01 % to −0.23 % in the K3 H(e,e'p) runs, −24 % to −34 % of ΣWeight in rho runs, up to −24 % in spect_mode runs, 0 when the hadron window lies inside the electron one |
| D-043 F5, D-051 F6 survival momentum | kaon and SIDIS runs with `doing_decay = 0` | the survival probability uses the hadron's own momentum, not the central momentum: the kaon **weights** (and K± SIDIS weights) and the `decdist` column |
| D-043 F6, D-051 F6 aerogel distance | decay-OFF kaon/SIDIS runs with an HMS hadron (path length with the missing square root); SHMS unaffected | survival probability |
| D-061 F1, F1e HRS | every HRS-L/HRS-R run | the intended transport lengths (Fortran: path length −3e232 cm, "decays" everywhere, kaon survival weight +Inf) and the 45° VDC rotation (Fortran: no rotation). HRS acceptance, FP path, decay |
| D-048 F2 Delta tail 1 | radiative Delta runs | the incoming electron radiates (the Fortran fails every such try): Ntried, wtcontr |
| D-048 F1 Delta D/A | Delta runs on D or A | they run (the Fortran always STOPs); there is no Fortran reference for them |
| D-033 F3 Benhar SF | A(e,e'p) with `use_benhar_sf` | the proper differential spectral function of each node instead of a bin probability divided by the event's own Pm² (divergent at small Pm): weights, wtcontr |
| D-051 F1, F2 SIDIS | SIDIS events with all PDFs 0, or D/A runs with `do_fermi` and a negative x_Bj | weight 0 instead of NaN; the run no longer STOPs |
| D-051 F4 π0 SIDIS | π0 SIDIS | the d̄ fragmentation function and the π0 mass |
| D-057 F1, F6, F9, F21 rho | rho runs | energy loss of the ρ with its own angles; efer of D/³He; yptar sign of far-side pions; a NaN pion (1 event per 1e5) |
| D-057 F2 rho decay | rho runs | no decay after a failed try: no weight change, but the **random stream** differs from the Fortran's |
| D-037 F4, F7 collimator / decay | pion/kaon runs with collimators or decay | a particle that ranges out in the collimator stops; the track after a decay ends at the right z (≤ 0.28 cm in the path length, no acceptance decision changed in the captures) |
| D-062 G1, G2, G3 π0 / calorimeter | π0 runs with `doing_decay = 1`; calorimeter as electron arm with decay; π0 photon path | photons (and an electron) never decay, so the run works (the Fortran STOPs); each photon path starts at 0 |
| D-068 π0 in HMS/SHMS with decay | π0 transported in the HMS/SHMS with `doing_decay = 1` | the π0 is stable (the Fortran STOPs) |
| D-067 T1, T2 target field | `using_tgt_field` with `targ_Bphi` not on an arm's side, or field on the hadron side with targ_Bangle < θp | the field angles of the arms |
| D-033 F7, D-028 F2 NaN guards | rare edge events (A(e,e'p) extreme edges; H elastic rounding) | NaN events are rejected or clamped instead of carrying NaN values |
| D-015 F9 raddeg | HMS/SHMS rotated apertures | π/180 exactly instead of 8 digits (no decision changed in the validation runs) |
| D-017 extreme trip | pudding can or cryo17 with a hadron angle > 45° | the hadron energy-loss limits (e.g. 1.3565 → 1.3591 MeV) |

### Changes in tree columns only (no effect on weights)

| id | column | change |
|---|---|---|
| D-051 F3 | SIDIS/rho `pfermi` | 0 instead of NaN for H and A ≥ 3 |
| D-043 F1, F2, F3 | kaon `saghai` | the Saghai cross section (the Fortran writes 0 in every event; its tables were also read off by one line and on a wrong s grid) |
| D-048 F3 | Delta `sigcm` | a meaningful pion fit value (the Fortran's is garbage, NaN in 0.4 % of the events) |
| D-067 T14 | π0 and rho runs with the target field | each value in its own column (the Fortran overwrites the polarized-angle columns) |
| D-006 kaon decay | `Mhadron` after a K → ππ decay in the first half of a magnet | the pion mass instead of the kaon mass |
| D-038 F10, D-033 F5, D-051 F9, D-057 F3 | internal values (`sigcm2` above W = 2 GeV, vertex PmPer, Mh2_final) | not written to any output |

### Changes in the text reports only

| id | file | change |
|---|---|---|
| D-038 F5, D-043 F13, D-048 F9 | `.hist` `CENTRAL.sigcc` of pion, kaon and Delta runs | computed at the real central kinematics (the Fortran's is ×4.9 to ×11 off). A report value; not used anywhere else |
| D-011, D-060, D-061 F8, D-037 F2 | `.geni` | SHMS "Q3 mid" shows q3_mid (Fortran: q2_mid); SOS/HRS slit hor/vert in label order; the SHMS "DC1, DC2" label; HRS "Successes" lines; a "Collimator stops (pi/mu)" line |
| D-021, D-037 F1 | `.geni` HMS/SHMS slit counters with collimators | each slit counter counted once per ray (Fortran: once per collimator step) |
| D-018 | `.hist` | integers written as integers (Fortran: T/F); angle offsets in mr; the Delta title |
| D-038 F15, D-051 F12, F14, D-057 F15, F18 | `.hist` titles; rho `.geni`/`.gen`/`.hist` P-arm lines | titles for coherent π−, π0 SIDIS, A ≥ 3 kaon SIDIS (the Fortran STOPs), D/³He rho; rho P-arm histograms filled with the decay pion |
| D-014, D-048 F6 | log | the de Forest cos φ warning is printed; the Delta warning text |

### Changes in what is accepted as input

| id | input | default mode |
|---|---|---|
| D-032 | arm code 6 | rejected (both modes) |
| D-038 F1 | H(e,e'π−)Δ++ (`which_pion = 3`) | runs (the Fortran STOPs) |
| D-038 F2 | π0 from the neutron on hydrogen | refused (the Fortran runs it with a neutron mass on a proton target) |
| D-037 F9 | `doing_decay` with a hadron that is not π/K | refused (the Fortran STOPs silently at the first decay) |
| D-025/D-032 | SHMS LSA tune | checked at every call (no valid input affected) |

## 2. Kept as in the Fortran: acceptance mode (class c)

The user has not decided this yet (D-031 item 1), so it is kept in both modes. On hydrogen, `doing_phsp` keeps the
elastic kinematics; radiation stays on; the target mass for the photon limits is 0
([5.8](../manual/user/05_reactions.md#58-phase-space-and-acceptance-runs-doing_phsp)).

## 3. Questions open for experts (the expert list)

These Fortran behaviours **might** be intentional. SIMC++ keeps them in both modes (`docs/FORTRAN_DIFFERENCES.md` §4).
In plain words:

* **Delta weight without cross section** (D-048 F4): Delta yields are not rates. A related question is whether a factor
  2 is missing when both proton-energy roots are physical (F5).
* **Radiative corrections** (D-019): a disabled tail can still be drawn with `one_tail`; switching off the hadron tail
  also removes the e-p interference term; with `intcor_mode = 0` part of a term is lost; the input `lambda` is
  overwritten (D-005); a few stale values for specific flag settings.
* **D(e,e'p)** has no hadron-angle jacobian (≤ 1.5 %), A(e,e'p) does (D-033 F4). He-3/He-4 use the C-12 spectral
  function (F1).
* **Spectrometer details**: HMS hodoscope 2 tested with hodoscope-1 limits, a Cherenkov foil counted twice, SHMS drift
  lengths that differ from the COSY values (D-015); the beer-can radius differs between two routines (D-017); the SOS last
  7.62 cm has no path length or decay (D-060); SOS/HRS reconstruct with the raster position, HMS/SHMS without (D-061 F3);
  the HRS aerogel distance (F6).
* **Collimators and decay**: the collimator scattering ignores `mc_smear`; the pion absorption parameters; where a decay
  inside a magnet is applied (D-037 F3, F5, F8).
* **Pion and kaon models**: coherent production uses nucleon models; the input final-nucleus mass is ignored for coherent
  production; a new Δ mass when tail 1 recomputes the event; single-precision MAID amplitudes (D-038 F3, F4, F9, F12); a
  fixed kaon decay model, the Saghai angles, the Saghai table records, the proton mass for Σ− (D-043).
* **SIDIS**: the missing-mass threshold formula and its use for kaons and with Fermi motion, the generated hadron range
  with radiation, R_h for deuterium (D-051 F7, F8, F10, F16).
* **Rho**: no Fermi motion for D/³He, s without Fermi terms, t at the pole mass, the decay frame, the xptar formula, the
  radiation of the neutral ρ (70 % of tail-3 tries fail) and a new ρ when tail 1 recomputes the event (D-057).
* **Calorimeter and π0**: the calorimeter windows reuse the previous trip's values (D-012 F4); `drift_to_cal` 0 is not
  checked (D-062 G5); a second photon pair in tail 1 (G8); a π0 in the HMS/SHMS is transported as a charged particle (D-065
  P1) and is bent in the target field (D-067 T13).
* **Target field**: the step limit depends on the arm side, the convergence test, the 90/E constant, a polarization
  vector without `targ_Bphi` (D-067 T4-T6, T15).
* **Bookkeeping**: with ngen < 0 wtcontr is divided by the tries; `doing_herho` only for A = 3; SAVEd values in a few
  report quantities (D-024, D-015 tmpfact, D-018).

## 4. Deliberate differences of SIMC++ (class dev)

* Default generator mt64 instead of RANLUX; RANLUX available ([8.2](../manual/user/08_running_and_reproducibility.md#82-random-number-generators-mt64-and-ranlux)); `random_state_file` ignored.
* STOP → exit status 3; rejected input → 2 (the Fortran always exits 0).
* JSON input as a superset of the `.inp`; plug-in selection; `--threads`.
* Outputs: the ROOT tree written directly, plus the summary JSON and `.gen.root`; no `.bin`.
* The model tables are read from the data directory at start-up and md5-checked by the tests; a missing file stops the run
  before it begins.
* The D-033 F1 spectral-function fallback warns and can be overridden from the JSON.
* rho on A ≥ 4 and rho with another reaction flag or with `doing_pizero` are refused at input.
