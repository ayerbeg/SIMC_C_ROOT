# Differences between SIMC++ and the Fortran SIMC

Status: task P7.5 (decisions D-031 and D-032, 2026-10-04); D/A(e,e'p) items added in P8.1/P8.2 (D-033). This file lists every known place where SIMC++ behaves
differently from the Fortran SIMC (`simc_gfortran`), or where it keeps a Fortran behaviour that looks odd. Each item
names the Fortran file and line, the C++ behaviour, the effect on outputs, and the test that covers it. The decisions
(D-xxx) are in `project_control/DECISIONS.md`. The Fortran reports (`reference/P*_report.md`) give the details.

## 1. Policy and the `fortran_compat` switch

* **Clear Fortran bugs are corrected by default** (D-031). Each correction carries a comment
  `// FORTRAN-BUG(D-xxx): original Fortran <file:line>: <code>; corrected because <reason>` in the C++.
* The run option **`fortran_compat`** restores the exact Fortran behaviour. It exists so that the bitwise tests
  (L1 tables, capture replays, regression baselines) can still prove the translation.
  * CLI: `simc++ <input> --fortran-compat`
  * JSON: `"run": {"fortran_compat": true}`
  * Code: `driver::RunOptions::fortran_compat` → `SetupOptions::fortran_compat` → `RunSetup::fortran_compat`. Every stage
    reads it from `StageContext::setup`; kernels get it explicitly (`SpectrometerAngles(..., fortran_compat)`,
    `ExtremeTripArgs::fortran_compat`, `proton::DeForestState::fortran_compat`, `optics::CosyEngine::setFortranCompat`,
    `McHms/McShms::setFortranCompat`, `optics::rotateHaxis/rotateVaxis(..., fortran_compat)`). There is no global.
  * The run log says which mode is used: ` SIMC++: fortran_compat = false (clear Fortran bugs corrected, ...)`.
* **Ambiguous cases** (it is not certain that the Fortran is wrong) keep the Fortran behaviour in both modes. They are
  marked `// FORTRAN-DIFF-CANDIDATE(D-xxx): ... Question: ...` and listed in §4 for the user and an expert.
* **Acceptance mode** (D-009, D-028) is an open question with the user (D-031 flag 1). It keeps the Fortran behaviour
  and is marked `// FORTRAN-OPEN(D-031 flag 1)`.
* **Faithful but odd behaviour** that is not a bug (precision, formats, guards, dead code) is marked `// FORTRAN-QUIRK`.

Classes used below:

| class | meaning | default (corrected) mode | `fortran_compat` |
|---|---|---|---|
| a | clear Fortran bug | corrected | Fortran |
| a-np | clear Fortran bug in a code path that is not ported yet | note only; fix when the path is ported | — |
| b | ambiguous, maybe intended | Fortran | Fortran |
| c | open question with the user (acceptance mode) | Fortran | Fortran |
| d | faithful Fortran behaviour, not a bug (documented) | Fortran | Fortran |
| dev | deliberate design deviation of SIMC++ (not about a Fortran bug) | SIMC++ | SIMC++ |

## 2. Class a: clear Fortran bugs, corrected by default

Effects are for the validation runs at seed 4242: P1.2 (`eep_h_elastic_{norad,rad}_A`, `acc_phsp_hms_shms_A`) and
P7.F1 (K2–K4). §6 has the L2 numbers.

| id | Fortran | original behaviour | C++ default | effect (outputs, size, configs) | test |
|---|---|---|---|---|---|
| D-006 pass_cuts | simc.f:237 | the hadron delta upper cut uses `SPedge%e%delta%max` | `SPedge%p%delta%max` (`driver::PassCuts`) | `pass_cuts` → npasscuts, `wtcontribute` (.hist wtcontr), resolution sums; with `hard_cuts` also the accepted events. Only events with recon hadron delta between the two limits are affected. K3 (HMS hadron at 423 MeV/c, limits 10 % vs 22 %): 42 of 200000 accepted events (norad) and 722 of 200000 (rad) no longer pass, .hist wtcontr −0.013 % / −0.230 % and the RECON SUMMARY resolutions change in the 4th digit (the Fortran trees give the same fractions: 27 and 755 events). P1.2, K2, K4: 0 events (no hadron recon delta in the window; in K2 the hadron limit is the larger one). | `fortran_diff_D006_pass_cuts` |
| D-006 kaon decay | shared/transp.f:158 | K → ππ decay in the first half of an element gives the product the kaon mass (`m_final = Mk`) | `Mpi`, as transp.f:251 and project.f:88 | kaons with `doing_decay` only (P10): mass and momentum of the decay product. No effect on H(e,e'p). | `fortran_diff_D006_kaon_decay_mass` |
| D-011 .geni Q3 | simc.f:519-520 | the SHMS "Q3 … mid" column prints `shmsSTOP_q2_mid` | prints `shmsSTOP_q3_mid` | `.geni` text only (the summary JSON keeps `q2_mid_dup` for the comparison with Fortran .geni files). 0 = 0 in the P1.2 runs. | `fortran_diff_D011_D018_report_formats` |
| D-011 .geni DC | simc.f:526 | label "DC1, DC2, Scin, Cal" with two values | label "DC1, DC2" | `.geni` text only (every SHMS run). | same |
| D-018 offsets | simc.f:815-818 | `xptar/yptar offset` printed in rad under the unit `mr` | printed in mr (×1000) | `.hist` text, only when a spectrometer angle offset is non-zero (0 in all validation runs). | same |
| D-018 L format | simc.f:897-898, 907, 927 | the integers `deForest_flag`, `use_first_cer`, `use_expon` written with the L edit descriptor (T/F; −1 and +1 both "T") | written as integers (i2, i2, i3) | `.hist` text of every run (2 or 3 lines). D-018 quoted l.896/906/925; the lines are 898/907/927. | same |
| D-018 delta title | simc.f:727-732 | the (e,e'p)π reaction title goes to unit 6, not to the report unit | written to the `.hist` | `.hist` text of Delta runs (P11). | same |
| D-014 deForest | physics_proton.f:84-86 | the `cos_phi` warning is tested after `cos_phi` was clamped to ±1, so it never prints | tested (and printed) before the clamp | warning text only; the cross section is unchanged. D/A(e,e'p): printed for the central event of the A(e,e'p) inputs of P8.F1 (cos_phi = 1 + 3.5e-8). | `fortran_diff_D014_deforest_cos_phi_warning` |
| D-015 F9 raddeg | shared/rotate_haxis.f:41, rotate_vaxis.f:40 | `raddeg = 0.017453292` (π/180 to 8 digits, rel. error 1.1e-8) | π/180 | the rotated aperture coordinates of the HMS dipole (±6°) and the SHMS apertures (±1.5°, …) move by ~1e-8 relative; the transported track is not changed. 0 aperture decisions changed in the 10 validation runs (trees bitwise identical, §7). | `fortran_diff_D015_rotate_raddeg` |
| D-017 extreme trip | target.f:501 | hadron section, pudding can / cryo17: `zz` from the ELECTRON angle `the%min` | from `thp%min` | `targ%Eloss(3)%max`, `teff(3)%max` (→ hadron edge limits, .hist "Eloss/teff lo-hi"). Only for `targ%can` 2/3 AND an angle > 45° (otherwise z is clamped to −L/2 either way). Example: can 2, hadron at 52–60°: 1.3565 → 1.3591 MeV. All validation runs use the beer can (can 1): no effect. | `fortran_diff_D017_extreme_trip_hadron_angle` |
| D-021 coll slits | hms/mc_hms_coll.f:83-104 | `hSTOP_slit_hor/vert/oct` incremented for each of the 20 steps outside the opening, also for rays that pass | at most once per ray and counter (the Fortran's unused `slit_*_flag` flags) | `.geni` HMS slit counters with `using_HMScoll` and pions/muons (P9). Example (2000 pions): 9746/3955/10944 → 637/263/713. Transport, RNG and acceptance unchanged. | `fortran_diff_D021_hms_coll_slit_counters` |
| D-028 F2 NaN | event.f:1641 | `dy = sqrt(1/cos_dtheta**2-1.-dx**2)` of a rounding-negative argument → NaN | argument clamped to ≥ 0 (dy = 0) | H elastic: the hadron yptar of a try can be NaN; such an event carries NaN in every hadron-arm recon/fp value but a finite weight (Fortran K2 norad: 1 of 200000 events). With seed 4242: 1 failed try in P1.2 rad and 1 in K3 rad (`.gen.root` geni_p_yptar +1 entry in the bin at 0), 0 accepted events in all 10 runs (the C++ K2 norad has no NaN event with this seed). | `fortran_diff_D028_spectrometer_angles_nan` (110 NaN rows of the L1 table) |
| D-033 F3 Benhar SF | sf_lookup.f:93 (`sf_lookup_diff`) | `SFd = SF/4/3.1415926535/(Pm*Pm)/5.0/20.0`: the interpolated BIN probability divided by the event's own Pm² and by fixed 5 MeV × 20 MeV/c widths: +Inf/NaN at Pm = 0, a 1/Pm² divergence below the first node (Pm < 9.9-10 MeV/c), widths wrong for benharsf_3mod/4 (dPm 19.733, dEm 1/5/10) | the differential SF of every node, sfval/(4π Pmval² dPm dEm) with the node's widths from the file, interpolated with the same bilinear scheme (`nuclear::SfTable::lookupDiff`) | A(e,e'p) with `use_benhar_sf` only: main%SF_weight, hence Weight and wtcontr. §8 has the measured effect (eep_fe56_bh_rad_A). | `fortran_diff_D033_F3_benhar_sf_diff` |
| D-033 F5 PmPer | event.f:926 | `vertex%PmPer = sqrt(max(0, vertex%Pm**2 - PmPar**2 - PmOop**2))` is computed before vertex%Pm is assigned (l.940/943), i.e. with the previous try's Pm | recomputed with this event's Pm after l.940/943 (`EepReactionBase::completeEv`) | D/A(e,e'p): vertex%PmPer only. It is not written to the ntuple (the tree has the recon PmPer) nor used by the weight: no effect on the outputs. | `fortran_diff_D033_F5_pmper_stale_pm` |
| D-033 F6 sf_lookup state | sf_lookup.f:160-184 | for Em == Emval(numEm) (and NaN Em) no branch assigns Em1/Em2/sf1/sf2: the SAVEd values of the previous call are used (NaN on the first call of a process) | Em == Emval(numEm) uses the last interval; no state is kept between calls; NaN Em gives 0 | Benhar SF (A(e,e'p) `use_benhar_sf`; A(e,e'π/K) generate_em, not ported yet): only events with Em exactly on the last node (measure zero with continuous Em). | `fortran_diff_D033_F6_sf_lookup_stale` |
| D-033 F7 heavy Mrec | event.f:944 | `vertex%Mrec = sqrt(vertex%Emiss**2-vertex%Pmiss**2)` without a guard: NaN Mrec/Em/Trec and success = .true. (a NaN Em even passes the VERTEXedge tests, comparisons with NaN being false) | complete_ev fails (success = false) when Emiss² − Pmiss² < 0 or NaN | A(e,e'p) at extreme edges only (He-4/C-12 table rows with ν → 0 or E_p = Ein); 0 events in the L2 runs (§8). | `fortran_diff_D033_F7_heavy_mrec_guard` |
| D-025/D-032 SHMS spectr | shms/mc_shms.f:335-347 | the "LSA tune" STOP is tested on the first call only | tested on every call (both modes) | none in valid runs. | `l1_shms_selftest` ("spectr != 5 → FortranStop") |
| D-032 arm code 6 | dbase.f:895-897/915-917 (STOP), but dbase.f:250/259, init.f:176/199, target.f:97/217, results_write.f:90, simc.f:510-511/905-906 still treat 6 as SHMS | 6 = retired SHMS LSA tune, the run STOPs in dbase | rejected when the input is read (InpReader, JSON `FromJson`, `BuildRunSetup`), both modes: "electron_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS"; exit status 2. One spectrometer table (`core/SpectrometerId.h`) is used everywhere; no other code knows 6. The shms plug-in key `"spectr"` was removed. | none for valid inputs. | `fortran_diff_D032_spectrometer_id`, `cli_arm_code_6_rejected` |

### Class a-np: clear bugs in paths that are not ported yet (fix when porting)

| id | Fortran | bug | C++ now | to do |
|---|---|---|---|---|
| D-018 delta STOP | dbase.f:780-788 | the (e,e'p)π announcement tests `doing_deutpi/doing_hepi` instead of `doing_deutdelta/doing_hedelta`, so valid D/A(e,e'p)π runs STOP | STOPs as the Fortran (both modes), `FORTRAN-BUG` note in SetupBuilder.cpp | P11: test the delta flags |
| D-018 Δ++ STOP | dbase.f:823-826 | "Can't do Delta+ production from neutron for pi0" also stops charged H(e,e'π−)Δ++ (`which_pion = 3` without `doing_pizero`; infiles/fpi12_deltapp.inp) | STOPs as the Fortran, note | P9: add the `doing_pizero` condition |
| D-006 tgt field | simc.f:139 | `ang_targ_parm = +1*sin(phi)*(targ_Bangle-theta)` (the e-arm line l.125 has `theta-targ_Bangle`) | `using_tgt_field` is NotPorted | when the target field is ported |
| D-011 SOS/HRS .geni | simc.f:477, 489, 501 | "Slit hor/vert" prints `slit_vert` first | the C++ writes no SOS/HRS block yet; the Fortran-.geni parser maps by meaning | P14 |
| D-012 F4 calo windows | target.f:80-131, 200-247 | no window branch for the calorimeter arms 7/8: stale `s_air, s_kevlar, s_mylar, forward_path` | `hallc_cryo` refuses arms 7/8 | P14 |

## 3. Class c: acceptance mode, open with the user (D-031 flag 1)

| id | Fortran | behaviour (kept) | effect | marker |
|---|---|---|---|---|
| D-009 elastic in phsp | dbase.f:217, 227-234; event.f:546 vs 709 | hydrogen + `doing_phsp` keeps `doing_hyd_elast`, so `complete_ev` uses elastic kinematics (not phase space), with sigcc = 1 and normfac = 1 | all H acceptance runs | `FORTRAN-OPEN` in Pipeline.cpp |
| D-009 radiation in phsp | dbase.f:227-228 | only `rad_flag` is zeroed; `using_rad` stays on: radiative tails are generated in acceptance mode | P1.2 phsp Egamma up to 241.9 MeV | `FORTRAN-OPEN` in SetupBuilder.cpp |
| D-028 F1 target mass | dbase.f:227-233 vs 293-295; init.f:411-413; radc.f:297 | phsp clears `doing_eep` before the masses are set, so `Mtar_struck = 0` and `Egamma_tot_max` is one proton mass too low | P1.2 phsp Egamma_max 241.858 instead of 1180.130 MeV; at 6.4 GeV (acc_phsp_K2) negative → 0 events | `FORTRAN-OPEN` in SetupBuilder.cpp |

## 4. Class b: ambiguous, Fortran behaviour kept (questions for the user and an expert)

| id | Fortran | behaviour (kept) | question | effect |
|---|---|---|---|---|
| D-033 F4 D jacobian | event.f:1020-1024 | D(e,e'p): the hadron 1/cos³(θ−θ0) jacobian is applied only for doing_heavy/pion/kaon/delta/semi, not for D, although the proton angles are generated in xptar/yptar and genvol contains dΩp (simc.f:388) | is the missing jacobian for D intended? | D(e,e'p) weights: the factor would be 1/cos³ ≈ 1 − 1.5(xptar²+yptar²), i.e. 0 at the centre to −1.5 % at the corners of the generated hadron solid angle (±60 × ±80 mr) |
| D-033 F1 SF fallback | init.f:839-851, dbase.f:598-616 | every A other than 2/12/56/197 silently uses c12.theory (benharsf_12.dat with use_benhar_sf): He-3, He-4, Ca-40 get the C-12 spectral function. SIMC++ keeps the default but prints a WARNING and lets the input name the file (`pipeline.reaction.spectral_function.file`) | which spectral function for He-3/He-4 (theory files exist only for 2/12/56/197)? | eep_he4_rad_A is "He-4 kinematics with the C-12 SF" |
| D-015 F7 HMS hodo 2 | hms/mc_hms_hut.f:537-551 | hodoscope-2 planes tested with the hodoscope-1 limits/offsets | copy-paste or the same active area on purpose? | HMS acceptance at the hut (scin trigger) |
| D-015 F7 HMS Cherenkov | hms/mc_hms_hut.f:523-529 | multiple scattering in the Cherenkov exit foil applied twice | double counting or a second foil? | HMS fp smearing (`mc_smear`) |
| D-015 F5 SHMS drifts | shms/mc_shms.f (zd_d1men 39.47, zd_d1out 60.68, zd_fp 307.95) | hard-coded drift lengths differ from the COSY classes (39.4715, 60.68058, 307.9514) | use the COSY values? | SHMS fp ≤ 2.2e-4 cm, path −0.0035 cm |
| D-017 can radius | target.f:368-376, 456-464 vs 108, 228 | beer-can corner of `extreme_trip_thru_target` uses 1.25", `trip_thru_target` 1.325" | which radius? | extreme Eloss/teff limits (beer can) |
| D-012 F4 th_can | target.f:236-244 | for \|z_can/R\| > 1 + 1e-6 the pudding-can `th_can` keeps its previous value (warning only) | reject, or clamp? | can 2 edge events |
| D-019 one_tail | radc.f:198-216, 234, 358 | with `one_tail ≠ 0` a disabled tail can be drawn: no photon, gen_weight = 1/hardcorfac | restrict `frac` to the enabled tails? | runs with one_tail ≠ 0 |
| D-019 doing_proton | init.f:693 | `doing_tail(3)` passed as `doing_proton` to lambda_dave: switching off tail 3 removes the ep interference from lambda(1:2) | intended? | one_tail runs |
| D-019 intcor_mode 0 | init.f:701-715 | `schwinger` does not set `dsoft_prime`; the SAVEd value (0) is used, g(4) loses its internal part (0.019443 vs 0.098637) | derivative term intended? | intcor_mode = 0 runs |
| D-019 rad_weight | radc.f:330-333, 403-406, 454-457 | `rad_weight*basicrad_weight` overwritten by `peaked_rad_weight` (rad_flag ≤ 1) | intended (peaked weight contains basicrad)? | radiative weights |
| D-019 tail-1 COMMON | radc.f:322-328, 512 | the `complete_ev` inside tail 1 recomputes hardcorfac/frac/g/c before the final division | post-radiation hardcorfac intended? | radiative weights |
| D-019 Q4 phi_ext | radc.f:566-642 | `rad_flag = 0` peaked weight independent of Egamma; phi_ext, dsoft_ext*, Friedrich terms unused | include phi_ext? | rad_flag = 0 |
| D-019 Q4 dsoft_int | radc.f:590-620 | rad_flag = 1 with extrad_flag = 3 uses the SAVEd dsoft_intmin/max (0) | compute them? | that flag combination |
| D-019 Egamma_min(1) | radc.f:283-286 | for ntail = 0 `Egamma_min(1)` keeps its stale COMMON value | set it? | D(e,e'p) rad_flag 3 (P8) |
| D-005 lambda input | init.f:692-695 | the input `lambda(1:3)` is overwritten every event | drop the input or honour it? | none (input ignored) |
| D-024 stale orig | event.f:400-415 | failed `generate_rad` leaves `orig` from the previous try | — (try rejected) | per-try outputs only |
| D-024 grnd dEbeam=0 | event.f:204 | a grnd is drawn for Ein even when dEbeam = 0 | — | random stream only |
| D-024 Eloss gaussians | event.f:193-205, trip_thru_target | energy-loss gaussians drawn even with using_Eloss = 0 | — | random stream only |
| D-015 F2 HMS resmult | hms/mc_hms_hut.f:292 | one grnd per hut entry even with smearing off | — | random stream only |
| D-015 tmpfact | simc.f:1334 | SAVEd `tmpfact`: when the electron arm stops before the hut, `ntup%resfac` gets the previous event's value | intended? | ntup resfac |
| D-015 F6 SHMS counters | shms/mc_shms_hut.f | S1 x and y both count in s1; 2x counts in s3, 2y in s2; `cal` counts two cuts (true position, fitted track 5 cm inside) | separate counters? | .geni diagnostics |
| D-011 recon xptar | simc.f:274, 277 vs 256, 262 | the .gen `recon` xptar histograms are filled NOT negated, `gen/geni` with −vertex xptar | which convention? | .gen / .gen.root histograms |
| D-015 F4 fry | simc.f:1467, 1487, 1733, 1753 vs 1472-1482 | mc_hms/mc_shms get `xtar_init` as fry, SOS/HRS get the raster fry | intended? | HMS/SHMS recon (fry terms) |
| D-024 H elastic W | event.f:723 block (pion/kaon only) | complete_ev never sets vertex W, epsilon, theta_pq, phi_pq for H elastic | intended? | vertex-only quantities |
| D-024 ngen<0 | simc.f:169-199, 341-350 | with ngen < 0 the .hist wtcontr and AVERAGE.sigcc are divided by tries | intended? | .hist with ngen < 0 |
| D-018 herho | dbase.f:213 | `doing_herho = A.eq.3` (the other he* flags use .ge.3) | intended? | rho on A > 3 (P13) |
| D-005 sigc_nbin | semi_physics.f:173 | default `sigc_nbin = 0` → division by zero | required input? | SIDIS (P12) |

## 5. Class d: faithful Fortran behaviour, not a bug (kept, `FORTRAN-QUIRK`)

| id | Fortran | behaviour | why kept |
|---|---|---|---|
| D-018 CTP float | CTP/thTestParse.c:392-393 | every real .inp literal is stored as a 32-bit float | defines the run inputs |
| D-010 literals | -fdefault-real-8; physics_proton.f | double literals, `pow` for real exponents, x+y+z dot products | bitwise numerics |
| D-024/D-018 powers | gfortran | `x**2` → x·x; `x**n` (integer n) → `__powidf2` (`optics::powi`) | bitwise numerics |
| D-015 F1 / D-021 lfit | cern/lfit.f, mc_*_hut.f | focal plane real*4, accumulators real*8, a hit at exactly 0.0 ignored (KEY = 0) | CERNLIB behaviour |
| D-015 F3 closure | COSY maps | forward→recon closure not exact (HMS xptar RMS 1 mr) | optics |
| D-015 fry | mc_hms_recon.f:118, mc_shms_recon.f:136 | fry = 0 replaced by 1e-30 | avoids 0**0 |
| D-019 gamma | radc.f:100-115 | the polynomial `gamma()` is dead code (callers bind to libm tgamma); wrong for x ≥ 2 | `std::tgamma` used, as the binary |
| P5.F1 Q6, Q8-Q10 | radc.f, brem.f | basicrad itail overwrite, schwinger's exp convention (only dhard used), brem/bremos kinematic simplifications, spen ≤ 101 terms, one-time warnings, a TAB in a message | no effect / physics model |
| D-012 F2 | enerloss_new.f, gauss1.f | gauss1 = 0 for u1 = 0.5 exactly; −2 ln 0 replaced by 1e5 | guard; rate ~0 with the 53-bit RNG |
| D-012 F3 | enerloss_new.f | invalid typeflags keep `x` (not used by SIMC); Eloss > KE clipped to KE − 1e-7 | guard |
| D-012 F5 | target.f, musc*.f | gauss1 truncations 3.5σ / 99σ / 10σ | explicit choices |
| D-012 F6 | gauss1.f | s = 1 exactly accepted (g = ±0) | negligible |
| D-012 F7 | event.f:1572-1616 | physics_angles assumes phi0 = π/2 or 3π/2 (warns otherwise) | all spectrometers satisfy it |
| D-012 F1 | event.f:1641 | dy ≈ 3e-9 instead of 0 for a true dy = 0 (the NaN part is corrected, §2) | rounding |
| D-024 θ = 0 | event.f:532-538 | elastic θ = 0 can round 1 ulp above Ein (failure); NaN at q = 0 does not STOP | no spectrometer at 0° |
| D-024 Pm | event.f:546-560 | H elastic: vertex Pm = 0, so PmPer = 0 | definition |
| P4.F1 Q4 | event.f:325 | the hadron physics_angles call uses the previous try's xptar/yptar (overwritten by complete_ev) | no effect |
| D-008 phsp | event.f:1132-1136, 1446-1448; simc.f:397 | phsp: recon hadron P set to the central momentum, sigcc = 1, normfac = 1 | phsp design |
| P4.F1 Q6 | event.f:1060-1140 | complete_recon_ev: W with Mp, phi_pq with acos (complete_ev: atan2), recon Ein assigned twice | definitions |
| D-005 .hist limits | init.f:209-216 | "Input Spectrometer Limits" printed after the slop widening | output definition |
| D-005 tree | make_root_tree.f | every branch Float_t; ≤ 1e7 events | output format |
| D-033 F6 logsf | sf_lookup.f:184 | the variable is called `logsf` but the Em interpolation is linear (no log is taken); SF < 1e-20 (also negative extrapolations, and one negative node of benharsf_4.dat) is set to 0 | the coded algorithm |
| D-033 F8 krel | event.f:977 | `ntup%krel` is written from the SAVEd local krel, never assigned on the (e,e'p) branches (0) | not an output of (e,e'p) |
| D-033 F8/F9 Benhar limits | init.f:335, dbase.f:592-621 | with use_benhar_sf, theory_init is not called: E_Fermi = 0 (VERTEXedge%Em%min = 0) and VERTEXedge%Pm%max = 790 MeV/c for every file (He-3/He-4 tables go to 976.8); the neutron columns of benharsf_4/12/56/197 normalise like the proton ones | definitions of the input files |
| D-033 F12 D warning | event.f:588-593 | the "two solutions" warning of D(e,e'p) is printed only for ntried ≤ 5000 (`IReaction::beginTry`) | message only |
| D-033 doing_deuterium | init.f:889 | theory_init overwrites doing_deuterium = (nrhoPm = 1 and E_Fermi < 1) from the file contents | definition |
| P8.2 targetfac | simc.f:94-95 | `3.75914d+6` is DOUBLE PRECISION, promoted to real*16 by -fdefault-real-8: targetfac is evaluated in quad precision (objdump: __divtf3/__multf3) and rounded to real*8; SIMC++ does the same (`__float128`); the double evaluation differs by 1 ulp for the D capture input | bitwise numerics |
| D-015 HMS cal | hms/mc_hms_hut.f:582-588 | calorimeter cut commented out: hSTOP_cal never incremented | deliberate in the Fortran |
| D-015 SHMS counters | shms/mc_shms*.f | 7 counters never incremented; 5 stop reasons unreachable; stop_id reused / 0 for slit stops | diagnostics |
| P6.F1 F8, F10-F12 | mc_shms(_hut).f, transp.f | `spec` real*4/real*8 mismatch (unused), unit 1 closed, TOF terms dropped, pathlen to the back of the calorimeter | no effect / definitions |
| P3.F1 | target.f:62-151 | `goto (10,20,30) narm` falls through to the beam for other narm | unreachable |
| limits_update | event.f limits_update | `orig%e%xptar` range updated twice | harmless duplicate |
| .geni histograms | simc.f:578-586 | geni buffers printed on the H%gen bin centres | the axes are the same (init.f:515-570) |
| calculate_central | simc.f:1268-1281 | `main0` is SAVEd: unset fields are 0 | definition |
| montecarlo | simc.f:1800-1842 | RECON z not written with the arm MC off; recon e E not Coulomb-corrected | definition |
| P1.6 | physics_proton.f | real(16) ABI of real*8 functions; fofa for timelike Q²; deForest at θpq = 0 by rounding | no effect / numerics |

## 6. Deliberate design deviations (class dev)

| id | Fortran | SIMC++ |
|---|---|---|
| D-002 RNG | RANLUX (call_ranlux.f; skips every 1000th number) | `mt64` (std::mt19937_64); `random_state_file` ignored; validation is statistical. Bitwise replays use `ScriptedRandom` with captured streams. |
| D-005 STOP | `STOP` exits with status 0 | `simc::FortranStop`, exit status 3; rejected inputs exit 2 |
| D-003 input | CTP .inp, `%` names only (obsolete `.` names silently ignored) | .inp (same CTP semantics, unknown names reported like CTP) or JSON; `inp2json` |
| D-032 | arm codes as raw integers, 6 half-supported | one `SpectrometerId` table; 6 rejected (see §2) |
| not ported | pion, kaon, Delta, rho, SIDIS, SOS/HRS/calorimeter arms, target field | `simc::NotPorted` at start-up (P9–P14) |
| D-033 F1 | the spectral-function file follows from nint(A) only | also selectable in the JSON (`pipeline.reaction.spectral_function`: type and file); a warning when the Fortran C-12 fallback is used |
| outputs | .hist/.geni/.gen, .bin, .root | the same text files, `.root` (h10, same 46 branches), plus `.summary.json/.txt` and `.gen.root`; no `.bin` |

## 7. Validation of the corrected default against the Fortran (L2, seed 4242, task P7.5)

`simc_compare` (D-011 verdict) of the default (corrected) mode against the Fortran references `reference/runs/<run>`
(A runs); outputs in `reference/comparisons/CPP_vs_F/P7.5_corrected_<run>/` (summary.txt, summary.json, compare.pdf).
"vs compat" compares the corrected run with the exact-Fortran-mode run of the same seed (baselines f9cea51_*).

| run | verdict (failed/N) | worst p | max \|pull\| | tree vs compat | text/histogram differences vs compat |
|---|---|---|---|---|---|
| eep_h_elastic_norad_A | PASS 1/206 | 1.3e-3 (phipq ks_w) | 1.03 | bitwise identical | .hist 2 format lines, .geni DC label |
| eep_h_elastic_rad_A | PASS 0/214 | 6.2e-2 | 1.71 | bitwise identical | format lines; geni_p_yptar +1 entry (NaN try, D-028) |
| acc_phsp_hms_shms_A | PASS 1/208 | 3.0e-3 (Pmx ks_w) | 1.23 | bitwise identical | format lines |
| eep_h_K2_norad_A | PASS 0/210 | 6.8e-2 | 1.82 | bitwise identical | format lines |
| eep_h_K2_rad_A | PASS 2/219 | 6.9e-3 | 2.15 | bitwise identical | format lines |
| acc_phsp_K2_A | PASS 0/12 (0 events in both, D-028 F1) | — | — | identical (empty) | format lines |
| eep_h_K3_norad_A | PASS 0/209 | 3.0e-2 | 1.74 | bitwise identical | wtcontr −0.013 %, resolutions (D-006 pass_cuts: 42 events) |
| eep_h_K3_rad_A | PASS 0/219 | 2.1e-2 | 1.88 | bitwise identical | wtcontr −0.230 %, resolutions (722 events); geni_p_yptar +1 (D-028) |
| eep_h_K4_norad_A | PASS 0/215 | 2.5e-2 | 2.23 | bitwise identical | format lines |
| eep_h_K4_rad_A | PASS 1/221 | 7.1e-3 | 2.69 | bitwise identical | format lines |

## 8. D/A(e,e'p): L2 in both modes and the measured effect of the corrections (task P8.2)

`simc_compare` (default config) of simc++ seed 4242, 200000 events, against the P8.F1 Fortran runs
`reference/runs/eep_<tag>_A` (and `_B` for D and C-12). Run folders `SIMC++_Project/simc_runs/2026100{4,5}_*_P8.2_L2_{compat,corrected}_<tag>/`
(`compare_vs_A/`, `compare_vs_B/`: summary.txt, summary.json, compare.pdf; `mode_diff_vs_compat.txt` in the corrected
folders: event-by-event comparison with the compat run, `validation/compare/mode_diff.cpp`).

| run | compat vs A: verdict (F/N), worst p | corrected vs A: verdict (F/N), worst p, wtcontr pull | corrected vs compat (same seed) |
|---|---|---|---|
| eep_d_rad | PASS 4/230, 4.0e-3 (vs B: PASS 1/229) | PASS 4/230, 4.0e-3, 0.74 (vs B: PASS 1/229) | trees identical; .hist wtcontr −0.504 % |
| eep_d_norad | PASS 1/227, 7.0e-3 | PASS 1/227, 7.0e-3, 1.43 | trees identical; wtcontr −0.607 % |
| eep_c12_rad | PASS 1/226, 7.1e-3 (vs B: PASS 0/226) | PASS 1/226, 7.1e-3, −0.09 (vs B: PASS 0/226) | trees identical; wtcontr −1.995 % |
| eep_fe56_rad | PASS 1/227, 4.7e-3 | PASS 1/227, 4.7e-3, 2.34 | trees identical; wtcontr −2.245 % |
| eep_fe56_norad | PASS 0/223, 2.0e-2 | PASS 0/223, 2.0e-2, 0.94 | trees identical; wtcontr −0.578 % |
| eep_au197_rad | PASS 3/226, 4.1e-4 | PASS 3/226, 4.1e-4, 1.94 | trees identical; wtcontr −2.468 % |
| eep_he4_rad (C-12 SF) | PASS 3/230, 2.5e-4 (hms.scin binomial, pull −3.7) | PASS 3/230, 2.5e-4, −0.14 | trees identical; wtcontr −1.789 % |
| eep_fe56_bh_rad | PASS 0/226, 1.2e-2 | PASS 0/226, 1.1e-2 (wtcontr), 2.55 | only Weight differs (199992 of 200000 events, ΣWeight −0.079 %); wtcontr −2.795 % |

* The corrected trees equal the compat trees bit for bit except the Benhar Weight: F5 (vertex PmPer) is not an output, F7
  and F6 occur in 0 events (same Ntried, same events), D-028 F2 in 0 events.
* The wtcontr difference is the pass_cuts correction (D-006, simc.f:237): the Fortran counts hadron recon δ up to the
  electron-arm limit (22 %) instead of 10 %. The compat/corrected ratios (+0.51, +0.61, +2.04, +2.30, +0.58, +2.53, +1.82,
  +2.80 %) equal the effect measured on the Fortran trees in P8.F1 (+0.52, +0.60, +2.0, +2.2, +0.53, +2.5, +2.1, +2.8 %),
  the remaining small differences being the statistics of the window 10-22 %. The corrected-vs-Fortran wtcontr pulls are
  therefore shifted by ≈ effect/relerr: up to 2.3-2.6 σ for Fe and Fe-Benhar, still below the p < 0.01 test threshold.
* Benhar (F3): the corrected weight is the per-node density interpolated (no 1/Pm² of the event); per recon-Pm bin the
  corrected/compat ratio is 0.960 (Pm < 50 MeV/c), 0.99 (50-150), 1.005-1.012 (150-350), 1.002 above; in total −0.08 %.
