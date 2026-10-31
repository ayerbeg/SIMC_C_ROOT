# Differences between SIMC++ and the Fortran SIMC

Status: task P7.5 (decisions D-031 and D-032, 2026-10-04); D/A(e,e'p) items added in P8.1/P8.2 (D-033); collimator
and decay-in-flight items added in P6.4 (D-037); pion items in P9 (D-038); kaon items in P10 (D-043); Delta items in P11 (D-048); SIDIS cross-section items in P12.1 (D-051 Fn = the P12.F1 findings; `P12.1 Qn` = further quirks of the port); rho items in P13.1/P13.2 (D-057 Fn = the P13.F1 findings, reference/P13.F1_report.md; D-056 maps the provisional P13.1 R ids). This file lists every known place where SIMC++ behaves
differently from the Fortran SIMC (`simc_gfortran`), or where it keeps a Fortran behaviour that looks odd. Each item
names the Fortran file and line, the C++ behaviour, the effect on outputs, and the test that covers it. The decisions
(D-xxx) are in `project_control/DECISIONS.md`. The Fortran reports (`reference/P*_report.md`) give the details.

> **Delta region (`doing_delta`, H/D/A(e,e'p)π): the event weight contains NO cross section** (D-048 F4, class b, kept in
> both modes). physics_delta.f:132 returns `peedelta = 1.0*jacobian*(gtpr*fac)`; the line with the model cross section
> (l.130) is commented out ("Need new xsec model.", event.f:1512). The Fortran (and SIMC++) Delta weight is phase space ×
> virtual-photon flux × the dt dφ_cm → dΩ_lab jacobian of the proton. Yields, `.hist` wtcontr and the `Weight` column of a
> Delta run are therefore NOT physical rates. Whether SIMC++ should offer a real Delta cross section is an open question
> for the user/an expert (§4).

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
| D-011 .geni SOS slit | simc.f:477 | the SOS "Slit hor/vert/corners" line prints `sSTOP_slit_vert, sSTOP_slit_hor, sSTOP_slit_oct` (horizontal and vertical swapped against the label and the HMS/SHMS lines) | prints slit_hor, slit_vert, slit_oct (`FormatGeni`, P14.1) | `.geni` text of SOS runs only (eep_h_test_pctnames: "1313 18 21" → "18 1313 21") | `l1_simcrun_ranlux_sos_eep_h_test_pctnames` (both modes) |
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
| D-038 F10 sigcm2 | physics_pion.f:133-154 | `ntup%sigcm2` (the MAID cross section) is assigned only for W < 2 GeV; for W ≥ 2 GeV it keeps the value of an earlier event (or of the central event) | `ntup%sigcm2 = 0` for W ≥ 2 GeV (MAID not evaluated) in `pion::peepi` (`main.ntup.sigcm2`) | pion runs only (P9): `ntup%sigcm2` is not a column of the pion ntuple (results_write.f:131-186 writes it only for kaons), and the cross section (`ntup%sigcm`, `main%sigcc`) is the same in both modes. L1 peepi table: 3139 of 5047 rows have W ≥ 2 GeV, 1137 of them with a non-zero stale value. | `fortran_diff_D038_F10_sigcm2_stale` |
| D-038 F1 Δ++ STOP | dbase.f:823-826 | `if(doing_hydpi .and. which_pion.eq.3) stop 'Can''t do Delta+ production from neutron for pi0'` does not test `doing_pizero`: charged H(e,e'π−)Δ++ (set up at dbase.f:361-364, weighted 0.55 at event.f:1479-1484; infiles/fpi12_deltapp.inp) can never run | the STOP only with `doing_pizero` (SetupBuilder); charged π−Δ++ on H runs | H(e,e'π−)Δ++ inputs only (none of the validation runs). The π0Δ0-from-the-neutron case still STOPs in both modes. | `fortran_diff_D038_F1_delta_pp_stop` (incl. a 50-event corrected run) |
| D-038 F2 π0 from n on H | dbase.f:141-142 | π− from hydrogen is stopped, π0 "from the neutron" (`which_pion = 1`, `doing_pizero`) is accepted on hydrogen and runs with Mtar_struck = Mn on a proton target | rejected at start-up (`std::invalid_argument`, exit status 2, "hydrogen has no neutron") | such inputs only (L1 cfg 15 of P9.F1). | `fortran_diff_D038_F2_pi0_neutron_on_h` |
| D-038 F5 CENTRAL.sigcc | simc.f:1236-1290 | calculate_central calls only complete_recon_ev(vertex0): the central peepi sees main0%epsilon = t = phi_pq = 0 (σ_L ∝ \|t\| = 0, ε = 0), so the `.hist` CENTRAL.sigcc is not the cross section at the central kinematics | main0 gets the central W, ε, θpq, t, tmin, Q², φpq (the complete_ev formulas on vertex0; `IReaction::prepareCentral`, `PionProduction`) | `.hist` CENTRAL.sigcc of pion runs only (not used by the weights or the normalisation). K1 π+ H: 0.839106E-06 → 0.925379E-05 (×11.03). §9 lists every L2 run. | `fortran_diff_D038_F5_central_sigcc` |
| D-038 F15 .hist which_pion 11 | simc.f:750-755 | the final-state line is printed for which_pion 0/1/10/2/3 but not for 11 (coherent π−) | `Final State is A + pi` also for 11 | `.hist` text of coherent π− runs. | `fortran_diff_D038_F15_hist_which_pion_11` |
| D-037 F1 SHMS coll slits | shms/mc_shms_coll.f:79, 86, 93 | `shmsSTOP_slit_hor/vert/oct` incremented for each of the 20 steps outside the opening (the SHMS copy of D-021) | at most once per ray and counter, as the HMS (P6.4; the P7.5 correction covered only mc_hms_coll) | `.geni` SHMS slit counters with `using_SHMScoll` and pions/muons. Capture p6f2_pi_shms_cap (3912 SHMS pions): 5646/20/5833 → 386/1/391. Transport, RNG and acceptance unchanged. | `fortran_diff_D037_F1_shms_coll_slit_counters` |
| D-037 F2 coll stops | simc.f:462-532 (hSTOP_coll mc_hms.f:218, sSTOP_coll mc_shms.f:518) | the collimator stops (absorbed pions) are written to no output: the `.geni` stop counts do not add up to the trials | `.geni` line `Collimator stops (pi/mu)` after the slit line of HMS/SHMS whenever `using_HMScoll`/`using_SHMScoll` is on (or the count is non-zero); the summary JSON has always had `coll` | `.geni` text only. Captures: 1253 of 3932 HMS pions (p6f2_pi_hms_cap), 258 of 3912 SHMS pions (p6f2_pi_shms_cap) become visible. The corrected slit counters (D-021) count rays that touched the collimator, not stops, so even then the lines add up only without them. | `fortran_diff_D037_F2_coll_stops_geni` |
| D-037 F4 ranged out | hms/mc_hms_coll.f:121-124, shms/mc_shms_coll.f:123-126 (with enerloss_new.f:75-77) | `if (epart.lt.sqrt(m2)) goto 500` is unreachable: enerloss_new clamps the loss to (E−m)−1e-7, so a particle that ranges out continues with p ≈ 4.6e-3 MeV/c (path up to 2e10 cm) | a clamped loss (`SimcMaterial::lastEnergyLossClamped`) stops the particle in the collimator (hSTOP_coll/sSTOP_coll) | only particles that range out in the tungsten (slow muons from kaon decays). 3 of 1602 rows of mc_hms_decay (720, 891, 990: Q1_in after 4.7e9-2.1e10 cm → coll); 0 calls in the two pion captures. | `fortran_diff_D037_F4_coll_ranged_out` |
| D-037 F7 decay drift | shared/project.f:113 | `tmpdrift = z_drift - z_decay/sqrt(1+dxdzs**2+dydzs**2)` with the slopes AFTER the decay: the track ends at z_drift + z_decay(1/√(1+s_old²) − 1/√(1+s_new²)) | the z already covered uses the parent slopes; the track ends at z_drift | decays in `project` (doing_decay, pions/kaons): positions, pathlen and downstream apertures of the decay product. project_decay table: max \|Δx\|,\|Δy\| 17 cm, \|Δpath\| 26 cm (synthetic edge rows), 0.016 cm / 0.28 cm in the real capture steps. Captures (replay of every recorded call): HMS 141 of 3932 pion calls change pathlen (max 0.28 cm), 43 accepted events change their recon (max 0.0014 % δ, 2.4e-5 xptar, 0.0013 cm ytar); SHMS 104 of 3912 (max 0.13 cm), 28 accepted (max 0.010 % δ); no acceptance decision changed. | `fortran_diff_D037_F7_project_decay_drift` |
| D-037 F9 decay non-π/K | shared/project.f:87-91, transp.f:160-165, 253-258; dbase.f:1147 | the decay of a particle that is neither a pion nor a kaon executes a bare `stop` (exit 0, no output); doing_decay is not tied to the reaction, so e.g. H(e,e'p) with doing_decay = 1 stops at the first proton decay | the set-up rejects doing_decay with a hadron that is not π/K (`std::invalid_argument`, exit status 2, "… is neither a pion nor a kaon …"); project/transp never check a non-π/K particle for decay (no random number). Compat: `simc::FortranStop`, exit 3 | configurations with doing_decay and a proton/ρ/π0 hadron only; none of the reference inputs. | `fortran_diff_D037_F9_decay_non_pi_k`, `cli_decay_non_pik` |
| D-043 F1 fint ABI | physics_kaon.f:262, 375, 492 (`real*4 fint`) vs cern/fint.f:10 (implicit REAL = REAL(8) under -fdefault-real-8, returned as binary128 with -ff2c) | the REAL*4 callers read the low 32 bits of the binary128 result: every fint value is +0.0, so eekeek/eekeeks give sigma_eep = 0 (ntuple column `saghai` = 0 for every event) and phspwght (dead code) 10000 | the value fint computes, rounded to REAL*4 as the declaration intends (`kaon::SaghaiOptions::fint_value`) | ntuple column `saghai` (ntup%sigcm1) of every kaon run; the weight uses sig_factorized and is unchanged. L1: eekeek/eekeeks/phspwght/peek `_fintabi` tables bitwise (1151 of 1244 eekeek rows non-zero) | `fortran_diff_D043_kaon_model`, `l1_*_fintabi` |
| D-043 F2 sigma0 file | dbase.f:662-677 with saghai_sigma0.dat | line 1 of the file is a truncated fragment, so every 4-line record is read one line off: zsrff1..zsiff3 = the 5-value kinematics line, zsrff4..zsiff6 = the true ff1..ff3, the true ff4..ff6 are lost | the fragment is skipped and each record is read from its own lines (`kaon::SaghaiTables::load(dir, true)`); the 18 stray records at the start of the file (s = 2.898 in the slots of s = 2.851, Q2 = 0) are kept where the file puts them and reported (class b below) | `saghai` column of Sigma0/Sigma-/coherent runs only (masked by F1 in the Fortran). 22800 of 22800 corrected ff1..3 values equal the compat ff4..6 | `fortran_diff_D043_kaon_model` |
| D-043 F3 s grid | physics_kaon.f:288-294 | eekeek interpolates on s nodes 2.6 + 0.3(i-1): nodes 9/10 = 5.0/5.3 GeV², but the 9th and 10th blocks read from saghai_proton.dat are s = 5.1/5.4 | nodes 9/10 = 5.1/5.4 (`SaghaiOptions::proton_grid_fixed`) | `saghai` column for s > 4.7 GeV² (W > 2.17 GeV) of K+Lambda runs: 465 of 1244 eekeek table rows, none with s ≤ 4.7 | `fortran_diff_D043_kaon_model` |
| D-043 F5 survival p | physics_kaon.f:160-162 | the survival probability of a kaon whose decay is not simulated uses beta*gamma of the CENTRAL momentum spec%p%P; p_kaon = spec%p%P(1+delta/100) is computed and unused (against the comment l.145-146) | beta*gamma of p_kaon (`PeekConfig::survival_p_kaon`) | EVENT WEIGHTS (and the `decdist` column = survival probability) of decay-OFF kaon runs; §10 has the L2 effect | `fortran_diff_D043_kaon_model` (1407 of 1428 decay-OFF peek rows change) |
| D-043 F6 zaero | physics_kaon.f:149-159 | no zaero branch for hadron_arm 5 (SHMS): the SAVEd value of an earlier call is used (0 in a SIMC run, −82.8/−183 after an arm 2/3/4 call); the path correction `zaero*(1+dx²+dy²)` lacks the square root | zaero = 0 for every arm without an aerogel distance, `zaero*sqrt(1+dx²+dy²)` (`PeekConfig::survival_path_fixed`) | none for HMS/SHMS runs (zaero = 0 in both modes); 60 order-dependent SHMS rows of the peek table | `fortran_diff_D043_kaon_model` |
| D-048 F1 Delta D/A STOP | dbase.f:779-788 (D-018) | the (e,e'p)π announcement tests `doing_deutpi` / `doing_hepi` (pion flags, never set for doing_delta): every D or A Delta run STOPs with `I don't have ANY idea what (e,e'p)pi we're doing!!!` (exit 0) | tests `doing_deutdelta` / `doing_hedelta`: D(e,e'p)π and A(e,e'p)π run (`SetupBuilder`) | D/A Delta inputs only. The D/A branches of event.f (generate l.336-371, complete_ev l.646-656, 949-971) and init.f (l.349-358) then run for the first time: **no Fortran reference exists for them**; they are validated by self-consistency only (energy-momentum balance of every event, `fortran_diff_D048_F1_delta_DA_stop`: D, He-3, C-12). | `fortran_diff_D048_F1_delta_DA_stop` |
| D-048 F1 Delta D/A data | dbase.f:563-565, 594 | the momentum distribution (deut.dat, he3.dat, ...) and the Benhar spectral function (generate_em, `Emval(numEm)`) are read for `doing_deutpi/hepi/deutkaon/hekaon` only, not for `doing_deutdelta/hedelta`, although generate and limits_init use them for the Delta (unreachable in the Fortran because of the STOP above) | the delta flags load the same files as the pion (corrected mode only) | D/A Delta runs only (they exist in corrected mode only) | same |
| D-048 F2 Delta tail 1 | radc.f:288-293 | the incident-electron tail limits are set for doing_heavy, hyd_elast, deuterium and `doing_pion .or. doing_kaon .or. doing_rho .or. doing_semi`, not for doing_delta: Egamma_min(1) = Egamma_max(1) = 0, basicrad returns weight 0, and EVERY try with ntail = 1 fails (frac(1) = 0.43-0.46 of the tries at the P11.F1 kinematics) | doing_delta uses CASE 4 (Egamma_min(1) = 0, Egamma_max(1) = gen%sumEgen%max − vertex%e%E): init.f:455-463 sets gen%sumEgen for doing_delta with the same "generated ELECTRON energy" branch as for the pion (`RadcRun::delta_tail1_case4`) | EVENT RATES of rad Delta runs: Ntried and wtcontr. §11 has the L2 effect. With tail 1 radiating, generate_rad calls complete_ev again (radc.f:322), which draws a new root choice (class b below). | `fortran_diff_D048_F2_delta_tail1_egamma` |
| D-048 F3 sig_blok args | physics_delta.f:99-100 | `sig_blok(thetacm,phicm,main%t/1.d6,vertex%q2/1.d6,s/1.d6,...)`: under -fdefault-real-8 `1.d6` is REAL(16), so binary128 temporaries are passed to the REAL*8 dummies and sig_blok reads their low 8 bytes: ntup%sigcm (column `sigcm`) is meaningless (0 in ~65 %, negative in ~19 %, NaN/Inf in ~0.4 % of the events) | the doubles t, Q², s [GeV²] are passed (`delta::DeltaConfig::sigblok_binary128_quirk = false`): `sigcm` = the Blok π fit at the proton's CM angles (still not a Δ model, see F4) | ntuple column `sigcm` of Delta runs only; the weight is unchanged (F4). peedelta table: 3069 of 3077 rows change, 14 non-finite compat values become finite | `fortran_diff_D048_F3_delta_sigblok_args`, `l1_peedelta` |
| D-048 F6 Delta warning | dbase.f:184-186 | quote typo: `'WARNING: Delta cross section model",` + continuation `" only set up for proton target!'` is ONE string, printed with `",`, 79 blanks and `"` | the sentence without the quote debris | log text of D/A Delta runs | `fortran_diff_D048_F6_delta_warning_text` |
| D-048 F9 CENTRAL.sigcc | simc.f:1236-1290 (as D-038 F5) | the central peedelta sees main0%epsilon = t = phi_pq = 0 | main0 gets the central kinematics (`PionProduction::prepareCentral`, inherited by `DeltaProduction`) | `.hist` CENTRAL.sigcc of Delta runs only: D1 rad 0.227051E+01 → 0.167557E+02 (×7.38) | `fortran_diff_D048_F9_central_sigcc_delta` |
| D-051 F1 sum_sq = 0 | semi_physics.f:508-513 | dsigdz = (...)/sum_sq/zhad without the `sum_sq > 0` guard of the old code (l.378-400, commented out): when every PDF is 0 (x_Bj > 1 clamped to 1, Q < Alambda) the cross section and the event weight are 0/0 = NaN | dsigdz = 0 as the old guard (`PeepiXConfig::sum_sq_guard`) | weight NaN → 0 for those events (the Fortran ΣWeight/.hist sums of such runs are NaN-contaminated if one occurs); 34 of the 3338 P12.F1 peepix rows | `fortran_diff_D051_sidis` (34 NaN rows → 0, no other row), `sidis_sanity` |
| D-051 F2 negative Fermi x_Bj | semi_physics.f:247-248 → Ctq5Pdf.f:80-83 | with do_fermi, x_Bj = Q²/2/(efer ν − pfer p̂·q) is negative for large pfer along q (deut.dat up to 1190 MeV/c): Ctq5Pdf STOPs the run (exit 0, no .hist); a D run with do_fermi = 1 cannot reach 1e5 events (P12.F1: stop at try 27848) | x_Bj ≤ 0 returns 0 before the PDF calls (same early return as the threshold, `PeepiXConfig::nonpositive_xbj_zero`); compat throws `FortranStop` | D/A runs with do_fermi = 1 complete; the affected events get weight 0 | `sidis_sanity` (compat STOP, corrected 0); P12.F1 stops.txt |
| D-051 F4 π0 FF/mass | semi_physics.f:490-495, 227 | the doing_pizero block ("just get averge of pi+ and pi-") assigns `d1 = (yf + yu) / 2.` twice and never `db`: the d̄ fragmentation function keeps `db = u1` of l.487 (the π+ fit yf for doing_hplus, the π− fit yu otherwise); and mhad = Mpi (charged) at l.227 while Mpi_gev (l.149) and dbase.f Mh use Mpi0 | db = (yf+yu)/2 like u1, d1, ub, and mhad = Mpi0 (`sidis::PeepiXConfig::pizero_fixed`) | the d̄ term of the π0 SIDIS cross section; π0 SIDIS is refused until the calorimeter arm (P14), so no run is affected yet | `fortran_diff_D051_sidis` (16 π0 rows of the P12.F1 peepix table: all change, no other row does), `sidis_sanity` |
| D-051 F6 survival p | semi_physics.f:622-625 (as D-043 F5) | the survival probability of a hadron whose decay is not simulated uses β·γ of the CENTRAL momentum spec%p%P; p_kaon = spec%p%P(1+δ/100) is computed and unused (against the comment l.607-608) | β·γ of p_kaon (`PeepiXConfig::survival_p_hadron`) | `decdist` of every decay-OFF SIDIS run; the WEIGHT of decay-OFF K± runs (event.f:1560 applies it for doing_semika only) | `fortran_diff_D051_sidis` (3268 of the 3281 decay-OFF peepix rows change, no other row), `sidis_sanity` |
| D-051 F6 zaero | semi_physics.f:611-621 (as D-043 F6) | no zaero branch for hadron_arm 5 (SHMS): the SAVEd value of an earlier call is used (0 in a SIMC run); the path correction `zaero*(1+dx²+dy²)` lacks the square root (HMS: zaero = −331.491 cm, so the path is ≈ 0.4 cm too short for dx = 0.05) | zaero = 0 for every arm without an aerogel distance, `zaero*sqrt(1+dx²+dy²)` (`PeepiXConfig::survival_path_fixed`) | decay-OFF SIDIS runs with the HMS as hadron arm (sqrt); none for SHMS runs (0 in both modes) | `fortran_diff_D051_sidis`, `sidis_sanity` |
| D-051 F3 pfermi NaN | results_write.f:208-209 | the SIDIS ntuple column `pfermi` = pfer/1000·|d|/d with d = p̂_fermi·q̂, without the `if (dummy.eq.0) dummy=1.e-20` of the pion block (l.143): 0/0 = NaN for every event with pfer = 0 (every H and A > 2 event) | the pion-block guard: pfermi = 0 (`FillSemiNtuple`) | ntuple column `pfermi` of H and A > 2 SIDIS runs (NaN → 0); nothing else | `fortran_diff_D051_F3_F9_semi_pfermi_pmper` (300 H events: compat NaN, corrected 0) |
| D-051 F9 PmPer | event.f:922 before 952-954 (as D-033 F5) | the semi complete_ev computes vertex%PmPer from the previous try's vertex%Pm, assigning Pm = Pmiss only afterwards | PmPer from this event's Pm (`SidisProduction::completeEv`) | vertex PmPer only (not an ntuple column of the SIDIS layout) | `fortran_diff_D051_F3_F9_semi_pfermi_pmper` (compat 0 / corrected 300 of 300 events use the current Pm) |
| D-051 F12 π0 title | simc.f:679-700 | the .hist header has no π0 line: a π0 SIDIS run is titled π+ or π− by doing_hplus (dbase.f announces (e,e'pi0)X) | title (e,e'pi0)X | .hist header (π0 SIDIS itself is refused until the calorimeter arm, P14) | `fortran_diff_D051_F12_F14_semi_hist_titles` |
| D-051 F14 A ≥ 3 kaon | simc.f:701-717 | the .hist header has no A ≥ 3 branch for kaon SIDIS: `STOP I don't have ANY idea what A(e,e'k)X we're doing!!!` after the whole event loop (exit 0, .hist truncated to 102 bytes; .gen/.geni written) | **run**: the A(e,e'k±)X title is written and the report completes. Choice: the run itself is fully defined (dbase.f:768-773 announces it, peepiX handles A > 2 with slacemcfit/Rhad_global and fDSS, the generation has no A-specific kaon branch); only the title branch is missing, so refusing the input would remove a working configuration | C-12 K± SIDIS runs complete with a full .hist | `fortran_diff_D051_F12_F14_semi_hist_titles` (compat STOP, corrected .hist + a 100-event run) |
| D-043 F13 CENTRAL.sigcc | simc.f:1236-1290 (as D-038 F5) | the central peeK sees main0%epsilon = t = theta_pq = phi_pq = 0 | main0 gets the central W, ε, θpq, t, tmin, Q², φpq (`PionProduction::prepareCentral`, inherited by `KaonProduction`) | `.hist` CENTRAL.sigcc of kaon runs only: K1 K+Λ H 0.143536E-06 → 0.704013E-06 (×4.90); §10 lists every L2 run | `fortran_diff_D043_F13_central_sigcc_kaon` |
| D-057 F6 rho efer | rho_physics.f:88-93 | peerho recomputes the COMMON efer on shell, sqrt(pfer²+Mtar_struck²), and writes it back to /pfermi/: the off-shell branch tests the PION flags doing_deutpi/doing_hepi, never set in a rho run (1001/1001 D capture calls: efer 936.047 → 938.272) | peerho uses the COMMON efer of the generation and leaves it unchanged (`rho::PeerhoConfig::efer_from_common`); which formula is right for A > 1 stays on the expert list (D-057) | D/³He rho weights (s, flux, CM boost); none for H (pfer = 0, efer = Mp in both). D: ×0.828 on a pfer = 200 MeV/c event; in a run pfer = 0 (F4), so only efer = Md − Mn vs Mp | `fortran_diff_D057_rho_kernels`, `l1_peerho` (compat) |
| D-057 F3 rho Mh2_final | generate_rho.f:80 | `Mh2_final = Mh` (a mass, not a mass squared), 5695/5695 calls | Mh2_final = Mh2 (`rho::GenerateRhoConfig::mh2_final_fixed`) | none observable (simc.f:174 resets it, rho_decay sets Mpi2 on every success) | `fortran_diff_D057_rho_kernels` |
| D-057 F1 rho stale hadron angle | event.f:280-283, 325, 1034; generate_rho.f:123-124 | generate draws no hadron slopes for doing_rho, so physics_angles converts the previous try's rho lab theta/phi (stored in vertex%p%xptar/yptar) and trip_thru_target(3) uses that angle for the rho energy loss and teff(3) (θ used: mean 1.06 rad vs the rho's 0.42 rad; ΔEloss(3) −6…+2.4 MeV, teff(3) −24…+58 %) | trip_thru_target(3) uses the rho's own lab angles (`RhoProduction`, corrected mode) | Eloss(3)/teff(3) of every event: the decay pion's energy-loss correction and multiple scattering in montecarlo | `fortran_diff_D057_F1_F2_rho` (3975 capture tries, teff(3) changes in all) |
| D-057 F2 rho_decay after failure | event.f:420-422 | rho_decay is called for every try that reaches l.420, also after a failed generate_rad: it decays the PREVIOUS try's orig (3.4 % of the H tries), draws 3–41 uniforms and may set Mh = Mpi | rho_decay only after a successful generation. No weight changes; the RANDOM STREAM differs from the Fortran's from the first such try on (those 1+2k uniforms are not drawn), and Mh/Mh2/ntup%rhotheta/orig keep their values on those tries | the random stream of a corrected run (no Fortran-identical replay in that mode); no physics quantity | `fortran_diff_D057_F1_F2_rho` (94 capture tries: compat draws the decay uniforms, corrected not; outcome identical) |
| D-057 F9 rho yptar sign | rho_decay.f:119-124 | for a beam-left arm (2/4/5) and a decay pion on the far side of the beam (pyf ≤ 0), th_inp = +(α + θ_spec): yptar gets the wrong sign, and orig%p%theta/phi follow | th_inp = −(α + θ_spec) (`rho::RhoDecayConfig::yptar_sign_fixed`) | such pions have |yptar| ≥ tan θ_SHMS, outside the acceptance (741 of 2782 decays in the first 3000 H tries): only failed tries and the .geni of such tries; with the sign fixed some far-side tries no longer fail the 'going backwards' test (they still miss the acceptance) | `fortran_diff_D057_rho_kernels` |
| D-057 F20 rho NaN | rho_physics.f:388-389 | the 'crazy sigma' guard (> 1e10 or < 0 → 0) lets NaN through (s < 0 edge row) | NaN → 0 (`rho::PeerhoConfig::nan_guard`) | edge cases only (no NaN in the L2 runs) | `fortran_diff_D057_rho_kernels`, `l1_peerho` |
| D-057 F21 rho NaN pion | rho_decay.f:111-112 | acos(pz/p/cos θ_oop) with an argument above 1 by rounding gives th_inp = NaN, yptar = NaN; mc_shms and pass_cuts let it through and the event contributes with every recon pion column NaN (1 per 1e5 events: rho_h_k1_rad_A record 57258, rho_he3_k1_rad_A 54227) | the argument is clamped to [−1, 1] (`rho::RhoDecayConfig::acos_clamp`) | about 1 event per 1e5 (its weight ~6e-8 of the sum); NaN tree entries disappear | `fortran_diff_D057_rho_kernels` |
| D-057 F15 rho .hist title | simc.f:720-725 | the .hist header says 'I am not set up for anything else yet!' for D and ³He rho runs, which run | ` ****--------  D(e,e'rho)  --------****` / `A(e,e'rho)` (`RunReport`) | .hist header of D/³He rho runs | `fortran_diff_D057_F15_F18_rho` |
| D-057 F18 rho P-arm report | simc.f:259-261, 283-285, 316-317, 324-325; event.f:22-24 | the .geni/.gen P-arm histograms, the .hist 'Hadron arm' resolution sums and the 'Limiting VERTEX P arm' ranges use vertex%p, which for doing_rho is the RHO with its lab theta/phi in the slope fields (rho_h_k1_rad: Ave.Error yptar −1580.9 mr, limits 589…2335 mr) | the decay pion orig%p (`GeneratedHadronArm`, corrected mode) | .geni, .gen (and gen.root) P-arm distributions, .hist resolution lines and limits; no weight | `fortran_diff_D057_F15_F18_rho` |

### Class a-np: clear bugs in paths that are not ported yet (fix when porting)

| id | Fortran | bug | C++ now | to do |
|---|---|---|---|---|
| D-006 tgt field | simc.f:139 | `ang_targ_parm = +1*sin(phi)*(targ_Bangle-theta)` (the e-arm line l.125 has `theta-targ_Bangle`) | `using_tgt_field` is NotPorted | when the target field is ported |
| D-011 HRS .geni | simc.f:489, 501 | "Slit hor/vert" prints `slit_vert` first | the C++ writes no HRS block yet; the Fortran-.geni parser maps by meaning | P14.2 |
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
| D-038 F4 coherent model | physics_pion.f:807-854 (exclfit), event.f:1496-1500 | coherent pion production (which_pion 10/11) uses the nucleon fit with Mtar_struck = M_A (W = √s photon-nucleus 4–13 GeV); exclfit keeps `mtar_gev = 0.938` (l.812); tgtweight = Z or N ("OK for coherent???"); sig_param_2021 is negative in 210 of 4564 table rows (no clamping) | is a nucleon model meant for coherent production, and should negative cross sections be clamped? | coherent pion weights (P9) |
| D-038 F3 coherent Mrec | dbase.f:283-288, 380-388 | for coherent production the input `targ%mrec_amu` is first replaced by the A−1 guess M − Mp, then by M − M_N(struck) + M_N(recoil); the input value is never used, even when it is the correct final nucleus (two replacement messages) | should the input final-nucleus mass be honoured? | coherent pion runs: Mrec_struck (He-3 coh: 2809.68 MeV) |
| D-038 F9 Δ re-generation | event.f:620-631, radc.f:324 | complete_ev overwrites the COMMON targ%Mrec_struck with a new Breit-Wigner Δ mass on every call; radiation tail 1 calls complete_ev again, so the event gets a second Δ mass and (π0) a second photon decay; the `.hist` prints the last Δ mass of the run (simc.f:844; SIMC++ too, `IReaction::finalTargetState`) | intended re-generation? | π Δ final states (which_pion 2/3) with radiation |
| D-038 F12 sig_multipole | physics_pion.f:197-398 (dead code: the call in peepi is commented out, l.157-163) | COMPLEX*8 (single-precision) amplitudes although the build is -fdefault-real-8; at exact (Q², ν) nodes or outside the table the amplitudes of the PREVIOUS call are used; the argument `nu` is overwritten (≥ 500 → 499.9); 197.327/mh fixed at the first call. SIMC++ reproduces all of it (`pion::SigMultipole`, `std::complex<float>` storage, Fortran kind rules per operation, SAVEd members) | single precision and the previous-call amplitudes intended? | none (not called by SIMC); L1 sig_multipole 612/612 rows bitwise |
| D-037 F3 coll scattering | hms/mc_hms_coll.f:115-118, shms/mc_shms_coll.f:117-120 | musc and the Landau-sampled enerloss_new (2 + 1 gauss1) of every collimator step are applied whatever ms_flag / mc_smear say | is the collimator scattering intended physics (the purpose of the stepping) or should it follow ms_flag? | pion/muon runs with using_*coll and mc_smear = 0 |
| D-037 F5 pion_coll_absorb | hms/pion_coll_absorb.f:16, 30-56 | own mpi = 139.56995 (constants.inc 139.57018); transmission 1 for Tpi ≤ 85 MeV, 0.938 per step just above; unused mate_t; the SHMS uses the HMS material (17 g/cm³, A = 171.57) | intended parameters? extrapolate below 85 MeV? | absorption probability of pions in both collimators |
| D-037 F8 transp decay | shared/transp.f:142-185, 232-276 | a first-half decay is applied at the element entrance (then the map acts on the daughter), a second-half decay at the exit with the parent's β, γ; decdist gets z_decay, pathlen zd+δz | an intended approximation (the COSY map cannot be split)? | decays inside magnets (doing_decay) |
| D-043 F2 data | saghai_sigma0.dat (content) | the file starts with 18 stray records (s = 2.898 GeV², Q² = 0.52/0.75) and has 1 record instead of 19 for s = 2.851, Q² = 0; line 72 has a truncated field (`-00E-02`, read as −0) | which amplitudes belong to s = 2.851, Q² = 0? | `saghai` column near the Σ threshold at Q² < 0.25 GeV² only |
| D-043 F7 decay model | shared/transp.f:149-160, 239-251; project.f:76-85 | K decay: fixed 70 % μν / 30 % ππ, p* = sqrt(Mk²/4 − Mπ²) = 203.6 MeV/c (equal-mass formula with the charged pion; K⁺→π⁺π⁰ is 205.1), the daughters never decay | intended approximation? | decay-ON kaon runs (the Mk of transp.f:158 is corrected, D-006) |
| D-043 F8 Saghai angles | physics_kaon.f:19, 79-85 | peeK's `phi` is never assigned (SAVEd, 0) and eekeek's `theta` argument is main%theta_pq (the lab p–q angle) | intended (comment "ALWAYS CALCULATING FOR PHI=0")? | `saghai` column |
| D-043 F9 coherent kaon | physics_kaon.f:82-86, 215; event.f:1505 | coherent production (which_kaon 10/11/12) uses the nucleon formulas with W of the photon–nucleus system (4–13.5 GeV); with Mrec_struck > 1150 MeV the coherent Λ uses the Σ branches; the cross section is the same for 10/11/12 (only tgtweight Z/N differs, "OK for coherent???") | is a nucleon model meant for coherent production? | coherent kaon weights |
| D-043 F10 Mp | physics_kaon.f:194-197 | sig_factorized uses the proton mass in ν/q_cm also for Σ⁻ (neutron) and coherent production | intended? | Σ⁻ and coherent weights |
| **D-048 F4 Delta weight** | physics_delta.f:130-132, event.f:1512 | **peedelta = 1.0*jacobian*(gtpr*fac): the cross-section line is commented out ("Need new xsec model."), so the Delta weight is phase space × photon flux × the transform_to_cm jacobian, with NO cross section.** The sig_blok value it computes (ntup%sigcm) is the Brauel p(e,e'π+)n fit evaluated with the proton's CM angles and the proton t (π− branch for which_pion 1/3/11): not a pπ0 or Δ model either | should SIMC++ offer a real Delta (e,e'p)π cross section, and which one? (no choice can be derived from the Fortran) | every Delta weight, wtcontr, `Weight`, `siglab` (14.8–351 in delta_h_d1_rad_A) |
| D-048 F5 two roots | event.f:682-684 | for doing_delta one grnd picks the second root of the proton-energy quadratic when grnd() > 0.5; no factor 1/0.5 enters gen_weight or the jacobian, so when both roots are physical each is generated at half rate with full weight (and the first-root check `p%E < 0` comes before the choice) | is a factor 2 missing? | Delta weights where both roots reach the acceptance (D2: 3 of 1000 capture events took the second root; D1: none) |
| D-048 F5b root re-choice | radc.f:322 with event.f:682-684 | with tail 1 radiating (corrected F2), generate_rad calls complete_ev again, which draws a NEW root choice (the radiated event may take the other root than the one the tail limits were computed with), as D-038 F9 for the Δ mass | intended re-generation? | corrected-mode rad Delta runs (tail-1 events) |
| D-048 F7 recoil mass | dbase.f:297-300 | targ%Mrec_struck = Mpi = 139.57018 (charged pion) although for γ*p → p X the missing pion is a π0 (134.9766) | π0 mass intended? | generated Delta kinematics; central recon MM 138.92 MeV |
| D-048 F8 which_pion | event.f:470-491, 612-623; physics_delta.f:99 | complete_ev reads which_pion also for doing_delta: with 2/3 it reads delta_relativistic_bw.inp and replaces the recoil pion by a Breit-Wigner Δ mass on every call (below threshold in the Δ region: every complete_ev fails, an L2 run would never finish); 1/3/11 select the π− branch of sig_blok (only `sigcm` changes) | should which_pion be ignored for doing_delta? | Delta runs with which_pion ≠ 0 |
| D-051 F7 threshold | semi_physics.f:278-284 | the "2π / πK threshold" cut uses mmpi2 = M² + 2Mν(1−z)(1−pt²[GeV²]) — a mass² multiplied by (1 − pt²) — and returns 0 below (M + m_h)² | which missing-mass formula was intended? | weight 0 or not for events at large z / large pt² |
| D-051 F8 semi thresholds | event.f:979-996 | complete_ev rejects a SIDIS try when (Mtar_struck + ν − E_h)² − P_miss² < (Mp + Mpi0)² — for every hadron (also kaons, where the lowest final state is ΛK), with Mtar_struck and no Fermi motion; the following z > 1 test is unreachable (z > 1 already fails the threshold); 71 % of all tries are z > 1 events | which threshold is intended for kaons / with Fermi motion? | generation efficiency and the accepted phase space near z → 1 (kept in both modes, `SidisProduction::completeEv`) |
| D-051 F16 hadron range | init.f:504-506 | with radiation gen%p%E reaches sumEgen%max − edge%e%E%min (hadron δ up to +116 %); 2.9 % of the rad events have orig δ > 80 % and 1612/100000 of them reconstruct inside −10…22 % (SHMS model) | intended generation range? | 0.68 % efficiency (rad) vs 4.8 % (norad); 1.2 % of ΣWeight from δ > 80 % (kept, set-up) |
| D-051 F10 model choices (R_h for D) | semi_physics.f:659-677 | Rhad_global (the nuclear hadron-attenuation fit) returns 1 only for A = 1: deuterium gets N₀ z^α (1−z)^β, 0.955 at z = 0.5 | intended for D? (D-051 F10 also lists the slacemcfit and Mtar_struck = Mp choices) | D(e,e′h)X weights (≈ −4.5 %) |
| D-043 F12 warning | physics_kaon.f:201-206 | the "Kaon C.M. momentum passed to sig_factorized does not agree" warning is printed for the central event of Σ⁰, Σ⁻ and A > 1 quasifree runs (vertex0 places the kaon at spec%p%P whatever the hyperon) | message only | log text |
| D-057 F4 rho nuclei | dbase.f:563, 594; event.f:336-372; sf_lookup.f:195-245 | no momentum distribution (nump = 0) and no ³He spectral function are loaded for doing_deutrho/herho: the Fermi block draws 4 uniforms but pfer = 0 always; ³He generate_em reads its empty table out of bounds and returns Em = 0 (2 uniforms); D: Em = 2.22494, efer = Md − Mn | should D/³He rho include Fermi motion / a spectral function? (a model limitation, kept in both modes; SIMC++ computes the ³He Em = 0 without the out-of-bounds read) | D/³He rho runs are a nucleon at rest (rho_d_k1_rad_A: wtcontr = 2.33 × H, i.e. Z + N with the same kinematics) |
| D-057 F5 rho s | generate_rho.f:72 | s = −Q² + efer² + 2ν·efer lacks the −2q·pfer − pfer² terms while the boost includes them | intended? (latent: pfer = 0 always, F4) | none today |
| D-057 F7 rho t | rho_physics.f:115, 233 | t and tmin of peerho use the nominal Mrho², not the event's Breit-Wigner mass (complete_ev's main%t uses it) | model defined at the pole mass on purpose? | the exp(−b·|t − tmin|) factor of every event |
| D-057 F8 rho decay frame | rho_decay.f:57-77 | the decay angles (cos θ* from sin²θ* + 2εR cos²θ*, φ* flat) are thrown w.r.t. the LAB z/x axes in the rho rest frame, not the helicity frame; ntup%rhotheta is that angle. (The header comment says εR, the code 2εR = SCHC; the envelope is valid) | intended approximation? | the detected pion's angular distribution (the 2εR cos² part) |
| D-057 F10 rho xptar | rho_decay.f:110, 142 | xptar = tan(asin(px/p)) instead of px/pz′ (up to 4.3e-4 inside ±0.1; orig%p%theta/phi differ from orig%up by ≤ 0.42 mr) | approximation accepted? | decay-pion xptar |
| D-057 F11 rho tail 3 | radc.f:386-387; init.f:231-232 | tail-3 radiation of the NEUTRAL rho, with limits from the hadron edges computed with Mrho: 70 % of tail-3 tries fail at the limits (91682 / 131589 H tries) | should the rho radiate, and with which limits? | rad efficiency (rad/norad wtcontr 0.854); accepted events have Egamma(3) up to 1.5 GeV |
| D-057 F12 rho tail 1 | radc.f:324 | a tail-1 photon re-calls complete_ev, which generates a NEW rho (new BW mass and CM angles; as D-038 F9) | intended? | rad events of tail 1 (393 of 1000 capture events) |
| P14.1 SOS 7.62 cm | mc_sos.f:337-340, 349-360 | the last 7.62 cm to the new exit window are a straight line without `project`: no pathlength (FP path −7.62 cm), no decay check over that segment; the code comment says "WE REALLY NEED TO FIX THIS UP" | include the 7.62 cm in pathlen / decay? | SOS FP path (TOF), decay ON runs |
| P14.1 SOS recon fry | simc.f:1349-1355, 1472, 1738 | mc_sos (and mc_hrs*) receive the raster `fry` (−rastery, 0 without correct_raster) as the recon vertical position; mc_hms/mc_shms receive xtar_init (the true x at z = 0) | intended difference (SOS/HRS matrices fitted with the raster), or a leftover? | SOS ytar/delta recon with an extended beam spot |

## 5. Class d: faithful Fortran behaviour, not a bug (kept, `FORTRAN-QUIRK`)

| id | Fortran | behaviour | why kept |
|---|---|---|---|
| D-018 CTP float | CTP/thTestParse.c:392-393 | every real .inp literal is stored as a 32-bit float | defines the run inputs |
| D-010 literals | -fdefault-real-8; physics_proton.f | double literals, `pow` for real exponents, x+y+z dot products | bitwise numerics |
| D-024/D-018 powers | gfortran | `x**2` → x·x; `x**n` (integer n) → `__powidf2` (`optics::powi`) | bitwise numerics |
| D-015 F1 / D-021 lfit | cern/lfit.f, mc_*_hut.f | focal plane real*4, accumulators real*8, a hit at exactly 0.0 ignored (KEY = 0) | CERNLIB behaviour |
| D-015 F3 closure | COSY maps | forward→recon closure not exact (HMS xptar RMS 1 mr) | optics |
| D-015 fry | mc_hms_recon.f:118, mc_shms_recon.f:136 | fry = 0 replaced by 1e-30 | avoids 0**0 |
| P14.1 SOS recon guard | mc_sos_recon.f:120-122 | every hut value (xfp, xpfp, yfp, ypfp, fry) with \|v\| ≤ 1e-30 is replaced by 1e-30 (HMS/SHMS: only fry) | `CosyEngine::recon(2, ...)` |
| P14.1 SOS hut start | mc_sos_hut.f:283-284, mc_sos.f:369 | the hut starts at zinit = −3.206267 cm, after the exit foil at −3.22 cm: the first drift is negative and clamped to 0.001 cm | faithful (`McSosHut.cpp`) |
| P14.1 SOS cal | mc_sos_hut.f:519-525, 535-539 | the calorimeter aperture and fiducial cuts are commented out: sSTOP_cal is never incremented | faithful |
| P14.1 SOS BM01 | mc_sos.f:234-235 | the BM01 entrance test is asymmetric in y (`yt > w_bm01` or `-yt > w_bm01-0.05*2.54`) | faithful (deliberate-looking 0.05 inch) |
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
| D-038 F11 krel | event.f:961-977 | pion: `ntup%krel` = the SAVEd local krel, set only for D/A quasifree production; 0 for H and coherent production (not an ntuple column of the pion layout) | not an output |
| D-038 F16 STOP status | dbase.f (all STOPs) | every STOP exits with status 0 | SIMC++ exits 3 (`FortranStop`) or 2 (rejected input), class dev D-005 |
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
| D-038 F13 | physics_pion.f:467-516 (sig_param04, not called) | pole at \|t\| = 0.02 GeV² (Inf/NaN), `alog` on a real*8 argument (= double log with -fdefault-real-8) | the routine's formula; dead code |
| D-038 F14 | physics_pion.f:659-669 (sigmaid) | nearest (lower) bin in Q², W and cos θ*, no interpolation; only θ rows 1-6 of the 23 stored per (Q², W) are reachable; σ_T divided by max(0.2, Q²) and by (W−1.132)/0.1 above 1.232 GeV | the coded model (MAID2007 tables of P. Bosted) |
| D-038 pi0 decay | pizero_decay.f:1-94 | 2 grnd per call (φ* then cos θ*); on `Mh/2 < Mgamma` (impossible, Mgamma = 0) the photons are not written | `pion::pizeroDecay` returns `success = false` and the caller keeps its values |
| P6.4 1.d-27 | hms/pion_coll_absorb.f:76, 82 | `sigA*1.d-27`: a DOUBLE PRECISION literal promoted to real*16 by -fdefault-real-8, so the product is formed in quad precision (objdump: __multf3) with the quad value of 1e-27; SIMC++ does the same (`__float128`); the double product differs by 1 ulp in 7 of 788 table rows | bitwise numerics |
| D-037 F10 save_coll | mc_hms.f:216, mc_shms.f:517 | `save_Hcollsuccess/save_Scollsuccess = 1.5` on a collimator stop, never reset or read (dead COMMON /save_coll/); mirrored as `McHms/McShms::saveCollSuccess()` for the L1 tables only | no effect |
| D-037 F12 thick_temp | mc_hms_coll.f:79-96, mc_shms_coll.f:81-98 | `thick_temp` is SAVEd and set only on a hit step (`thick = thick_temp` is read on every step, used only on hit steps); kept as a member | no effect |
| D-037 F6 | shared/transp.f:158 | (already class a, D-006 kaon decay) now covered bitwise by transp_decay_hms/shms (55/103 rows of the Mk case) | — |
| D-037 F13/F14 | gauss1.f, call_ranlux.f; CTP | real*16 return ABI of grnd/gauss1 (harness only); CTP float inputs (D-018) | — |
| D-043 F4 binary128 | physics_kaon.o | DOUBLE COMPLEX z1a..z8a, ur, ui, dble() and the literal 1.d+06 are kind 16 under -fdefault-real-8: the Saghai amplitude sums run in binary128 (`__float128`, glibc `hypotf128` for abs), px(2) = q22/1.d+06 is rounded quad → float; cern/fint.f computes in REAL(8) (implicit typing), not in single precision as CERNLIB | bitwise numerics (L1 eekeek/eekeeks `_fintabi` bitwise) |
| D-043 F14 krel | event.f:961-977 | as D-038 F11: ntup%krel never set on the H/coherent kaon branch | not an output |
| D-043 F16 STOP | dbase.f, physics_kaon.f:539 | every STOP exits with status 0 | SIMC++ exits 3 (`FortranStop`) |
| P11.1 tfcos | physics_delta.f:74 | `tfcos-1. .lt. 1.d-8` compares in binary128 (1.d-8 is a REAL(16) literal, __lttf2); SIMC++ does the same (`__float128`, both modes) | bitwise numerics |
| D-048 F9 stale fields | event.f:961-977, physics_delta.f | ntup%krel = the SAVEd local krel, never set for H (0); ntup%sigcm2 is never set by peedelta (not columns of the Delta ntuple); tgtweight = 1 on the Delta branch (event.f:1442) also for D/A | not outputs / the coded weight |
| D-048 F11 STOP | dbase.f:779-788 | the Delta STOP exits with status 0 | SIMC++ exits 3 (`FortranStop`, compat only; corrected: no STOP) |
| P12.1 Q1 early return | semi_physics.f:280-284 | below the mmpi2 threshold peepiX returns 0 before main%davejac, ntup%sigcm and the survival probability (decdist) are set: those ntuple columns keep the previous event's values (ntup%xfermi is already written) | weight 0; faithful stale-value semantics (D-024) |
| D-051 F5 dead inputs | semi_physics.f:142, 525; dbase.f:1096 | the input pt_b_param is read and ignored (`b = pt_b_param` is commented out); the pt² slope is b = 1/(0.120 z² + 0.200) GeV⁻² | kept as the Fortran (D-051 F5); SIMC++ prints a WARNING when pt_b_param ≠ 0 (`sidis::ignoredInputNotes`) |
| D-051 F13 CENTRAL.sigcc | simc.f:1236-1290 | the central peepiX call uses vertex0 with z and p_T² from complete_recon_ev: meaningful (unlike pion D-038 F5); its survival factor uses main0 FP path 0; the log prints 'multiplicity' = σ_had of that call | kept (SidisProduction::prepareCentral is a no-op); the SIDIS central log lines (simc.f:1288-1303) are not printed by SIMC++ |
| D-051 F5 doing_cent | event.f:1523; semi_physics.f:157-210 | the central cross section call peepiX(..., .TRUE.) is commented out: the inputs sigc_flag, sigc_nbin, sigc_kin_min/max/ind are read and never used | the branch is ported (`sidis::peepiX(..., doing_cent = true)`) and never called; SIMC++ prints a WARNING when they are set; main%sigcent / ntuple `sigcent`, `centjac` are always 0 (P12.2) |
| D-051 F5 semi_dilution | event.f:1524; semi_dilution.f | dead code: the call is commented out (`ntup%dilu = 1.0`, not an ntuple column) | ported (`sidis::semiDilution`), not called; unit-tested only |
| P12.1 Q2 sqesub | semi_physics.f:288-296; F1F2IN21_v1.0.f:1309-1401 | the first peepiX call runs SQESUB, which only reads the deuteron smearing tables f1f2tables/f1*-cdbonn.dat into SAVEd arrays used by the deuteron routines; F1F2IN21 is called with A = 1 | no effect on any output; not ported (the Fortran run needs the files, SIMC++ does not) |
| P12.1 Q3 resmod cut | F1F2IN21_v1.0.f:664, 953 | `if((w-mp).LT.wdif(1)) sig = 0.0` with wdif(1) = w−(mp+mpi): never true | no effect (the xpr thresholds do the cut) |
| P12.1 Q4 Ctq5Pdf Q | cteq5/Ctq5Pdf.f:77-80 | Q < Alambda: 'Q out of range' and the CALLER's Q (peepiX Qgev) is set to Alambda for the remaining PDF calls of the event | Q² < 0.051 GeV² only |
| P12.1 Q5 binary128 | fdss/fdss.o | the D literals of fdss.f are REAL(16): (1.D0−X)**4 and **7 are `_gfortran_pow_r16_i4` calls (repeated squaring in binary128), the products fFINT·(1−X)ⁿ·X^p are binary128 rounded once to REAL(8), FAC·(1.D0−D) in fFINT is a binary128 product; Ctq5Pdf.o and semi_physics.o have no binary128 arithmetic (only the exact -ff2c return conversions) | bitwise numerics (`__float128` in src/physics/sidis/Fdss.cpp) |
| P12.1 Q6 π survival | event.f:1560; semi_physics.f:610-627 | with doing_decay = 0 peepiX writes the survival probability to decdist also for pions, but complete_main applies it to the weight for doing_semika only | faithful (the decdist column of decay-OFF π runs is a survival probability, not a distance) |
| P12.1 Q7 xbj clamp | semi_physics.f:255-258 | x_Bj > 1 prints 'XBj is too large!' and uses x_Bj = 1 (every such event) | faithful |
| P12.1 Q8 PartonX x ≤ XV(0) | cteq5/Ctq5Pdf.f:148-156, 190-194 | the clamp Jx = max(Jx, 0) is only done in the first call with x < Xmin (SAVEd First): later an x ≤ XV(0) = 0 (x = 0, or a NaN x) gives Jx = −1 and POLINT reads XV(−1) and UPD(J0) outside the arrays, i.e. Al and (for JFL = 0, the b̄, with Jq = 0) QL(MXQ) of COMMON /CtqPar1/ | reproduced by a contiguous COMMON storage (`sidis::Cteq5Pdf`); x = 0 needs Q² = 0, never reached by peepiX; the P12.F1 x = 0 row is bitwise |
| D-057 F7b rho main%w/t | rho_physics.f:113, 117 | peerho overwrites main%w with √s in GeV and main%t with t in GeV² (complete_ev set them in MeV, MeV²) | faithful: the rho ntuple/report see the GeV values as in the Fortran |
| D-057 F2b rho_decay success | rho_decay.f:1-160 | rho_decay never sets success = .true.: the incoming success is returned on every non-failure path | faithful (`rho::rhoDecay(..., success_in, ...)`); the call after a failed try is corrected (F2, §2) |
| D-057 F13 rho weight | event.f:1016-1022; simc.f:375-391 | no hadron jacobian, genvol = dΩe · 4π · ΔEe (flat CM generation, peerho = gtpr·dσ/dΩcm); the truncated Breit-Wigner sampling (±500 MeV: 90.5 % of the BW) carries no normalisation factor; davejac is computed but not applied | faithful |
| D-057 F14 rho ntuple | results_write.f:187-229 | 59 columns (OUTPUT_FORMAT.md §1f): `t` is the PION t (Mh2 = Mpi2), zi/pt2i/sigcent/xfermi always 0, Mhadron in MeV, pfermi NaN for H (D-051 F3: 0 in the corrected mode) and ±0 for D/³He, jacobian/centjac = davejac/johnjac, phipqi = the rho's φpq; with using_tgt_field `mmnuc` is never filled (not ported) | faithful (pfermi: D-051 F3) |
| D-057 F17 rho CENTRAL.sigcc | simc.f:1236-1290 | CENTRAL.sigcc = peerho for the SHMS central particle treated as a rho of nominal mass with main0%phi_pq = 0 (as D-038 F5) | kept in both modes (`RhoProduction::prepareCentral` does nothing) |
| D-057 F19 rho pass_cuts | simc.f:237 (D-006) | the electron delta max in the hadron test changes Σ Weight of the rho runs by −24 … −34 % | corrected by D-006 (§2) |
| D-057 F22 rho_decay dead code | rho_decay.f:54, 132-136 | dlen/ctaurho are computed and unused; the 'going backwards' test is unreachable after the pzprime test | faithful |
| D-057 diag222 | rho_physics.f:252, 262, 293, 390 | the diagnostics (R < 0, sigt < 0, gtpr ≤ 0, 'crazy sigma_rho') go to unit 222 (fort.222; none written in the L2 runs) | kept in the model (`RhoModel::diagnostics()`, printed with 'rho (fort.222):' when verbose); no fort.222 file |
| D-057 literals | rho_physics.f:280, 384 | 1.d+06 and 1.d3 are REAL(16) under -fdefault-real-8: sig219/1.d+06 and davesig/1.d3 are binary128 divisions | `__float128` divisions (bitwise numerics) |

## 6. Deliberate design deviations (class dev)

| id | Fortran | SIMC++ |
|---|---|---|
| D-002 RNG | RANLUX (call_ranlux.f; skips every 1000th number) | `mt64` (std::mt19937_64); `random_state_file` ignored; validation is statistical. Bitwise replays use `ScriptedRandom` with captured streams. |
| D-005 STOP | `STOP` exits with status 0 | `simc::FortranStop`, exit status 3; rejected inputs exit 2 |
| D-003 input | CTP .inp, `%` names only (obsolete `.` names silently ignored) | .inp (same CTP semantics, unknown names reported like CTP) or JSON; `inp2json` |
| D-032 | arm codes as raw integers, 6 half-supported | one `SpectrometerId` table; 6 rejected (see §2) |
| not ported | HRS/calorimeter arms (π0 with the calorimeter, also π0 SIDIS, rho with an HRS hadron arm), target field | `simc::NotPorted` at start-up (P14.2/P14.3); the SOS (arm 2) is ported (P14.1) |
| D-057 F16 rho A ≥ 4 | dbase.f:209-214, 875-876: no doing_*rho flag for nint(A) ≥ 4, the run STOPs (exit 0) after printing E_bind = −1e10; hadron_arm 6/7/8 STOP in rho_decay | refused at input in both modes (std::invalid_argument; A ≥ 4: SetupBuilder; arms 7/8: the "rho" model; 6 is rejected for every run, D-032) |
| P13.1 rho + other reaction | dbase.f:136-215 takes the first branch of its if-chain (pion, kaon, delta, semi) and leaves doing_rho set, so event.f:420-421 decays that hadron as a rho | doing_rho together with doing_pion/kaon/delta/semi is refused at start-up (std::invalid_argument, both modes); doing_rho with doing_pizero too (results_write.f:230-236 overwrites the rho columns) |
| P13.2 RANLUX | grnd = RANLUX luxury 3 (cern/ranlux.f, call_ranlux.f: 999 of every block of 1000), gauss1.f | besides the default `mt64` (D-002), the IRandom `"rng": {"type": "ranlux"}` reproduces the Fortran stream bit for bit from random_seed (or `"state"`: the 25 integers of *_start_random_state.dat): with `--fortran-compat` a SIMC++ run can then be identical to a Fortran run (both P13.F1 rho captures: every try, tree, .hist/.geni/.gen) |
| P12.1 SIDIS files | Ctq5Pdf.f:SetCtq5 opens cteq5/cteq5m.tbl, fdss.f opens fdss/KANLO.GRID (kaons) from the working directory at the first peepiX call | the "sidis" model reads the copies in `data/sidis` (md5-checked, `l1_sidis_data_files`; JSON `physics.data_dir`) at start-up; a missing file is an error at start-up |
| D-033 F1 | the spectral-function file follows from nint(A) only | also selectable in the JSON (`pipeline.reaction.spectral_function`: type and file); a warning when the Fortran C-12 fallback is used |
| outputs | .hist/.geni/.gen, .bin, .root | the same text files, `.root` (h10, same 46 branches), plus `.summary.json/.txt` and `.gen.root`; no `.bin` |
| P10.1 kaon files | dbase.f:624-679 reads weightc.dat, weightd.dat, saghai_proton.dat, saghai_sigma0.dat from the working directory | the "kaon" model reads the copies in `data/kaon` (md5-checked, `l1_kaon_data_files`; JSON `physics.data_dir`) at start-up |
| P9.1 MAID files | sigmaid reads maid*.dat from the working directory at its first call (W < 2 GeV), sig_multipole MULTIPOL.file | the "pion" model reads the copies in `data/pion` (md5-checked, `l1_pion_data_files`; JSON `physics.data_dir`) when it is initialised: a missing file is an error at start-up, not in the middle of a run |
| P9.1 peepi writes main | peepi(vertex, main) writes main%thetacm, phicm, pcm, davejac, johnjac, wcm and ntup%sigcm* | `IPhysicsModel::sigccMain(vertex, main&)` does the same; `sigcc(vertex, const main&)` evaluates on a copy. COMMON /pfermi/ is `MainEvent::pfermi`, ntup%sigcm1/2 are `NtupVars::sigcm1/2` |

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

## 9. Exclusive pion production: L2 in both modes (task P9.2)

`simc_compare` (default config) of simc++ seed 4242, 200000 successes, against the P9.F1 Fortran runs
`reference/runs/pion_<tag>_A` (and `_B` for h_k1_rad_dec0). Run folders
`SIMC++_Project/simc_runs/20261007_*_P9.2_L2_{compat,corrected}_<tag>/` (`compare_vs_A/`, `compare_vs_B/`: summary.txt,
summary.json, compare.pdf; `mode_diff_vs_compat.txt` in the corrected folders: event-by-event comparison with the compat
run of the same seed). Binary: commit 3b27d5c.

| run | compat vs A: verdict (F/N), worst p | corrected vs A: verdict (F/N), failing tests | corrected vs compat (same seed) |
|---|---|---|---|
| h_k1_rad_dec0 | PASS 1/257, 9.4e-3 (vs B: PASS 1/257) | FAIL 2/257: CENTRAL.sigcc (F5) + 1 statistical | trees identical; wtcontr identical |
| h_k1_rad_dec1 (π decay) | PASS 2/273, 8.7e-3 | FAIL 3/273: CENTRAL.sigcc (F5) + 2 statistical | 2915 of 200000 events differ in the SHMS recon only (D-037 F7 decay drift); same events, Weight identical |
| h_k1_norad_dec0 | PASS 0/246, 2.2e-2 | FAIL 1/246: CENTRAL.sigcc (F5) | trees identical |
| d_k1_pim_rad_dec1 | PASS 6/279, 7.3e-4 (Q2 ks_w) | FAIL 7/279: CENTRAL.sigcc (F5) + 6 statistical | 2778 events differ in the SHMS recon (F7); Weight identical |
| h_k2_rad_dec1 | PASS 1/253, 6.4e-3 | FAIL 3/253: CENTRAL.sigcc (F5), wtcontr (D-006, pull 43.8) + 1 | wtcontr −13.12 %; sequences diverge after an F7 acceptance change (ΣWeight +0.0009 %) |
| d_k2_pim_rad_dec1 | PASS 5/259, 2.9e-3 | FAIL 7/259: CENTRAL.sigcc, wtcontr (pull 42.6) + 5 | wtcontr −13.59 %; ΣWeight −0.0016 % |
| he3coh_k2_rad_dec1 | PASS 4/253, 6.2e-3 | FAIL 6/253: CENTRAL.sigcc, wtcontr (pull 40.8) + 4 | wtcontr −12.53 %; ΣWeight +0.0013 % |

* **compat**: every run PASSes. The D(e,e'π−) K1 run has the largest count (6 of 279, Q2/q/Q2i/hsyptari, worst p 7.3e-4);
  two more C++ seeds (4243: 3/279, 4244: 0/279, all four variables p > 0.07) show it is a fluctuation of seed 4242.
* **corrected**: the only failures beyond the compat statistics are the two expected deterministic effects:
  * CENTRAL.sigcc (D-038 F5): the central cross section at the real central ε, t, φ: K1 H 0.839106E-06 → 0.925379E-05,
    K1 H norad 0.838925E-06 → 0.922498E-05, K1 D 0.385660E-06 → 0.654696E-05, K2 H 0.449334E-05 → 0.605319E-04, K2 D
    0.822758E-06 → 0.475557E-04, He-3 coh 0.800815E-07 → 0.146748E-05 (a `.hist` number only).
  * K2 wtcontr (D-006, simc.f:237): the hadron (HMS) δ upper cut becomes 10 % instead of the SHMS 22 %. Fortran/corrected
    = +15.42 % (H), +15.69 % (D), +14.30 % (He-3); compat/corrected at the same seed +15.10 %, +15.72 %, +14.33 %; P9.F1
    measured +14.8, +15.6, +14.1 % on the Fortran trees (passcuts_l2.py). Corrected × that factor reproduces the Fortran
    wtcontr within 1.7 σ (H), 0.2 σ (D), 0.6 σ (He-3).
  * Without these two tests the corrected runs count 1/256, 2/272, 0/245, 6/278, 1/251, 5/257, 4/251 failures.
* D-037 F7 (decay drift) moves the recon of ~1.4 % of the K1 events (π decays in the SHMS) by ≤ 1e-2 % δ; in K2 (HMS) it
  changes an acceptance decision early in the run, after which the random sequences diverge: the distributions stay equal
  (ΣWeight within 2e-5).

## 10. Kaon production: L2 in both modes and the effect of each correction (task P10.2)

`simc_compare` (default config) of simc++ seed 4242, 200000 successes, against the P10.F1 Fortran runs
`reference/runs/kaon_<tag>_A` (and `_B` for h_k1_rad_dec0). Run folders
`SIMC++_Project/simc_runs/20261011_*_P10.2_L2_{compat,corrected}_<tag>/` (`compare_vs_A/`, `compare_vs_B/`, `mode_diff_vs_compat.txt`
in the corrected folders). Binary: commit 1275770 (numerically identical to the final commit). Batch:
`validation/scripts/p10_kaon_l2.sh`.

| run | compat vs A: verdict (F/N), worst p | corrected vs A | wtcontr Fortran / compat / corrected | CENTRAL.sigcc corrected/Fortran |
|---|---|---|---|---|
| h_k1_rad_dec0 | PASS 2/268, 3.0e-3 (vs B: PASS 6/268) | FAIL 136/271 | 2.1195 / 2.1162 / 2.8930 | ×4.905 |
| h_k1_rad_dec1 | PASS 3/279, 1.8e-3 | FAIL 106/282 | 2.7957 / 2.7978 / 3.7731 | ×4.905 |
| h_k2_rad_dec0 | PASS 5/253, 2.5e-3 | FAIL 134/256 | 1.2661 / 1.2710 / 1.0985 | ×5.726 |
| h_k2_rad_dec1 | PASS 0/259, 1.8e-2 | FAIL 54/262 | 1.7592 / 1.7568 / 1.4798 | ×5.726 |
| h_k1_sig0_rad_dec0 | PASS 0/268, 1.0e-2 | FAIL 139/271 | 6.8830 / 6.9053 / 8.1603 | ×6.780 |
| d_k1_sigm_rad_dec1 | PASS 0/285, 1.2e-2 | FAIL 166/288 | 9.8722 / 9.8255 / 11.2915 | ×6.781 |
| he3_k1_rad_dec1 | PASS 1/285, 1.5e-3 | FAIL 120/288 | 1.9239 / 1.9094 / 2.4619 | ×4.903 |
| he3coh_k2c_rad_dec1 | PASS 8/258, 4.4e-3 | FAIL 44/261 | 0.65205 / 0.65695 / 0.55462 | ×8.319 |
| h_k1_norad_dec0 | FAIL 10/266, 1.9e-3 (seed 4242; seeds 4243, 4244: PASS 0/266) | FAIL 134/269 | 2.4935 / 2.4981 / 3.5807 | ×4.900 |

* **compat**: every configuration agrees with the Fortran. h_k1_norad_dec0 at seed 4242 has 10 failures, all in correlated
  target-z / SHMS y-slope variables (hsytari, ssytari, ssyptari, ssypfp; worst p 1.9e-3, no gross failure); seeds 4243 and
  4244 give PASS 0/266 with those variables at p > 0.05: a fluctuation of seed 4242 (as P9.2 d_k1_pim).
* **corrected**: in the decay-OFF runs the corrected and compat runs contain the SAME events (same Ntried; only Weight,
  decdist and saghai differ: `mode_diff_vs_compat.txt`); in the decay-ON runs the decay corrections change acceptance decisions
  and the random sequences diverge. The failures are the deterministic effects of the corrections, measured one by one with a
  scratch build in which each correction can be switched off (same seed, 200k; corrected / corrected-without-X − 1 on wtcontr):

| run | D-006 pass_cuts (simc.f:237) | D-043 F5 survival p_kaon | D-006 K→ππ mass (transp.f:158) | D-037 F7 drift |
|---|---|---|---|---|
| h_k1_rad_dec0 | +37.5 % | +3.35 % (ΣWeight all events −0.72 %) | – | – |
| h_k1_rad_dec1 | +35.9 % | – | −1.13 % | +0.01 % |
| h_k2_rad_dec0 | −16.5 % | +0.23 % (ΣWeight +2.26 %) | – | – |
| h_k2_rad_dec1 | −15.9 % | – | +0.16 % | +0.09 % |
| h_k1_sig0_rad_dec0 | +18.9 % | +1.41 % (ΣWeight −2.82 %) | – | – |
| d_k1_sigm_rad_dec1 | +15.9 % | – | −1.37 % | +0.02 % |
| he3_k1_rad_dec1 | +29.7 % | – | −1.16 % | +0.11 % |
| he3coh_k2c_rad_dec1 | −15.8 % | – | +0.13 % | +0.23 % |
| h_k1_norad_dec0 | +43.8 % | +4.06 % (ΣWeight −0.03 %) | – | – |

  * D-006 pass_cuts alone (corrected-without-F5 / compat in the decay-OFF runs): K1 H 1.3227, K2 H 0.8623, Σ0 1.1654, norad
    1.3775; P10.F1 measured on the Fortran trees (passcuts_l2.py) 1/(1−0.245) = 1.3245, 1/1.158 = 0.8636, 1/(1−0.142) = 1.1655,
    1/(1−0.274) = 1.377: the correction is the only difference. (The values in the table, with F5 on, are larger for K1 because
    F5 raises the weight of the high-δ kaons that the corrected cut keeps.)
  * F5 reweights each decay-OFF event by exp(−L·Mk/cτ·(1/p_kaon − 1/p_central)): per recon-Pm bin the corrected/compat ratio runs
    from 1.10 to 0.75 in K1 H, so every weighted distribution changes shape (most of the 134-139 failing tests). The `decdist`
    column (survival probability) follows.
  * D-006 kaon mass: the K→ππ product of a first-half decay in transp gets Mπ instead of Mk: −1.1…−1.4 % of the K1 (SHMS)
    acceptance (Ncontribute/Ntried −1.1 %, 3.6σ, and the SHMS aperture counters), +0.1…0.2 % in the K2 (HMS) runs.
  * D-037 F7: ≤ 0.23 %, within the statistics of the diverged sequences.
  * F1/F2/F3: only the `saghai` column (0 for every Fortran event, the Saghai model value in corrected mode). F6: no effect with
    HMS/SHMS (zaero = 0). F13: CENTRAL.sigcc ×4.9–8.3 (a `.hist` number only).

## 11. Delta region (e,e'p)π: L2 in both modes and the effect of each correction (task P11.1)

`simc_compare` (default config) of simc++ seed 4242, 200000 successes, against the P11.F1 Fortran runs
`reference/runs/delta_h_<tag>_A` (and `_B` for d1_rad). Run folders `SIMC++_Project/simc_runs/20261016_*_P11.1_L2_{compat,corrected}_h_<tag>/`
(`compare_vs_A/`, `compare_vs_B/`, `mode_diff_vs_compat.txt`; the corrected rad folders also hold `scratch_noF2/`: a scratch build with
only D-048 F2 switched off). Batch: `validation/scripts/p11_delta_l2.sh`. Binary: commit ca2d338 (numerically identical to the final commit).

| run | compat vs A: verdict (F/N), worst p | corrected vs A | Ntried Fortran / compat / corrected | wtcontr ×1e7 Fortran / compat / corrected |
|---|---|---|---|---|
| h_d1_rad | PASS 0/245, 1.4e-2 (vs B: PASS 2/245) | FAIL 167/245 (vs B 168/245) | 17612941 / 17644171 / 10261664 | 4.3116 / 4.3059 / 7.5453 |
| h_d1_norad | PASS 0/243, 1.7e-2 | FAIL 5/243 (deterministic only) | 10143696 / 10134383 / 10134383 | 9.5367 / 9.5306 / 9.4863 |
| h_d2_rad | PASS 7/260, 1.0e-3 (seeds 4243, 4244: PASS 0/260) | FAIL 177/260 | 19433914 / 19397195 / 10681851 | 3.8394 / 3.8515 / 7.0874 |
| h_d2_norad | PASS 0/260, 1.0e-2 | FAIL 5/260 (deterministic only) | 10598717 / 10582401 / 10582401 | 8.7355 / 8.7347 / 8.7672 |

* **compat**: every configuration agrees with the Fortran. d2_rad at seed 4242 has 7 failures in correlated HMS (electron) y-slope / Q²
  variables (worst p 1.0e-3); seeds 4243 and 4244 give PASS 0/260: a fluctuation of seed 4242. Open point: the SHMS `hb_out` stop
  fraction of d1_rad (SHMS = electron arm) is lower in all three C++ seeds (0.002601, 0.002641, 0.002603 per try) than in the two Fortran
  samples (0.002643, 0.002660): about −1.4 %, ~3σ on the pooled samples (seed 4244 vs B: FAIL 10/245, hb_out pull 3.3 plus hsdelta/phad
  shape). The montecarlo replay of both captures is bitwise (`l1_montecarlo_replay_delta_*`), so the algorithm is the Fortran one;
  more Fortran seeds would settle it (as D-040/D-041 for the SHMS dc2 counter).
* **corrected, norad**: the corrected and compat runs contain the SAME events (same Ntried, `mode_diff`: only `sigcm` differs, Weight
  identical). The 5 failures are deterministic: CENTRAL.sigcc (D-048 F9: 0.227270E+01 → 0.167664E+02, ×7.38) and the 4 tests of `sigcm`
  (D-048 F3: the Blok fit with the intended arguments instead of the binary128 debris). wtcontr changes by −0.47 % (D1) / +0.37 % (D2):
  D-006 pass_cuts (P11.F1 measured +0.49 % / −0.31 % on the Fortran trees), within 1.7σ of the Fortran.
* **corrected, rad = D-048 F2 (tail-1 Egamma limits)**: the scratch build without F2 (all other corrections on) gives the SAME Ntried as
  compat (17644171, 19397195) and agrees with the Fortran apart from the same deterministic F3/F9 tests (D1 5/245, D2 12/260 = the 5 +
  the 7 seed-4242 fluctuations of compat). F2 alone therefore causes the whole rad difference:
  * Ntried for 200000 successes: ×0.5816 (D1), ×0.5507 (D2), i.e. 41.8 % / 44.9 % of the Fortran tries were lost to the never-radiating
    tail 1 (P11.F1: every ntail = 1 try fails, frac(1) = 0.43-0.46);
  * acceptance per try 0.01134 → 0.01949 (D1, ×1.719), 0.01031 → 0.01872 (D2, ×1.816); wtcontr ×1.758 (D1), ×1.836 (D2);
  * rad/norad wtcontr: 0.452 → 0.795 (D1), 0.440 → 0.808 (D2) (pion K1: 0.878);
  * the shapes change too (radphot, missmass, phad, Weight ...: tail-1 events now exist, 43 % of the radiated events in the
    fortran_diff run), hence the ~170 failing tests.
* D/A Delta (F1): no Fortran reference exists; validated by self-consistency only (`fortran_diff_D048_F1_delta_DA_stop`: D, He-3, C-12,
  300 events each: the vertex energy-momentum balance holds for every event to 1e-13, weights finite and positive, Fermi momentum drawn).
* **Reminder (D-048 F4): none of these yields contains a Delta cross section.**

## 12. SIDIS (e,e'h)X: L2 in both modes and the effect of each correction (task P12.2)

`validation/scripts/p12_sidis_l2.sh` (commit 1d31aa7, simc_compare/simc_mode_diff of the same commit; seed 4242, 100000
successes = the reference ngen) against the P12.F1 references `reference/runs/sidis_*_A` (and `_B` for H π+ K1 rad).
F = tests with p < 0.01 out of N (expected 0.01 N); "ratio" = C++ corrected / C++ compat.

| config | compat verdict (vs A) | compat F/N | compat wtcontr A / C++ (pull) | corrected verdict | corrected wtcontr (ratio) | corrected tree vs compat (simc_mode_diff) |
|---|---|---|---|---|---|---|
| H π+ K1 rad | PASS (vs B: PASS, F = 7) | 8/270 | 95.92 / 97.52 (−2.7) | FAIL (wtcontr only, + decdist/pfermi columns) | 138.26 (1.418) | Weight identical in all 100000 events; decdist and pfermi differ |
| H π+ K1 norad | PASS | 6/253 | 106.69 / 106.48 (+0.4) | FAIL (same) | 151.68 (1.424) | Weight identical; decdist, pfermi |
| H π− K1 rad | PASS | 8/270 | 47.89 / 48.70 (−2.7) | FAIL (same) | 68.01 (1.396) | Weight identical; decdist, pfermi |
| D π+ K1 rad | PASS | 0/270 | 167.53 / 165.23 (+2.3) | FAIL (wtcontr, decdist) | 234.93 (1.422) | Weight identical; decdist in 99999 events |
| H K+ K1 rad | PASS | 5/270 | 7.00 / 7.13 (−3.2) | FAIL (every weighted test) | 10.86 (1.522) | Weight and decdist differ in 95903 events: ΣWeight × 1.068, B/A 1.29 at low Pm → 1.04 at high Pm |
| H π+ K2 rad | PASS | 2/270 | 243.33 / 243.76 (−0.3) | FAIL (wtcontr, decdist, pfermi) | 355.47 (1.458) | Weight identical; decdist, pfermi |
| C-12 π+ K1 rad | PASS | 3/270 | 23.99 / 24.15 (−1.1) | FAIL (wtcontr, decdist, pfermi) | 34.06 (1.410) | Weight identical; decdist, pfermi |
| D π+ K1 rad, do_fermi = 1 | — (the Fortran STOPs at try 27848, D-051 F2) | — | — | runs: 100000 successes in 14783024 tries | — | no reference |

Explanation of the corrected-mode differences (every one is a documented correction; no other difference):
* **wtcontr × 1.40–1.46**: D-006 (simc.f:237 pass_cuts tests the SHMS hadron δ against the electron-arm limit
  SPedge%e%delta%max): the P12.F1 report measured −28.9 … −32.2 % of ΣWeight in the Fortran, i.e. × 1/(1−0.30) ≈ 1.43 in
  the corrected mode. The tree (all successes) is unchanged by D-006.
* **decdist** (every decay-OFF run): D-051 F6 — the survival probability from the hadron momentum and the √ path; for pions it
  is written but not applied to the weight (event.f:1560), so Weight is unchanged.
* **K+ weights**: the same F6 survival enters the kaon weight: ΣWeight × 1.068 (larger for lower kaon momenta, B/A 1.29 at
  recon Pm 1.7 GeV/c), so every weighted K+ distribution moves; the compat K+ run passes. Combined with D-006: × 1.522.
* **pfermi**: D-051 F3 — 0 instead of 0/0 = NaN for H and C-12 (simc_compare skips the all-NaN Fortran column).
* D-051 F1 (NaN → 0 when every PDF is 0), F2, F9, F12, F14 change no event of these runs' trees (F9 is not an ntuple column).
* Compat: the recurring `missmass` / `decdist` / `wtcontr` p ≈ 0.002–0.006 failures of H π± K1 rad come from one C++ sample
  (π+ and π− share the random stream); against seed B the same C++ run gives F = 7 with different variables, and the cap_h /
  cap_d replays are bitwise, so they are fluctuations of the C++ sample, not a model difference.
* P12.3 (D-053): the wtcontr sigma of simc_compare was the relative error of the whole tree yield; it is now that of the
  pass_cuts subset (validation/README.md). Re-compared (simc_runs/20261024_012401_P12.3_wtcontr_recompare/table.tsv): H π+
  K1 rad compat wtcontr pull −2.74 → −1.93 (σ 0.415 → 0.586), H π− −2.74 → −1.94, D π+ +2.28 → +1.61, H K+ −3.22 → −2.27;
  the director's seeds 5201-5204 give −0.07, −0.23, −2.01, −0.36.

## 13. Exclusive rho: L2 in both modes and the effect of the corrections (task P13.2)

SIMC++ 75f8dbd, seed 4242, ngen 100000 (as the references), vs the P13.F1 Fortran runs A (simc_compare with --mode-b;
runs in simc_runs/20261029_014541_P13.2_L2_*). wtcontr = the .hist wtcontr (A: Fortran; B: SIMC++).

| config | compat verdict (F/N) | wtcontr A / compat (pull) | corrected verdict (F/N, gross) | wtcontr corrected / compat | wall (s) compat / corr |
|---|---|---|---|---|---|
| h_k1_rad | PASS 3/279 (vs B: wtcontr pull +0.24) | 52.966 / 53.054 (−0.19) | 2/280, 2 gross | 70.001 (×1.319) | 6911 / 6851 |
| h_k1_norad | FAIL 13/275 (see below) | 61.987 / 61.769 (+0.44) | 13/276, 2 gross | 84.109 (×1.362) | 2871 / 2916 |
| h_k1_pim_rad | PASS 3/279 | 52.966 / 53.054 (−0.19) | 2/280, 2 gross | 70.001 (×1.319) | 6378 / 6353 |
| d_k1_rad | PASS 6/279 | 123.388 / 123.558 (−0.16) | 4/279, 2 gross | 163.480 (×1.323) | 6386 / 6491 |
| he3_k1_rad | PASS 2/279 | 57.206 / 57.763 (−1.13) | 2/279, 2 gross | 75.154 (×1.301) | 6388 / 6359 |
| h_k1_dec1_rad | PASS 3/289 | 50.824 / 51.009 (−0.42) | 3/289, 2 gross | 66.954 (×1.313) | 6527 / 6554 |
| h_k2_rad | PASS 2/280 | 21.286 / 21.407 (−0.63) | 4/281, 2 gross | 32.404 (×1.514) | 7443 / 7423 |

* The two gross failures of every corrected run are the intended corrections: wtcontr (D-006, the hadron delta upper cut
  of pass_cuts: ×1.30–1.51, the Fortran's own −24 … −34 % of D-057 F19 inverted) and the constant `pfermi` column
  (D-051 F3: NaN for every H event in the Fortran, 0 corrected). The other corrections (F1, F2, F3, F6, F9, F15, F18, F20,
  F21) change no tree distribution beyond the statistical tests (F = 2–4 corrected vs 2–6 compat); F2 changes the random
  stream (the corrected and compat trees differ event by event after the first failed generate_rad), F6 has no effect at
  pfer = 0 (D/³He: efer = Md − Mn vs Mp; corrected D AVERAGE.sigcc 5.78e-7 vs Fortran 5.75e-7, −1.2σ).
* h_k1_norad compat F = 13: the failing variables are all functions of the electron energy (hsdelta, hsxpfp, hsyfp,
  hsypfp, hsdeltai, nu, W, epsilon, siglab, Weight, sighad; p = 0.002–0.01) plus phipq, i.e. one correlated fluctuation.
  Two fresh seeds pass (1111: F = 3, P = 0.52; 2222: F = 9, P = 2.0e-3), and the compat run with the RANLUX stream and the
  Fortran seed 12345 (IRandom "ranlux") is BITWISE IDENTICAL to rho_h_k1_norad_A (100000 events, 59 columns, .hist,
  .geni, .gen; simc_runs/20261030_191250_P13.2_norad_check_ranlux_fseed): the seed-4242 result is a statistical
  fluctuation (of the C++ sample or of Fortran A), not a code difference.
