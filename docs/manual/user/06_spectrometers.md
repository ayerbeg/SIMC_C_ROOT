# 6. Spectrometers

## 6.1 Spectrometer numbers

`electron_arm` and `hadron_arm` take **one number per spectrometer** (decision D-032, `include/simc/core/SpectrometerId.h`):

| code | spectrometer | side of the beam | plug-in | stop counters (`.geni`, summary) |
|---|---|---|---|---|
| 1 | HMS | right | `hms` | `hms` |
| 2 | SOS | left | `sos` | `sos` |
| 3 | HRS-R (Hall A) | right | `hrsr` | `hrsr` |
| 4 | HRS-L (Hall A) | left | `hrsl` | `hrsl` |
| 5 | SHMS | left | `shms` | `shms` |
| 7 | calorimeter (NPS-like), HMS side | right | `calo` | `calo` |
| 8 | calorimeter, SHMS side | left | `calo` | `calo` |

* **6 is rejected** when the input is read, in both modes: `simc++: input error: ...: electron_arm = 6 was the retired SHMS
  LSA tune; use 5 = SHMS` (exit status 2). The Fortran STOPs too (dbase.f:895-897), but it still treats 6 as SHMS in some
  places.
* Any other number STOPs as in the Fortran: `I dont know what phi should be for the electron arm` (exit status 3).
* The side of the beam fixes the azimuth of the spectrometer (beam right φ = 3π/2, left φ = π/2) and which ntuple prefix
  holds which particle ([7.1.2](07_outputs.md#712-arm-prefixes-hs-and-ss)). The angle `spec%e%theta`/`spec%p%theta` is
  taken as an absolute value.

## 6.2 The magnetic spectrometers

HMS, SHMS, SOS, HRS-R and HRS-L are the ported Fortran single-arm Monte Carlos (`hms/mc_hms.f`, `shms/mc_shms.f`,
`sos/mc_sos.f`, `hrsr/`, `hrsl/`). Each one does:

* COSY transport of the particle through each magnetic element, with the aperture checks at the element entrance, middle
  and exit;
* decay in flight (π, K) when `doing_decay = 1`;
* multiple scattering in the materials and the wire-chamber (drift chamber) resolution when `mc_smear = 1`;
* the detector-hut apertures (chambers, hodoscopes/scintillators, calorimeter);
* the reconstruction of delta, yptar, xptar and ytar from the smeared focal plane with the reconstruction matrix (with
  the raster position if `correct_raster = 1`).

The COSY maps are the Fortran files, copied in the data directory's `matrices/simc_fortran/` (md5-checked). Where a try is
lost is counted in the stop counters (`.geni`, [7.3](07_outputs.md#73-geni-gen-and-genroot)).

Collimators: `using_HMScoll`/`using_SHMScoll` step pions and muons through the HMS/SHMS collimator (absorption, energy
loss, scattering). The default mode prints the collimator stops as a `.geni` line `Collimator stops (pi/mu)` (the Fortran
prints them nowhere) and counts each slit counter at most once per ray (D-021, D-037 F1, F2).

Hall A HRS: in the Fortran the HRS transport lengths and the 45° rotation of the VDC acceptance are lost by a
precision bug (REAL(16) constants passed as REAL(8)). The FP path length becomes −3e232 cm, every ray in the hut has
"decayed", the kaon survival weight is +Inf, and the VDC test is done without its rotation. The default mode uses the
intended values (D-061 F1, F1e). This changes the HRS acceptance and HRS decay results. Use `--fortran-compat` to get the
Fortran HRS.

## 6.3 The calorimeter arms 7 and 8

The calorimeter arm (`calo`, `calo/mc_calo.f`) is a drift to a calorimeter face at `drift_to_cal` cm from the target, with
the calorimeter aperture. It has no magnetic field, no smearing and no reconstruction (RECON = the values at the spectrometer entrance, after the target). It is
used for:

* **π0 → γγ** (`doing_pizero` with hadron arm 7 or 8): both photons are tracked to the face. `pizero_ngamma = 1` accepts
  the event when at least one photon hits, `2` when both hit. The tree gets the photon positions on the face
  (`xcal_gamma1/2`, `ycal_gamma1/2`, −1e10 for a miss) and the photon four-momenta. In the default mode the photons are
  massless and never decay, and each photon's path starts at 0 (D-062 G1, G3; the Fortran STOPs any π0 run with
  `doing_decay = 1`).
* **an electron in the calorimeter** (`electron_arm` 7 or 8). The default mode never "decays" the electron (D-062 G2).

Kept from the Fortran (expert list): `drift_to_cal` is not checked (0 puts the face at the target, D-062 G5), and the
calorimeter arms reuse the target-window values of the previous trip through the target (D-012 F4). The NPS geometry is
not defined in the Fortran.

## 6.4 `spect_mode`: switching the spectrometer Monte Carlo off

This answers the question "do I simulate the acceptance, or do I throw events at the spectrometer?":

| `spect_mode` | electron arm | hadron arm |
|---|---|---|
| 0 (default) | full Monte Carlo | full Monte Carlo |
| −1 | **off** | full Monte Carlo |
| −2 | full Monte Carlo | **off** |
| 1 | **off** | **off** |

* **Arm on**: the event is transported through the spectrometer: apertures, decay, multiple scattering, resolution and
  reconstruction. Only events that the spectrometer accepts survive. This is the normal "simulate the measured yield" mode.
* **Arm off**: no spectrometer is simulated for that arm. The event is "thrown at the spectrometer" and kept whatever
  the spectrometer would have done. The only limits for that arm are the generation limits (the `SPedge` window, without
  the slop of the spectrometer) and the `SPedge` pass_cuts window that decides whether an event counts in wtcontr. In the
  tree, the arm's `*delta`, `*yptar`, `*xptar` columns then hold the values at the spectrometer entrance (after the target
  energy loss and multiple scattering, no resolution). Its focal-plane columns `*xfp`, `*xpfp`, `*yfp`, `*ypfp` and its
  `*ytar` are **0.0 in every event**: they are not measurements, so leave them out of any analysis of such a run. Its
  `.geni` counters are 0.

So `spect_mode = 1` with `doing_phsp = 1` gives the pure generated phase space, `spect_mode = 0` with `doing_phsp = 1` the
acceptance of both spectrometers ([5.8](05_reactions.md#58-phase-space-and-acceptance-runs-doing_phsp)), and −1/−2 the
acceptance of one spectrometer with the other arm open.

Any other value means "both arms on", as in the Fortran. The default mode prints a WARNING for it (D-065 S2). With an arm
off, the D-006 correction of the hadron delta window (pass_cuts) can change wtcontr by up to −24 % with respect to the
Fortran (`docs/FORTRAN_DIFFERENCES.md` §16): the window is then the only limit on the hadron delta. All spect_mode
settings reproduce the Fortran exactly with `--fortran-compat` and RANLUX.

## 6.5 Polarized-target field (`using_tgt_field`)

With `using_tgt_field = 1`, the particles are tracked through the magnetic field of a polarized target (`trg_track.f`):
from the vertex to the field-free region on the way to the spectrometer, and back for the reconstruction. The field map
comes from `tgt_field_file` (`'trg_field_map.dat'` is in the data directory's `field/`; `'0'` = no field; `' '` = a
uniform 5 T test field). The orientation comes from `targ_Bangle` and `targ_Bphi`. The tree gets 8 more columns with the
angles between q, the scattering plane and the polarization (`th_tarq`, `phitarq`, `beta`, `phis`, `phic` and the vertex
`betai`, `phisi`, `phici`; [7.1.7](07_outputs.md#717-polarized-target-columns-using_tgt_field)).

* The calorimeter arms 7/8 cannot be used with the field (Fortran STOP).
* Default-mode corrections (D-067 T1, T2, T14): the "opposite side" test of the field direction, the sign of the hadron-arm
  angle for targ_Bangle < θp, and the tree columns of π0 and rho runs with the field (the Fortran overwrites the angle
  columns with the photon values).
* Kept (expert list): the step limit of the tracking depends on the sign of the arm (2.5 % of the SHMS events in the
  shipped example stop early), the convergence test, the bending constant 90/E instead of 89.88/E (0.14 % too strong),
  a π0 bent as a charged particle, and a polarization vector that ignores `targ_Bphi` and the sign of `targ_Bangle`.

The shipped `poltar_exclusivepion_example.inp` and the converted `inputfile_protoype.inp` reproduce the Fortran exactly
with RANLUX and `--fortran-compat`.

## 6.6 π0 in the HMS or SHMS

With `doing_pizero` and hadron arm 1 or 5, the Fortran transports the π0 itself through the magnetic spectrometer as if it
were charged (apertures, COSY maps, multiple scattering), and RECON = the values at the spectrometer entrance. SIMC++
does the same in both modes (expert list, D-065 P1), so read such runs as "π0 kinematics inside the HMS/SHMS acceptance",
not as a detector simulation. With `doing_decay = 1` the Fortran STOPs at the first decay test; the default mode treats
the π0 as stable.
