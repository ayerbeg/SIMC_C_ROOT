# 6. Spectrometers

## 6.1 Spectrometer numbers

`electron_arm` and `hadron_arm` take **one number per spectrometer** (`include/simc/core/SpectrometerId.h`):

| code | spectrometer | side of the beam | plug-in | stop counters (`.geni`, summary) |
|---|---|---|---|---|
| 1 | HMS | right | `hms` | `hms` |
| 2 | SOS | left | `sos` | `sos` |
| 3 | HRS-R (Hall A) | right | `hrsr` | `hrsr` |
| 4 | HRS-L (Hall A) | left | `hrsl` | `hrsl` |
| 5 | SHMS | left | `shms` | `shms` |
| 7 | calorimeter (NPS-like), HMS side | right | `calo` | `calo` |
| 8 | calorimeter, SHMS side | left | `calo` | `calo` |

* **6 is rejected** when the input is read: `simc++: input error: ...: electron_arm = 6 was the retired SHMS LSA tune; use
  5 = SHMS` (exit status 2).
* Any other number stops the run: `I dont know what phi should be for the electron arm` (exit status 3).
* The side of the beam fixes the azimuth of the spectrometer (beam right φ = 3π/2, left φ = π/2) and which tree prefix
  holds which particle ([7.1.2](07_outputs.md#712-arm-prefixes-hs-and-ss)). The angles `spec%e%theta`, `spec%p%theta` are
  taken as absolute values.

## 6.2 The magnetic spectrometers

For the HMS, SHMS, SOS, HRS-R and HRS-L, each event goes through:

* the transport through each magnetic element with COSY optics maps, with aperture checks at the element entrance,
  middle and exit;
* decay in flight (π, K) when `doing_decay = 1`;
* multiple scattering in the materials and the wire-chamber resolution when `mc_smear = 1`;
* the detector-hut apertures (chambers, hodoscopes/scintillators, calorimeter);
* the reconstruction of delta, yptar, xptar and ytar from the smeared focal plane with the reconstruction matrix (with the
  raster position if `correct_raster = 1`).

The optics maps are in the data directory's `matrices/simc_fortran/` (md5-checked). Where a try is lost is counted in the
stop counters ([7.3](07_outputs.md#73-geni-gen-and-genroot)).

Collimators: `using_HMScoll` / `using_SHMScoll` step pions and muons through the HMS / SHMS collimator (absorption, energy
loss, scattering). The `.geni` then has a `Collimator stops (pi/mu)` line, and each slit counter counts a ray at most once.

## 6.3 The calorimeter arms 7 and 8

The calorimeter arm (`calo`) is a straight drift to a calorimeter face at `drift_to_cal` cm from the target, with the
calorimeter aperture. It has no magnetic field, no smearing and no reconstruction (RECON = the values at the spectrometer
entrance, after the target). It is used for:

* **π0 → γγ** (`doing_pizero` with hadron arm 7 or 8): both photons are tracked to the face. `pizero_ngamma = 1` accepts the
  event when at least one photon hits, `2` when both hit. The tree gets the photon positions on the face (`xcal_gamma1/2`,
  `ycal_gamma1/2`, −1e10 for a miss) and the photon four-momenta. Photons are massless and never decay.
* **an electron in the calorimeter** (`electron_arm` 7 or 8).

Set `drift_to_cal`: it is not checked, and 0 puts the face at the target. The target windows seen by particles going to
the calorimeter are not modelled specifically (the values of the previous particle's trip through the target are reused).

## 6.4 `spect_mode`: switching the spectrometer Monte Carlo off

This answers the question "do I simulate the acceptance, or do I throw events at the spectrometer?":

| `spect_mode` | electron arm | hadron arm |
|---|---|---|
| 0 (default) | full Monte Carlo | full Monte Carlo |
| −1 | **off** | full Monte Carlo |
| −2 | full Monte Carlo | **off** |
| 1 | **off** | **off** |

* **Arm on**: the event is transported through the spectrometer (apertures, decay, multiple scattering, resolution,
  reconstruction). Only events that the spectrometer accepts survive. This is the normal "simulate the measured yield" mode.
* **Arm off**: no spectrometer is simulated for that arm. The event is "thrown at the spectrometer" and kept whatever the
  spectrometer would have done. The only limits for that arm are the generation limits (the `SPedge` window, without the
  spectrometer's slop) and the `SPedge` window that decides whether an event counts in wtcontr. In the tree, the arm's
  `*delta`, `*yptar`, `*xptar` hold the values at the spectrometer entrance (after the target energy loss and multiple
  scattering, no resolution). Its focal-plane columns `*xfp`, `*xpfp`, `*yfp`, `*ypfp` and its `*ytar` are **0.0 in every
  event**: they are not measurements, so leave them out of any analysis of such a run. Its `.geni` counters are 0.

So `spect_mode = 1` with `doing_phsp = 1` gives the pure generated phase space, `spect_mode = 0` with `doing_phsp = 1` the
acceptance of both spectrometers ([5.8](05_reactions.md#58-phase-space-and-acceptance-runs-doing_phsp)), and −1/−2 the
acceptance of one spectrometer with the other arm open.

Any other value means "both arms on"; a WARNING is printed.

## 6.5 Polarized-target field (`using_tgt_field`)

With `using_tgt_field = 1`, the particles are tracked through the magnetic field of a polarized target, from the vertex
to the field-free region on the way to the spectrometer, and back for the reconstruction. The field map comes from
`tgt_field_file` (`'trg_field_map.dat'` is in the data directory's `field/`; `'0'` = no field; `' '` = a uniform 5 T test
field). Its orientation comes from `targ_Bangle` and `targ_Bphi`. The tree gets 8 more columns with the angles between q,
the scattering plane and the polarization ([7.1.7](07_outputs.md#717-polarized-target-columns-using_tgt_field)).

* The calorimeter arms 7/8 cannot be used with the field (the run stops).
* Approximations to know: the polarization vector is in-plane towards beam left (`targ_Bphi` and the sign of `targ_Bangle`
  are not used for it); the tracking step limit depends on the side of the arm (about 2.5 % of the SHMS events of the
  polarized-target pion example stop early); the bending constant is 90/E (0.14 % stronger than c = 89.88).

## 6.6 π0 in the HMS or SHMS

With `doing_pizero` and hadron arm 1 or 5, the π0 itself is transported through the magnetic spectrometer as if it were
charged (apertures, optics maps, multiple scattering), and RECON = the values at the spectrometer entrance. Read such runs
as "π0 kinematics inside the HMS/SHMS acceptance", not as a detector simulation. The π0 does not decay in flight.
