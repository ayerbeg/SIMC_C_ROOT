# 5. Reactions

## 5.1 How the reaction is chosen

As in the Fortran (`dbase.f:122-235`), the reaction follows from the flags and from `targ%A`. SIMC++ then picks the
reaction and physics-model plug-ins (`src/driver/Pipeline.cpp`, `PipelineConfig::Defaults`). The first matching row wins:

| flags | target | reaction plug-in | physics model | ntuple layout ([7.1](07_outputs.md#71-the-event-tree-h10)) |
|---|---|---|---|---|
| `doing_phsp = 1` | any | `phase_space` | none (σ = 1) | as the reaction it replaces |
| none of the flags below | A = 1 | `hydrogen_elastic` | `sigep` | (e,e'p), 46 columns |
| none | A = 2 | `deuterium_eep` | `deforest` | (e,e'p) |
| none | A ≥ 3 | `nucleus_eep` | `deforest` | (e,e'p) |
| `doing_pion` | any | `pion_production` | `pion` | pion, 53 (65 with π0) |
| `doing_kaon` | any | `kaon_production` | `kaon` | kaon, 55 |
| `doing_delta` | any | `delta_production` | `delta` | pion layout, 53 |
| `doing_semi` + `doing_pion` or `doing_kaon` | any | `sidis_production` | `sidis` | SIDIS, 56 (68 with π0) |
| `doing_rho` | A = 1, 2, 3 | `rho_production` | `rho` | rho, 59 |

`doing_semi` takes priority over `doing_pion`/`doing_kaon`: they then select the SIDIS hadron (dbase.f:125-132). With
`using_tgt_field`, 8 polarized-target columns are added ([7.1.7](07_outputs.md#717-polarized-target-columns-using_tgt_field)).
The run prints the reaction title, as the Fortran does (for example ` ****--------  H(e,e'pi+)X  --------****`).

Common to all reactions: the beam (energy spread `dEbeam`, size `gen%xwid/ywid`, raster), the target (energy loss,
multiple scattering, Coulomb correction), the radiative tails (`using_rad`), the two spectrometer arms and the
normalisation. The generated volume `genvol` is the product of the generated ranges: electron solid angle, plus the
hadron solid angle and the electron energy for all reactions except H(e,e'p), plus the hadron energy for A(e,e'p) and
SIDIS (simc.f:374-392).

## 5.2 H(e,e'p) elastic, D(e,e'p) and A(e,e'p)

* **H(e,e'p)** (`hydrogen_elastic`, model `sigep`): elastic scattering. Only the electron angles are generated; the
  proton follows from the kinematics. The cross section is the elastic e-p cross section with the form-factor fit of
  `physics_proton.f`.
* **D(e,e'p) and A(e,e'p)** (`deuterium_eep`, `nucleus_eep`, model `deforest`): quasi-elastic knock-out with the de Forest
  off-shell cross section (`deForest_flag`: σcc1, σcc2 or σcc1 on-shell) multiplied by a spectral function:
  * the **theory files** (default): `h2.theory`, `c12.theory`, `fe56.theory`, `au197.theory` in the data directory's
    `spectral/`, chosen from A (2, 12, 56, 197). Any other A uses the C-12 file. SIMC++ prints a WARNING then, while the
    Fortran does so silently (D-033 F1);
  * the **Benhar spectral functions** (`use_benhar_sf = 1`): `benharsf_*.dat`, weighted by `transparency`.

  In a JSON input, `"pipeline": {"reaction": {"spectral_function": {"type": "benhar", "file": "benharsf_12.dat",
  "transparency": 0.6}}}` chooses the file explicitly (`--data-dir` sets its directory). `cuts%Em` limits the generated
  missing energy.

Known limitations (Fortran behaviour, kept in both modes): the D(e,e'p) weight has no hadron-angle jacobian, while
A(e,e'p) has one (≤ 1.5 % at the corners of the generated solid angle, D-033 F4); He-3/He-4 use the C-12 spectral
function unless you give a file. Corrections in the default mode: the Benhar differential spectral function (D-033 F3;
the Fortran divides by the event's own Pm², which diverges at small Pm), and NaN recoil masses at extreme edges (D-033 F7).

## 5.3 Exclusive pion production (`doing_pion`)

`pion_production` with the `pion` model (`physics_pion.f:peepi`): H(e,e'π+)n, D/A(e,e'π±) quasifree (Fermi motion from
the momentum distributions), Δ final states and coherent production, selected with `which_pion`
([4.6](04_inputs.md#46-parameter-reference)). The cross section is the Fortran parametrisation; below W = 2 GeV the
MAID2007 tables are also evaluated (data directory `pion/`).

* **π0** (`doing_pizero = 1`): the π0 decays to γγ (`pizero_decay`). With the calorimeter arm 7/8 as hadron arm, the two
  photons are tracked to the calorimeter, and `pizero_ngamma` (1 or 2) says how many must hit
  ([6.3](06_spectrometers.md#63-the-calorimeter-arms-7-and-8)). The tree has 12 more photon columns.
* **Decay in flight** (`doing_decay = 1`, `ctau`): π → μν in the spectrometer; the `Mhadron` column is the mass of the
  detected particle.
* Restrictions taken from the Fortran: π− from hydrogen, coherent production on H/D, coherent π− on ³He and coherent
  production on ⁴He STOP. The π0 "from the neutron" on hydrogen is refused (default mode, D-038 F2). H(e,e'π−)Δ++ runs in
  the default mode; the Fortran always STOPs it (D-038 F1).

Model limitations, kept in both modes (expert list): coherent production uses the nucleon fit with the photon-nucleus W;
the Δ final states draw a new Breit-Wigner mass when radiative tail 1 recomputes the event; the MAID amplitudes are
single precision and on a coarse grid. The default mode computes the `.hist` CENTRAL.sigcc at the real central
kinematics (the Fortran value is ×11 off for the K1 π+ run, D-038 F5). This is a printed number only, not used in the weights.

## 5.4 Kaon electroproduction (`doing_kaon`)

`kaon_production` with the `kaon` model (`physics_kaon.f:peeK`): K+Λ, K+Σ0, K+Σ− (on the neutron), quasifree on D and A,
and coherent hypernucleus production (`which_kaon` 10/11/12). The weight uses the factorized model (`sigcm`/`factor`
columns); the Saghai model is evaluated for the `saghai` column (tables in the data directory's `kaon/`).

* Decay: with `doing_decay = 1` the kaon decays in flight (70 % μν, 30 % ππ). With `doing_decay = 0` the weight is
  multiplied by the survival probability, which is written to `decdist`.
* Default-mode corrections that change numbers: the survival probability uses the kaon's own momentum, not the central
  momentum (D-043 F5, changes the weights of decay-OFF runs); the `saghai` column is no longer always 0 (D-043 F1, column
  only); the Saghai tables are read correctly (F2, F3, column only); the CENTRAL.sigcc report value (F13).

Model limitations (expert list): fixed decay branching and p*; coherent production uses the nucleon formulas; the proton
mass is used for the Σ− kinematics.

## 5.5 (e,e'p)π in the Delta region (`doing_delta`)

`delta_production` with the `delta` model: H(e,e'p)π0 with the proton detected and the pion missing (`targ%Mrec_struck`
= the charged-pion mass, as in the Fortran).

**The event weight contains no cross section.** In the Fortran the cross-section line is commented out ("Need new xsec
model.", `physics_delta.f:132`). The weight is phase space × virtual-photon flux × the jacobian of the proton angles. The
yields, the `.hist` wtcontr and the `Weight` column of a Delta run are therefore **not physical rates**. Use them for
acceptance and kinematic studies only. This is the same in both modes and is an open question for the experts (D-048 F4).

* D and A targets run in the default mode only. The Fortran STOPs every D/A Delta run because of a wrong flag test
  (D-048 F1). No Fortran reference exists for them; they were validated by energy-momentum balance only.
* Radiation: the Fortran never radiates the incoming electron for Delta runs (every try with tail 1 fails); the default
  mode gives tail 1 its limits (D-048 F2). This changes Ntried and wtcontr of radiative Delta runs.
* The `sigcm` column (a pion fit, not used in the weight) is meaningless in the Fortran (binary128 argument bug, NaN in
  0.4 % of events) and correct in the default mode (D-048 F3).
* `which_pion` 2/3 also act for Delta runs (a Δ Breit-Wigner replaces the recoil pion), as in the Fortran.

## 5.6 Semi-inclusive production (`doing_semi`)

`sidis_production` with the `sidis` model (`semi_physics.f:peepiX`): (e,e'π±)X with CTEQ5M parton distributions and
P. Bosted's 2021 pion fragmentation fit, (e,e'K±)X with fDSS fragmentation, π0 as the average of π+ and π−. D uses the
isospin-symmetric nucleon sum; A ≥ 3 adds the SLAC EMC fit and the hadron-attenuation fit. The tables are in the data
directory's `sidis/`.

| inputs | hadron |
|---|---|
| `doing_semi=1 doing_pion=1 doing_hplus=1` / `0` | π+ / π− |
| `doing_semi=1 doing_pion=1 doing_pizero=1` | π0 (with the calorimeter arm; `pizero_ngamma` required) |
| `doing_semi=1 doing_kaon=1 doing_hplus=1` / `0` | K+ / K− |
| `doing_semi=1` alone | STOP: `I don't have ANY idea what (e,e'm)X we're doing!!!` |

Other inputs: `do_fermi` (D and A only), `doing_decay`/`ctau` (decay OFF: the survival probability goes to `decdist`,
and multiplies the weight for kaons only). `pt_b_param` and `sigc_*` are read and ignored, as in the Fortran.

Default-mode corrections: events whose PDFs are all 0 get weight 0 instead of NaN (D-051 F1); D runs with `do_fermi = 1`
no longer STOP on a negative x_Bj (the event gets weight 0, D-051 F2); the `pfermi` column is 0 instead of NaN for H and
A ≥ 3 (D-051 F3); the π0 d̄ fragmentation function and mass (F4); the survival probability (F6); A ≥ 3 kaon SIDIS
completes its `.hist` (the Fortran STOPs at the end, F14). `doing_semi` with `doing_phsp` is not supported (refused).

Model questions kept (expert list): the 2π/πK threshold formula, the same threshold for kaons, the generated hadron
energy range with radiation, R_h for deuterium.

## 5.7 Exclusive ρ0 production (`doing_rho`)

`rho_production` with the `rho` model (`rho_physics.f:peerho`): p(e,e'ρ0)p, ρ0 → π+π−, with one pion detected
(`doing_hplus` chooses π+ or π−). The ρ mass is drawn from a Breit-Wigner, the decay angles from the ρ decay
distribution. The `Mrho` and `Thrho` columns record them.

* Targets: H; D and ³He with a nucleon at rest. No Fermi motion and no ³He spectral function are loaded, as in the
  Fortran (a model limitation, D-057 F4). A ≥ 4 is refused at input in both modes (the Fortran STOPs).
* Hadron arm 1, 2, 3, 4 or 5; the calorimeter is refused. `doing_rho` together with another reaction flag or with
  `doing_pizero` is refused.
* The efficiency is very low (about 3.5e-4 per try with 4π generation), so rho runs are long (hours for 1e5 events).
* Default-mode corrections: the energy loss of the ρ uses its own angles (D-057 F1); the decay is not called after a
  failed try (F2: this changes the random stream, so a corrected rho run cannot be compared event by event with a
  Fortran run); the yptar sign for far-side pions (F9); NaN guards (F20, F21); the `.hist` title and resolution lines (F15, F18).

## 5.8 Phase-space and acceptance runs (`doing_phsp`)

`doing_phsp = 1` makes the cross section 1 and normfac = 1, for acceptance studies. The weight is then the generation
weight only. Two Fortran properties are kept in both modes, because the question is open with the user (D-031 item 1, D-009):

* **On hydrogen the kinematics stay elastic.** `doing_phsp` keeps `doing_hyd_elast`, so the proton is still
  computed from the elastic kinematics. The run is an elastic acceptance run, not a flat phase space.
* **Radiation stays on.** Only `rad_flag` is set to 0, while `using_rad` remains, so radiative tails are generated. The
  target mass used for the photon limits is 0 (D-028 F1). At high energy this can give zero accepted events (for
  example `acc_phsp_K2` at 6.4 GeV). Set `using_rad = 0` for a pure acceptance run.

`doing_phsp` switches off `doing_pion/kaon/delta/rho` for the physics but keeps the hadron mass they selected.

To measure the acceptance of one spectrometer rather than throwing events at both, see `spect_mode`
([6.4](06_spectrometers.md#64-spect_mode-switching-the-spectrometer-monte-carlo-off)).
