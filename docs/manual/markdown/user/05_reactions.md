# 5. Reactions

## 5.1 How the reaction is chosen

The reaction follows from the reaction flags and from `targ%A`. SIMC++ then picks the reaction and physics-model plug-ins
(`PipelineConfig::Defaults`, `src/driver/Pipeline.cpp`). The first matching row wins:

| flags | target | reaction plug-in | physics model | tree layout ([7.1](07_outputs.md#71-the-event-tree-h10)) |
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

`doing_semi` takes priority: `doing_pion`/`doing_kaon` then only select the SIDIS hadron. With `using_tgt_field`, 8
polarized-target columns are added ([7.1.7](07_outputs.md#717-polarized-target-columns-using_tgt_field)). The log shows
the reaction title (for example ` ****--------  H(e,e'pi+)X  --------****`).

Common to all reactions: the beam (energy spread `dEbeam`, size `gen%xwid/ywid`, raster), the target (energy loss,
multiple scattering, Coulomb correction), the radiative tails (`using_rad`), the two spectrometer arms and the
normalisation. The generation volume `genvol` is the product of the generated ranges: the electron solid angle; also the
hadron solid angle and the electron energy for every reaction except H(e,e'p); also the hadron energy for A(e,e'p) and
SIDIS.

## 5.2 H(e,e'p) elastic, D(e,e'p) and A(e,e'p)

* **H(e,e'p)** (`hydrogen_elastic`, model `sigep`): elastic scattering. Only the electron angles are generated; the
  proton follows from the kinematics. Elastic e-p cross section with a proton form-factor fit.
* **D(e,e'p) and A(e,e'p)** (`deuterium_eep`, `nucleus_eep`, model `deforest`): quasi-elastic knock-out with the de Forest
  off-shell cross section (`deForest_flag`: σcc1, σcc2 or σcc1 on-shell) times a spectral function:
  * **theory files** (default): `h2.theory`, `c12.theory`, `fe56.theory`, `au197.theory` in the data directory's
    `spectral/`, chosen from A (2, 12, 56, 197). Any other A uses the C-12 file, with a WARNING;
  * **Benhar spectral functions** (`use_benhar_sf = 1`): `benharsf_*.dat`, weighted by `transparency`.

  In a JSON input, `"pipeline": {"reaction": {"spectral_function": {"type": "benhar", "file": "benharsf_12.dat",
  "transparency": 0.6}}}` chooses the file explicitly (`--data-dir` sets its directory). `cuts%Em` limits the generated
  missing energy.

Limitations: the D(e,e'p) weight has no hadron-angle jacobian while A(e,e'p) has one (a ≤ 1.5 % effect at the corners of
the generated solid angle); He-3 and He-4 use the C-12 spectral function unless you give a file.

## 5.3 Exclusive pion production (`doing_pion`)

`pion_production` with the `pion` model: H(e,e'π+)n, D/A(e,e'π±) quasifree (Fermi motion from momentum distributions), Δ
final states and coherent production, selected with `which_pion` ([4.6](04_inputs.md#46-parameter-reference)). The cross
section is a parametrisation of the pion electroproduction data; below W = 2 GeV the MAID2007 tables (data directory
`pion/`) are also evaluated.

* **π0** (`doing_pizero = 1`): the π0 decays to γγ. With the calorimeter arm 7/8 as hadron arm both photons are tracked
  to the calorimeter, and `pizero_ngamma` (1 or 2) says how many must hit ([6.3](06_spectrometers.md#63-the-calorimeter-arms-7-and-8)).
  The tree gets 12 photon columns.
* **Decay in flight** (`doing_decay = 1`, `ctau`): π → μν in the spectrometer; the `Mhadron` column is the mass of the
  particle that reached the detectors.
* Not allowed (the run stops or is refused): π− from hydrogen, coherent production on H or D, coherent π− on ³He,
  coherent production on ⁴He, π0 "from the neutron" (`which_pion = 1`) on hydrogen.

Limitations: coherent production uses the nucleon fit with the photon-nucleus W; for the Δ final states a new
Breit-Wigner Δ mass is drawn when the incoming-electron radiation recomputes the event; the MAID amplitudes are single
precision and taken from the nearest grid point.

## 5.4 Kaon electroproduction (`doing_kaon`)

`kaon_production` with the `kaon` model: K+Λ, K+Σ0, K+Σ− (on the neutron), quasifree on D and A, and coherent hypernucleus
production (`which_kaon` 10/11/12). The weight uses the factorized model (`sigcm` and `factor` columns); the Saghai model
is evaluated for the `saghai` column. Tables in the data directory's `kaon/`.

* With `doing_decay = 1` the kaon decays in flight (70 % μν, 30 % ππ). With `doing_decay = 0` the weight is multiplied by
  the kaon's survival probability to the detectors, which is also written to `decdist`.

Limitations: a fixed decay model (branching ratios and decay momentum); coherent production uses the nucleon formulas;
the proton mass is used in the Σ− kinematics of the factorized model.

## 5.5 (e,e'p)π in the Delta region (`doing_delta`)

`delta_production` with the `delta` model: p(e,e'p)π with the proton detected and the pion missing (the missing particle
has the charged-pion mass). D and nuclei are supported in the default mode.

**The event weight contains no cross section** (no Delta model is implemented). The weight is phase space × virtual-photon
flux × the jacobian of the proton angles. The yields, the `.hist` wtcontr and the `Weight` column of a Delta run are
therefore **not physical rates**; use Delta runs for acceptance and kinematic studies only. The `sigcm` column holds a
pion-production fit evaluated at the proton's CM angles, for information; it is not in the weight.

`which_pion` 2 or 3 also acts on Delta runs: a Breit-Wigner Δ mass then replaces the recoil pion.

## 5.6 Semi-inclusive production (`doing_semi`)

`sidis_production` with the `sidis` model: (e,e'π±)X with CTEQ5M parton distributions and a 2021 fit of the pion
fragmentation functions (P. Bosted), (e,e'K±)X with fDSS fragmentation, π0 as the average of π+ and π−. Deuterium uses the
isospin-symmetric nucleon sum; A ≥ 3 adds the SLAC EMC fit and a hadron-attenuation fit. Tables in the data directory's
`sidis/`.

| inputs | hadron |
|---|---|
| `doing_semi=1 doing_pion=1 doing_hplus=1` / `0` | π+ / π− |
| `doing_semi=1 doing_pion=1 doing_pizero=1` | π0 (with the calorimeter arm; `pizero_ngamma` required) |
| `doing_semi=1 doing_kaon=1 doing_hplus=1` / `0` | K+ / K− |
| `doing_semi=1` alone | the run stops: `I don't have ANY idea what (e,e'm)X we're doing!!!` |

Other inputs: `do_fermi` (Fermi motion, D and A only), `doing_decay`/`ctau` (with decay OFF the survival probability goes
to `decdist` and multiplies the weight for kaons only). `pt_b_param` and `sigc_*` are read and ignored. Events where the
parton distributions vanish, or where Fermi motion gives x_Bj ≤ 0, get weight 0. `doing_semi` with `doing_phsp` is refused.

Limitations: the 2π / πK threshold uses an approximate missing-mass formula, also for kaons and with Fermi motion; with
radiation the generated hadron-energy range is wide; the hadron-attenuation factor for deuterium is not 1.

## 5.7 Exclusive ρ0 production (`doing_rho`)

`rho_production` with the `rho` model: p(e,e'ρ0)p, ρ0 → π+π−, with one pion detected (`doing_hplus` chooses π+ or π−). The ρ
mass is drawn from a Breit-Wigner and the decay angles from the ρ decay distribution; the `Mrho` and `Thrho` columns record
them.

* Targets: H; D and ³He with a nucleon at rest (no Fermi motion, no ³He spectral function). A ≥ 4 is refused.
* Hadron arm 1, 2, 3, 4 or 5; the calorimeter is refused. `doing_rho` together with another reaction flag or with
  `doing_pizero` is refused.
* The ρ is generated over 4π, so the efficiency is low (about 3e-4 accepted events per try): runs of 1e5 events take hours
  (use `--threads`).

Limitations: s without the Fermi terms, t at the ρ pole mass, decay angles in the lab-axis frame rather than the helicity
frame, radiation of the neutral ρ (most of its tail-3 tries fail), and a new ρ drawn when the incoming-electron radiation
recomputes the event.

## 5.8 Phase-space and acceptance runs (`doing_phsp`)

`doing_phsp = 1` sets the cross section to 1 and normfac to 1, for acceptance studies: the weight is then the generation
weight only. Two properties to know:

* **On hydrogen the kinematics stay elastic.** The proton is still computed from the elastic kinematics, so the run is an
  elastic acceptance run, not a flat phase space.
* **Radiation stays on** if `using_rad = 1`, with the photon limits computed for a target mass of 0. At high energy this
  can give zero accepted events. Set `using_rad = 0` for a pure acceptance run.

`doing_phsp` switches off the pion, kaon, Delta and rho physics but keeps the hadron mass that these flags selected.

To measure the acceptance of one spectrometer, or to throw events at a spectrometer without simulating it, see
`spect_mode` ([6.4](06_spectrometers.md#64-spect_mode-switching-the-spectrometer-monte-carlo-off)).

(How each reaction behaves in the legacy-compatible mode, and what the default mode corrects, is listed in
`docs/translation/user_differences.md`.)
