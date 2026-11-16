# 4. Inputs

`simc++` takes one input file, a `.inp` or a `.json` (chosen by the file extension):

```
simc++ <file.inp | file.json> [options]
```

Both formats set the same 106 input variables, with the same names, units and defaults ([4.6](#46-parameter-reference)).
The JSON format can also choose the plug-ins and the random-number generator, and set the run options.

## 4.1 The `.inp` format

### 4.1.1 Syntax

`.inp` files are made of named blocks of `name = value` lines (the "CTP parm" syntax of SIMC input files):

```
; comment
begin parm experiment
  ngen = 50000			;  POS: # of successes; NEG: # of tries
  EXPER%charge = 1.0		;  total charge (mC)
  doing_pion = 1
  extra_dbase_file = 'extra_nucpi'   ; more variables in extra_nucpi.inp
end parm experiment
begin parm debug
  debug(1) = 0
end parm debug
```

* Only the lines between `begin parm <name>` and `end parm` are read. The block name does not matter.
* `;` starts a comment, except inside quotes. Tabs count as spaces.
* Variable names are **case-insensitive** (`ebeam` = `Ebeam`). Structure members are separated by `%` (`targ%thick`).
* `name(i)` is 1-based, `name[i]` 0-based. `name = a, b, c` fills consecutive array elements.
* Values may be expressions (`+ - * / ( )`). A real value given to an integer variable is rounded.
* Flags such as `doing_pion` are integers in the file; a value > 0 means true.
* **Every real value is stored as a 32-bit float first.** For example `ctau = 780.4` becomes 780.4000244140625. This keeps
  the inputs identical to those of the legacy program. `inp2json` writes the rounded values.

The example inputs of this manual are in [examples/](examples/) (installed in `$SIMCPP_PREFIX/share/simc/docs/manual/markdown/user/examples/`).

### 4.1.2 `extra_dbase_file`

`extra_dbase_file = 'name'` loads a second file after the main one; its values override the main file's. `.inp` is
added when the name has no `.`. The file is looked for in the **directory of the main input file**, or in the directory
given with `--infiles-dir DIR`.

### 4.1.3 Old `.` names

Old input files write structure names with `.` instead of `%` (`targ.A`, `spec.e.P`, `SPedge.e.delta.min`, ...). These
names are **not** input variables any more: each one is reported, its value is ignored, and the run uses the defaults,
which usually makes it stop. SIMC++ lists the ignored names. With the example [examples/legacy_names.inp](examples/legacy_names.inp):

```
$ simc++ $EX/legacy_names.inp -o legacy > legacy.log 2> legacy.err; echo "exit status $?"; cat legacy.err
exit status 3
simc++: WARNING: 46 input variable(s) are not SIMC variables and are IGNORED (legacy SIMC behaviour; their defaults are used):
  EXPER.charge, spec.e.P, spec.e.theta, spec.p.P, spec.p.theta, targ.A
  ...
  Legacy '.' names (targ.A, spec.e.P, ...) must be converted to '%' names with simc_gfortran/infiles/convert_inputfile.sh.
simc++: STOP during setup (legacy SIMC behaviour): STOP bad targ.can value
simc++: note: the input has 46 ignored variable(s) (listed above), so their defaults were used — the likely cause of this STOP. ...
```

Convert such a file by replacing the `.` of these names with `%` (keep a copy of the original). This `sed` command does
it for every structure name of the inputs (`targ`, `spec`, `SPedge`, `gen`, `cuts`, `EXPER`):

```bash
sed -E ':a; s/^(\s*(targ|spec|SPedge|gen|cuts|EXPER)[A-Za-z_%]*)\./\1%/; ta' legacy_names.inp > converted.inp
```

Convert the `extra_dbase_file` too. (The SIMC distribution has an equivalent script, `infiles/convert_inputfile.sh`.)
The converted file runs normally (`simc++ converted.inp`).

### 4.1.4 Unknown names

Any name that is not one of the 106 variables, such as a typo or an old name, is reported once per name
(`<name> not registered, registering as int`) and ignored. `inp2json --strict` turns these names into an error (exit
status 3).

## 4.2 The JSON format

A SIMC++ JSON input has up to three sections:

```json
{
  "parameters": { ... the input variables ... },
  "pipeline":   { ... optional: plug-ins, random-number generator, threads ... },
  "run":        { "name": "myrun", "output_dir": "out", "fortran_compat": false }
}
```

### 4.2.1 `parameters`

* The key path is the variable name split on `%`: `targ%thick` becomes `{"targ": {"thick": 723.0}}`, and
  `SPedge%e%delta%min` becomes `{"SPedge": {"e": {"delta": {"min": -10.0}}}}`. Arrays are JSON arrays
  (`"lambda": [0, 0, 0]`, `"debug": [false, false, false, false, false, false]`).
* Values are in the **`.inp` units** (MeV, MeV/c, degrees, mrad, mg/cm², %). `inp2json --table` prints the values after
  the internal unit conversions.
* Flags accept `true`/`false` or numbers (> 0 is true).
* Keys are matched case-insensitively. Keys that start with `_` (for example `"_comment"`) are ignored. Unknown keys
  produce `simc++: warning: unknown parameter '<key>' ..., ignored`.
* A missing key keeps the default ([4.6](#46-parameter-reference)).

The easiest way to write a JSON input is to convert a `.inp`:

```bash
inp2json -o run.json $EX/h_elastic_rad.inp
simc++ run.json -o out_json
```

The converted file holds all 106 variables (with `extra_dbase_file` already merged) and a `_source` block that records
the converted files and any names that were not recognised. A `.inp` and its conversion give identical outputs
([3.4](03_quick_start.md#34-the-same-run-from-json); CTest `cli_json_vs_inp_*`).

Other `inp2json` options: `--infiles-dir DIR` (where to find `extra_dbase_file`), `--no-extra` (do not load it),
`--table[=FILE]` (also print the converted values), `--mapping` (print the `.inp` → JSON table of all variables), `--strict`,
`-q`.

Example configurations come with SIMC++ in `data/config/` (installed in `$SIMCPP_PREFIX/share/simc/data/config/`):
`h_elastic_norad.json`, `h_elastic_rad.json` and `acceptance_phsp_hms_shms.json`, with the three sections written out.

### 4.2.2 `run`

| key | command-line equivalent | default | meaning |
|---|---|---|---|
| `name` | `--run-name NAME` | the input file name without its extension | stem of the output files `<name>.root`, `<name>.hist`, ... |
| `output_dir` | `-o DIR`, `--output-dir DIR` | `.` | output directory (created if needed) |
| `fortran_compat` | `--fortran-compat` | `false` | `true`: the legacy-compatible mode ([8.3](08_running_and_reproducibility.md#83-default-and-legacy-compatible-mode)) |

A command-line option wins over the JSON value. A `.inp` input runs in the default mode unless you give `--fortran-compat`.

## 4.3 The `pipeline` section: choosing plug-ins

Without a `pipeline` section, SIMC++ chooses the components from the input flags. `--print-pipeline` prints the resolved
choice when the run starts (the run then continues). Each entry is a plug-in name, or an object
`{"type": "<name>", <options>}`:

| key | default (from the flags) | registered plug-ins |
|---|---|---|
| `rng` | `mt64`, seeded with `random_seed` | `mt64` (std::mt19937_64), `ranlux` (the legacy SIMC generator); options `"seed"`, and for `ranlux` `"state"` (25 integers) ([chapter 8](08_running_and_reproducibility.md)) |
| `target` | `hallc_cryo` | `hallc_cryo` (cryotarget and solid targets: windows, energy loss, multiple scattering) |
| `beam` | `standard` | `standard` |
| `reaction` | from the flags ([chapter 5](05_reactions.md)) | `hydrogen_elastic`, `deuterium_eep`, `nucleus_eep`, `pion_production`, `kaon_production`, `delta_production`, `sidis_production`, `rho_production`, `phase_space` |
| `physics` | from the flags | `sigep`, `deforest` (and `deforest_cc1`, `deforest_cc1_onshell`, `deforest_cc2`), `pion`, `kaon`, `delta`, `sidis`, `rho` |
| `radiator` | `simc_rad` | `simc_rad` |
| `spectrometers` | `{"electron": ..., "hadron": ...}` from `electron_arm`/`hadron_arm` | `hms`, `shms`, `sos`, `hrsr`, `hrsl`, `calo`; the example `ideal_box` ([chapter 6](06_spectrometers.md)) |
| `sinks` | `root_tree` (only if `Nntu > 0`), `run_summary`, `gen_hist`, `fortran_hist` | the list **replaces** the default list ([chapter 7](07_outputs.md)) |
| `threads` | 1 | number of event-loop threads; `--threads N` wins ([chapter 8](08_running_and_reproducibility.md)) |

Options of the built-in plug-ins (all optional):

| plug-in | options |
|---|---|
| all stages | `"verbose": false` silences the stage messages (`--quiet` sets it for every stage) |
| `reaction` of D/A(e,e'p) | `"spectral_function": {"type": "theory" \| "benhar", "file": "<file in --data-dir>", "transparency": x}` ([5.2](05_reactions.md#52-heep-elastic-deep-and-aeep)) |
| `pion`, `kaon`, `sidis` models | `"data_dir"`: directory of the model tables (default: the data directory's `pion/`, `kaon/`, `sidis/`) |
| `hms`, `shms`, `sos`, `hrsr`, `hrsl` | `"matrix_dir"` (optics maps), and expert overrides of the per-arm switches: `"ms_flag"`, `"wcs_flag"` (multiple scattering and wire-chamber smearing, default `mc_smear`), `"decay"` (default `doing_decay` for the hadron arm), `"ctau"`, `"use_sieve"`; `hms`: `"using_HMScoll"`; `shms`: `"using_SHMScoll"`, `"use_coll"`, `"skip_hb"`, `"use_2sieve"` |
| `calo` | `"drift_to_cal"`, `"ctau"` |
| sinks | `root_tree`: `"file"`, `"tree"`; `run_summary`, `gen_hist`: `"file"`; `fortran_hist`: `"hist"`, `"geni"`, `"gen"` (file names) |

An example: transport the hadron through the ideal box acceptance of the example plug-in instead of the real
spectrometer. `--pipeline FILE` merges a pipeline file over the input's own pipeline:

```
$ simc++ run.json -o out_ideal --ngen 200 --pipeline $SIMCPP_PREFIX/share/simc/examples/ideal_box/hadron_arm.json --print-pipeline --quiet
  ...
  "spectrometers": {
    "electron": { "type": "" },
    "hadron": { "delta_max": 8.0, "delta_min": -8.0, "type": "ideal_box", "xptar_max": 0.05, ... }
  },
```

(`"type": ""` means "the default for this arm code".) Two arms with the same plug-in and the same options share one
spectrometer object, so they share one set of stop counters.

## 4.4 Run-time messages about the inputs

| message | meaning |
|---|---|
| `<name> not registered, registering as int` | `.inp` name unknown: ignored |
| `simc++: WARNING: N input variable(s) are not SIMC variables and are IGNORED ...` | summary of the above, with the hint to convert old names |
| `simc++: warning: unknown parameter '<key>' ...` | JSON key unknown: ignored |
| `WARNING: forcing target mass to Mp!!!` and similar | the set-up adjusted an input (for example the H/D masses) |
| `NOTE: random_state_file '...' is not used by SIMC++` | see `random_state_file` below |
| a WARNING when `pt_b_param` or `sigc_*` are set | these SIDIS inputs are read and ignored |
| a WARNING for `spect_mode` outside {0, 1, −1, −2} | any other value means "both arms on" ([6.4](06_spectrometers.md#64-spect_mode-switching-the-spectrometer-monte-carlo-off)) |
| `simc++: input error: ...` (exit 2) | a rejected input, for example arm code 6 |
| `simc++: STOP during setup ...` (exit 3) | an input check of the set-up failed (the message comes from the original code; [10](10_faq.md#simc-stop-during-setup-legacy-simc-behaviour-stop--exit-3) lists the common ones) |

## 4.5 Units and conventions

* Energies and momenta in MeV and MeV/c; angles in degrees (spectrometer angles, target angle, field angles);
  spectrometer slopes (`xptar`, `yptar`) in mrad; delta in %; lengths in cm; target thickness in mg/cm²; density in
  g/cm³; charge in mC.
* Spectrometer angles are taken as their absolute value. The side of the beam is fixed by the spectrometer number
  ([6.1](06_spectrometers.md#61-spectrometer-numbers)).
* Target coordinates: +x = beam right, +y = up, +z = downstream.
* Defaults: 0 / false / blank, except the `SPedge%*` limits and `cuts%Em` (±1e10).

## 4.6 Parameter reference

All 106 input variables, grouped by topic. **Name** is the `.inp` name. The JSON path is the same name with `%` replaced
by nested objects ([4.2.1](#421-parameters)). **Default** is the value used when the variable is absent. **Used by**
says which runs read it: *all*, *eep* (H/D/A(e,e'p)), *D/A eep*, *pion*, *kaon*, *delta*, *semi* (SIDIS), *rho*,
*phsp*, *rad* (with `using_rad`), *calo* (arms 7/8), *field* (`using_tgt_field`). The test
`docs_user_manual_parameters` fails if an input variable of SIMC++ is missing here.

### Run control

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `ngen` | events | 0 | > 0: number of **accepted** events to generate; < 0: number of **tries** (\|ngen\|). 0 writes nothing. `--ngen N` overrides it. With ngen < 0 the `.hist` wtcontr and AVERAGE.sigcc are divided by the tries | all |
| `random_seed` | integer | 0 | seed of the random-number generator. −1: the current time. For `mt64` a value ≤ 0 means the fixed seed 5489; for `ranlux` it seeds RANLUX at luxury level 3 (0 = RANLUX's default seed). `--seed N` overrides it | all |
| `random_state_file` | path | blank | **Not used** (a note is printed); to start RANLUX from a saved state use `"rng": {"type": "ranlux", "state": [25 integers]}` | — |
| `extra_dbase_file` | file stem | blank | second input file, loaded after the main one ([4.1.2](#412-extra_dbase_file)) | all |
| `Nntu` | flag | 0 | > 0: write the event tree `<run>.root`. Without it no `.root` file is written | all |
| `EXPER%charge` | mC | 0 | beam charge of the simulated run; the luminosity and every normalised yield are proportional to it | all |
| `debug` | 6 flags | all F | diagnostic print-outs: `debug(1)` prints the bremsstrahlung calculation; `debug(2..6)` are accepted but print nothing | — |
| `hard_cuts` | flag | F | apply the `SPedge` limits (and, for (e,e'p), `cuts%Em%max`) as cuts on the reconstructed events: events outside are not accepted and are not written. Without it the limits only decide which events count in the `.hist` wtcontr (pass_cuts) | all |

### Reaction selection and physics models

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `doing_phsp` | flag | F | phase-space (acceptance) run: the cross section is 1 and normfac = 1. On hydrogen the elastic kinematics are kept ([5.8](05_reactions.md#58-phase-space-and-acceptance-runs-doing_phsp)) | phsp |
| `doing_pion` | flag | F | exclusive pion production; with `doing_semi`: SIDIS pions | pion, semi |
| `which_pion` | code | 0 | 0: π+ (π0 with `doing_pizero`) on a proton, quasifree; 1: π− (π0 from the neutron), on a neutron; 2: π+Δ0 (π0Δ+); 3: π−Δ++ (π0Δ0 from the neutron); 10: coherent π+ A→A'; 11: coherent π− | pion |
| `doing_pizero` | flag | F | the pion is a π0 (decays to γγ, detected in the calorimeter arm 7/8 or transported as a π0 through the HMS/SHMS, [6.6](06_spectrometers.md#66-π0-in-the-hms-or-shms)) | pion, semi |
| `pizero_ngamma` | count | 0 | for π0 runs: 1 = at least one photon, 2 = both photons must hit the calorimeter. Must be 1 or 2 with `doing_pizero` (else the run stops) | pion, semi (π0) |
| `doing_kaon` | flag | F | kaon electroproduction; with `doing_semi`: SIDIS kaons | kaon, semi |
| `which_kaon` | code | 0 | 0: K+Λ; 1: K+Σ0; 2: K+Σ− (on a neutron); 10/11/12: coherent production of the Λ/Σ0/Σ− hypernucleus | kaon |
| `doing_delta` | flag | F | (e,e'p)π in the Delta region; the proton is detected. **The weight has no cross section** ([5.5](05_reactions.md#55-eepπ-in-the-delta-region-doing_delta)) | delta |
| `doing_semi` | flag | F | semi-inclusive (e,e'h)X; needs `doing_pion` or `doing_kaon` (and `doing_hplus`) | semi |
| `doing_hplus` | flag | F | SIDIS and rho: the detected hadron is positive (1) or negative (0) | semi, rho |
| `doing_rho` | flag | F | exclusive ρ0 → π+π−, one pion detected (H, D, ³He) | rho |
| `doing_decay` | flag | F | decay in flight of the hadron (π, K) in the spectrometer. Off: kaon runs multiply the weight by the survival probability | pion, kaon, semi, rho |
| `ctau` | cm | 0 | cτ of the hadron for the decay / survival probability (780.4 for pions, 371.3 for kaons) | pion, kaon, semi, rho |
| `do_fermi` | flag | F | SIDIS Fermi motion (D and A only; switched off for hydrogen with a WARNING) | semi |
| `use_benhar_sf` | flag | F | A(e,e'p): use the Benhar spectral function tables instead of the theory files | D/A eep |
| `transparency` | factor | 0 | A(e,e'p) with `use_benhar_sf`: the nuclear transparency in the weight (SF_weight = Z × transparency × SF). With the Benhar tables it must be set (0 gives weight 0) | D/A eep |
| `deForest_flag` | code | 0 | off-shell e-p cross section of D/A(e,e'p): 0 = σcc1, 1 = σcc2, −1 = σcc1 with on-shell kinematics; \|value\| > 1 STOPs | D/A eep |
| `pt_b_param` | GeV⁻² | 0 | SIDIS pT slope: **read and ignored** (the slope is b = 1/(0.120 z² + 0.200)); a WARNING is printed when it is set | — |
| `sigc_flag` | code | 0 | SIDIS "central cross section" binning (0 = z, 1 = pT²): **unused** (the central cross-section calculation is disabled) | — |
| `sigc_nbin` | count | 0 | unused (as `sigc_flag`) | — |
| `sigc_kin_min` | z or GeV² | 0 | unused (as `sigc_flag`) | — |
| `sigc_kin_max` | z or GeV² | 0 | unused (as `sigc_flag`) | — |
| `sigc_kin_ind` | GeV² or z | 0 | unused (as `sigc_flag`) | — |

### Beam and spectrometer settings

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `Ebeam` | MeV | 0 | beam energy | all |
| `dEbeam` | % of Ebeam | 0 | full width of the flat beam-energy spread | all |
| `electron_arm` | code | 0 | spectrometer of the electron: 1 HMS, 2 SOS, 3 HRS-R, 4 HRS-L, 5 SHMS, 7/8 calorimeter. **6 is rejected** (the retired SHMS LSA tune); other values stop the run (`I dont know what phi should be for the electron arm`) ([6.1](06_spectrometers.md#61-spectrometer-numbers)) | all |
| `hadron_arm` | code | 0 | spectrometer of the hadron, same codes | all |
| `spec%e%P` | MeV/c | 0 | central momentum of the electron arm | all |
| `spec%e%theta` | deg | 0 | central angle of the electron arm (absolute value) | all |
| `spec%p%P` | MeV/c | 0 | central momentum of the hadron arm | all |
| `spec%p%theta` | deg | 0 | central angle of the hadron arm (absolute value) | all |
| `spect_mode` | code | 0 | 0: both spectrometer Monte Carlos on; −1: electron arm off (hadron only); −2: hadron arm off (electron only); 1: both off ([6.4](06_spectrometers.md#64-spect_mode-switching-the-spectrometer-monte-carlo-off)) | all |
| `use_first_cer` | flag | 0 | "use first Cherenkov in SHMS": printed in the `.hist`, not used by the simulation | — |

### Beam spot, raster and target position

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `gen%xwid` | cm | 0 | horizontal beam size, one σ (gaussian) | all |
| `gen%ywid` | cm | 0 | vertical beam size, one σ | all |
| `targ%fr_pattern` | code | 0 | fast raster: 1 = old "bedpost" rectangle, 2 = circle, 3 = flat rectangle; other values: no raster | all |
| `targ%fr1` | cm | 0 | raster x half-width (patterns 1, 3) or inner radius (2) | all |
| `targ%fr2` | cm | 0 | raster y half-width (patterns 1, 3) or outer radius (2) | all |
| `targ%xoffset` | cm | 0 | target x offset (+x = beam right) | all |
| `targ%yoffset` | cm | 0 | target y offset (+y = up) | all |
| `targ%zoffset` | cm | 0 | target z offset (+z = downstream) | all |

### Target

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `targ%A` | nucleons | 0 | target mass number. Chooses H (1), D (2) or a nucleus (≥ 3) for every reaction | all |
| `targ%Z` | charge | 0 | target charge. Z < 2.4 means a liquid (cryogenic) target | all |
| `targ%mass_amu` | amu | 0 | target mass (forced to Mp / Md for H / D with a WARNING) | all |
| `targ%mrec_amu` | amu | 0 | recoil (A−1 system) mass for (e,e'p); final-nucleus mass for coherent production. Replaced by the program's estimate when it is more than 100 MeV off | all |
| `targ%rho` | g/cm³ | 0 | target density (with `targ%thick` it gives the target length) | all |
| `targ%thick` | mg/cm² | 0 | target thickness | all |
| `targ%angle` | deg | 0 | angle of a solid target (0 = perpendicular to the beam); forced to 0 for cryotargets; sin(angle) > 0.85 STOPs | all |
| `targ%abundancy` | % | 0 | isotopic purity of the target. **Must be set** (100 for a pure target): the luminosity is divided by abundancy/100, so 0 gives zero luminosity and zero weights | all |
| `targ%can` | code | 0 | cryotarget cell: 1 = beer can, 2 = pudding can, 3 = cryo17. Required for Z < 2.4 (`bad targ.can value` STOP otherwise) | all |

### Polarized target and target field

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `using_tgt_field` | flag | F | track the particles through the polarized-target magnetic field ([6.5](06_spectrometers.md#65-polarized-target-field-using_tgt_field)) | field |
| `tgt_field_file` | path | blank | field map. `'0'` = no field; `' '` = a uniform 5 T test field; otherwise a map file looked for in `--data-dir`, the current directory, then the data directory's `field/` (`trg_field_map.dat` is included). Unset with `using_tgt_field = 1`: refused at start-up | field |
| `targ_Bangle` | deg | 0 | angle of the target field (+ = beam left); also the polarization direction used for the polarized-target angle columns of the tree ([7.1.7](07_outputs.md#717-polarized-target-columns-using_tgt_field)) | field |
| `targ_Bphi` | deg | 0 | side of the beam the field points to (90 = beam left, 270 = beam right) | field |
| `targ_pol` | factor | 0 | target polarization (+ = along B); 0 makes the `th_tarq` column NaN | field |

### Spectrometer offsets

Offsets applied to the event before the spectrometer Monte Carlo of each arm (misalignments).

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `spec%e%offset%x` | cm | 0 | electron arm x offset | all |
| `spec%e%offset%y` | cm | 0 | electron arm y offset | all |
| `spec%e%offset%z` | cm | 0 | electron arm z offset | all |
| `spec%e%offset%xptar` | mrad | 0 | electron arm xptar (vertical slope) offset | all |
| `spec%e%offset%yptar` | mrad | 0 | electron arm yptar (horizontal slope) offset | all |
| `spec%p%offset%x` | cm | 0 | hadron arm x offset | all |
| `spec%p%offset%y` | cm | 0 | hadron arm y offset | all |
| `spec%p%offset%z` | cm | 0 | hadron arm z offset | all |
| `spec%p%offset%xptar` | mrad | 0 | hadron arm xptar offset | all |
| `spec%p%offset%yptar` | mrad | 0 | hadron arm yptar offset | all |

### Spectrometer limits and missing-energy cuts

The `SPedge` limits describe the spectrometer acceptance that you want to simulate. SIMC++ widens them by the
spectrometer "slop" (resolution, multiple scattering, energy loss) to build the **generation limits**.
Reconstructed events inside the unwidened limits count in the `.hist` wtcontr (pass_cuts); with `hard_cuts` they are
also the only accepted events. The defaults (±1e10) are wide open: **set them** in every input. A delta minimum
≤ −100 % is changed to −99.99 %.

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `SPedge%e%delta%min` | % | −1e10 | electron arm delta lower limit | all |
| `SPedge%e%delta%max` | % | 1e10 | electron arm delta upper limit | all |
| `SPedge%e%yptar%min` | mrad | −1e10 | electron arm yptar lower limit | all |
| `SPedge%e%yptar%max` | mrad | 1e10 | electron arm yptar upper limit | all |
| `SPedge%e%xptar%min` | mrad | −1e10 | electron arm xptar lower limit | all |
| `SPedge%e%xptar%max` | mrad | 1e10 | electron arm xptar upper limit | all |
| `SPedge%p%delta%min` | % | −1e10 | hadron arm delta lower limit | all |
| `SPedge%p%delta%max` | % | 1e10 | hadron arm delta upper limit (the legacy-compatible mode uses `SPedge%e%delta%max` here, see [8.3](08_running_and_reproducibility.md#83-default-and-legacy-compatible-mode)) | all |
| `SPedge%p%yptar%min` | mrad | −1e10 | hadron arm yptar lower limit | all |
| `SPedge%p%yptar%max` | mrad | 1e10 | hadron arm yptar upper limit | all |
| `SPedge%p%xptar%min` | mrad | −1e10 | hadron arm xptar lower limit | all |
| `SPedge%p%xptar%max` | mrad | 1e10 | hadron arm xptar upper limit | all |
| `cuts%Em%min` | MeV | −1e10 | A(e,e'p): lower end of the generated missing-energy range (widened by the slop, at least 0). Ignored (wide open) for pion, kaon and Delta runs and when min = max | eep |
| `cuts%Em%max` | MeV | 1e10 | A(e,e'p): upper end of the generated Em range; with `hard_cuts`, (e,e'p) events with recon Em above it are rejected | eep |

### Simulation switches

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `mc_smear` | flag | F | multiple scattering in the target and spectrometer materials, and wire-chamber smearing (the focal-plane resolution) | all |
| `using_Eloss` | flag | F | energy loss in the target and windows | all |
| `correct_Eloss` | flag | F | correct the reconstructed momenta for the most probable energy loss (forced off when `using_Eloss` is off) | all |
| `correct_raster` | flag | F | reconstruct with the raster position (the vertical beam position enters the reconstruction) | all |
| `using_Coulomb` | flag | F | Coulomb correction of the electron energy (forced off for Z = 1) | all |
| `using_HMScoll` | flag | F | step pions/muons through the HMS collimator (absorption, scattering, the `.geni` collimator line) | HMS arm |
| `using_SHMScoll` | flag | F | the same for the SHMS collimator | SHMS arm |
| `drift_to_cal` | cm | 0 | distance from the target to the calorimeter face (arms 7/8). 0 puts the face at the target (not checked) | calo |
| `using_cit_generation` | flag | F | accepted for compatibility, not used | — |

### Radiative corrections

| Name | Units | Default | Meaning | Used by |
|---|---|---|---|---|
| `using_rad` | flag | F | radiative corrections: generate the radiative tails | rad |
| `rad_flag` | code | 0 | radiative option: 0 = best formulas, (ntail, Eγ) basis; 1 = BASICRAD only, (ntail, Eγ) basis; 2, 3 = BASICRAD in the (Eγ1, Eγ2, Eγ3) basis. Set to 0 by `doing_phsp`; other values STOP | rad |
| `extrad_flag` | code | 0 | external radiation: 1 = basic formulas; 2 = basic × φ; 3 = Friedrich approximation; 0 = default (3 for rad_flag 0, 1 otherwise); < 0 STOPs | rad |
| `one_tail` | code | 0 | which tails radiate: 0 all; 1 incoming e; 2 scattered e'; 3 hadron; −1 all but the incoming e; −2 all but e'; −3 all but the hadron. \|one_tail\| > 3 with `using_rad` STOPs | rad |
| `intcor_mode` | code | 0 | internal-correction formula: 0 = Schwinger; otherwise brem (or bremos with `use_offshell_rad`). Normally 1 | rad |
| `use_expon` | flag | 0 | 1 = exponentiate the corrections (brem); normally 0 | rad |
| `use_offshell_rad` | flag | F | use the off-shell bremos formula for the hard correction | rad |
| `lambda` | 3 numbers | 0, 0, 0 | radiative λ(1:3) for rad_flag 4: **overwritten** for every event by the computed values, so the input has no effect | — |
| `Egamma_gen_max` | MeV | 0 | > 0.01: use this fixed maximum photon energy instead of the computed limits | rad |
| `dE_edge_test` | MeV | 0 | widens the energy edges of the generation and radiation limits by this amount (a test of the edges) | all |

## 4.7 Settings that have no `.inp` variable

These have no `.inp` counterpart:

| setting | JSON | command line | default |
|---|---|---|---|
| legacy-compatible mode | `run.fortran_compat` | `--fortran-compat` | off |
| output name and directory | `run.name`, `run.output_dir` | `--run-name NAME`, `-o DIR` | input stem, `.` |
| plug-ins, random-number generator | `pipeline` | `--pipeline FILE` | from the flags ([4.3](#43-the-pipeline-section-choosing-plug-ins)) |
| event-loop threads | `pipeline.threads` | `--threads N` | 1 |
| data directory | — | `--data-dir DIR`, env `SIMC_DATA_DIR` | [2.5](02_installation.md#25-where-the-data-files-are-found) |
| directory of `extra_dbase_file` | — | `--infiles-dir DIR` | the main file's directory |
| number of events, seed | `parameters.ngen`, `parameters.random_seed` | `--ngen N`, `--seed N` (override the input; the `.hist` prints the overridden values) | from the input |
| quiet stages | `"verbose": false` per stage | `--quiet` | verbose |
| print the resolved pipeline | — | `--print-pipeline` | off |
| print the data directory and exit | — | `--print-data-dir` | — |
| usage | — | `-h`, `--help` | — |
