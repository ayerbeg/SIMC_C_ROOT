# 1. Introduction

## 1.1 What SIMC++ is

SIMC++ is a C++20 / ROOT 6 translation of **SIMC**, the Hall C coincidence Monte Carlo (`simc_gfortran`, written in
Fortran). Like SIMC, it simulates a coincidence reaction A(e,e'X) over a phase space matched to the spectrometer
acceptances. It includes:

* the beam and the target: beam energy spread, beam size and raster, energy loss and multiple scattering in the
  target and its windows;
* the reaction: event generation and the physics model cross section;
* radiative corrections, with the three radiative tails of the incoming electron, the scattered electron and the
  hadron;
* the transport through the spectrometers: COSY optics maps, apertures, decay in flight, multiple scattering, the
  focal-plane resolution and the reconstruction;
* the normalisation: luminosity, generation volume and the event weight.

SIMC++ is not a detector simulation (GEANT4) and not a general background generator. These limits are the same as for
SIMC (see the "SIMC is NOT" list in `simc_gfortran/README.md`).

The program is `simc++`. It reads the Fortran input files (`.inp`) without changes, and it also reads a JSON
configuration that can choose the components ("plug-ins") of the simulation. Its outputs are the Fortran text
reports (`.hist`, `.geni`, `.gen`), the same ROOT tree `h10` that `make_root_tree` produces from the Fortran `.bin`
file, and a run summary in JSON. [Chapter 7](07_outputs.md) describes them.

## 1.2 Relation to the Fortran SIMC

The Fortran code is the physics specification. Each routine of the Fortran (`event.f`, `physics_*.f`, `radc.f`,
`target.f`, the spectrometer Monte Carlos in `hms/`, `shms/`, `sos/`, `hrsl/`, `hrsr/`, `calo/`, `trg_track.f`, ...) has a C++
counterpart that implements the same algorithm and cites the Fortran file and line it comes from. Around these ported
routines, every stage of the event chain (random numbers, target, beam, reaction, physics model, radiation,
spectrometers, outputs) is a plug-in selected by name. A Fortran input selects the same components the Fortran would use,
so you do not need to know about plug-ins to run a Fortran input.

What is ported: everything a shipped `simc_gfortran/infiles` input can request.

| reaction (input flags) | targets | chapter |
|---|---|---|
| elastic and quasi-elastic (e,e'p) (no reaction flag) | H, D, A (theory or Benhar spectral functions) | [5.2](05_reactions.md#52-heep-elastic-deep-and-aeep) |
| exclusive pion production (`doing_pion`), π+, π−, π0, Δ final states, coherent | H, D, A | [5.3](05_reactions.md#53-exclusive-pion-production-doing_pion) |
| kaon electroproduction (`doing_kaon`), K+Λ, K+Σ0, K+Σ−, coherent hypernuclei | H, D, A | [5.4](05_reactions.md#54-kaon-electroproduction-doing_kaon) |
| (e,e'p)π in the Delta region (`doing_delta`) | H (D, A in the default mode) | [5.5](05_reactions.md#55-eepπ-in-the-delta-region-doing_delta) |
| semi-inclusive (e,e'π±/π0)X and (e,e'K±)X (`doing_semi`) | H, D, A | [5.6](05_reactions.md#56-semi-inclusive-production-doing_semi) |
| exclusive ρ0 → π+π− (`doing_rho`) | H, D, ³He | [5.7](05_reactions.md#57-exclusive-ρ0-production-doing_rho) |
| phase-space / acceptance runs (`doing_phsp`) | any | [5.8](05_reactions.md#58-phase-space-and-acceptance-runs-doing_phsp) |

The spectrometers are the HMS, SHMS, SOS, HRS-R, HRS-L and the calorimeter arms 7/8 (NPS-like), with π0 → γγ detection.
The polarized-target field (`using_tgt_field`) is also ported ([chapter 6](06_spectrometers.md)).

## 1.3 Two modes: corrected (default) and `--fortran-compat`

The Fortran contains a number of clear bugs. Examples are a wrong array element in the cut that decides whether an
event counts in the `.hist` yield, constants read with the wrong precision, and NaN values written to the ntuple.
The user decided (decision D-031 of the project) that SIMC++ **corrects clear bugs by default**. The exact Fortran
behaviour is still available:

| mode | how | what you get |
|---|---|---|
| corrected (default) | nothing to do | the Fortran physics with the clear bugs fixed. Each fix is listed in `docs/FORTRAN_DIFFERENCES.md` §2 and summarised in [chapter 10](10_differences_from_fortran.md) |
| Fortran-compatible | `simc++ ... --fortran-compat`, or `"run": {"fortran_compat": true}` in a JSON input | exactly what the Fortran computes, bugs included |

Ambiguous cases, where the Fortran might be intentional, are **not** changed in either mode. They are listed for the
experts in `docs/FORTRAN_DIFFERENCES.md` §4 ([10.3](10_differences_from_fortran.md#103-questions-open-for-experts-the-expert-list)).
The run log always says which mode is active:

```
 SIMC++: fortran_compat = false (clear Fortran bugs corrected, docs/FORTRAN_DIFFERENCES.md)
```

## 1.4 What "validated" means

The project uses four levels of evidence. Each level is backed by tests that run in CTest
([chapter 2](02_installation.md#25-running-the-tests)). The numbers below are recorded in the project validation ledger
and in `docs/FORTRAN_DIFFERENCES.md` §7-§17.

1. **L1, routine level, bit for bit.** Each ported routine (cross sections, form factors, COSY transport of a given
   ray, radiative-correction kernels, energy loss, apertures, ...) is evaluated on input tables that were produced by
   running the Fortran routine itself. The C++ must reproduce every output column **bitwise** in the
   Fortran-compatible mode (CTest labels `l1`, `replay`).
2. **Capture replays.** The complete C++ program is fed the random numbers recorded from a Fortran run. It must
   consume them in the same order and reproduce the per-event values, the tree and the text reports exactly.
3. **Full-run identity.** SIMC++ contains a bit-exact port of the Fortran random-number generator (RANLUX). With
   `--fortran-compat` and `"rng": {"type": "ranlux"}`, a complete SIMC++ run of a Fortran input gives the **same files as
   the Fortran run**: the `.hist` (apart from the two time-stamp lines), the `.geni`, the `.gen` and every entry and
   branch of the `h10` tree. This was shown for **all 39 shipped inputs that the Fortran can run** (some of them after
   renaming the legacy input names, see [4.1.3](04_inputs.md#413-legacy--names-and-convert_inputfilesh)), plus the
   spect_mode, π0 and polarized-target probes (CTest label `regression_shipped`; [section 8.4](08_running_and_reproducibility.md#84-reproducing-a-fortran-run-exactly)).
4. **L2, statistical.** Independent runs of the default (corrected) mode with the SIMC++ default generator (mt64) are
   compared with Fortran reference runs using `simc_compare`. The tool runs pulls on the `.hist` normalisation and
   stop counters and weighted χ²/KS tests on every tree variable ([chapter 9](09_comparing_with_fortran.md)). The
   differences that are found are the documented corrections ([chapter 10](10_differences_from_fortran.md)).

The regression tests keep the outputs of a set of standard runs fixed from one version to the next, in both modes.

## 1.5 Limitations

* **Delta-region runs have no cross section in the weight.** The Fortran line that would apply the model cross section
  is commented out ("Need new xsec model.", `physics_delta.f:132`). The `doing_delta` weight is phase space × virtual
  photon flux × a jacobian, so Delta yields are not physical rates. This is the same in both modes. It is an open
  question for the experts.
* **Acceptance mode on hydrogen.** `doing_phsp` on a hydrogen target still uses elastic kinematics, and radiation stays on.
  This is the Fortran behaviour. It is open with the user (D-031 item 1), so SIMC++ keeps it in both modes
  ([5.8](05_reactions.md#58-phase-space-and-acceptance-runs-doing_phsp)).
* **Model limitations kept from the Fortran.** Examples: the coherent pion and kaon production use the nucleon models;
  the ρ0 on D/³He has no Fermi motion; a π0 sent to the HMS/SHMS is transported like a charged particle. These are
  listed in plain words in [10.3](10_differences_from_fortran.md#103-questions-open-for-experts-the-expert-list).
* **Random numbers.** The default generator is `mt64` (std::mt19937_64), not the Fortran RANLUX. A default run is
  therefore statistically equivalent to the Fortran run, not identical. Use the `ranlux` generator to reproduce the
  Fortran stream ([chapter 8](08_running_and_reproducibility.md)). `random_state_file` is not used.
* **No `.bin` file.** SIMC++ writes the ROOT tree directly. There is no Fortran `.bin` file or PAW ntuple.
* **Run-time diagnostics.** The Fortran `debug(2..6)` print-outs are not reproduced ([chapter 4](04_inputs.md#46-parameter-reference)).
