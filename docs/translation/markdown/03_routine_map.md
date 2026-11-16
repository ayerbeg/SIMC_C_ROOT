# 3. Fortran → C++ map

## 3.1 The generated table

[fortran_routine_map.md](fortran_routine_map.md) lists **every routine** of the Fortran files that the Makefile links
into `simc` (172 routines) with its line range, the C++ files that cite it, its status, and the L1 tables named after
it. It is generated from the code, not written by hand:

```bash
python3 -I docs/translation/tools/gen_fortran_map.py      # rewrites docs/translation/fortran_routine_map.md
```

The script reads the Fortran Makefile (`my_objs` + `simc.f` + the module files), finds every SUBROUTINE / FUNCTION /
PROGRAM statement, and searches `src/` and `include/` for lines that contain `<file>.f:<routine>` (or the file name and
the routine name on one line). Routines that have no C++ counterpart on purpose carry a curated reason
(`NOT_PORTED` in the script: dead code, CERNLIB error handling, the HBOOK ntuple, F1F2IN21 branches that SIMC never
reaches with A = 1, ...). At P15.D1: 132 routines cited by C++ code, 40 with a reason, 0 "not cited". **Re-run the
script after a port** and check that no routine is "not cited".

`project_control/TRANSLATION_MAP.md` is the original P0.3 audit table (2026-09-29, H(e,e'p) scope only); it is kept as a
historical record and is superseded by the generated table.

## 3.2 Curated map by subsystem

The main entry points (file → C++ function or class). Line numbers are those of the `// Fortran:` comment at the
C++ definition.

| Fortran | C++ |
|---|---|
| simc.f main l.53-160 / 169-352 / 354-620 | `driver::SimcRun::initialise` / `step` / `finish` (src/driver/SimcRun.cpp) |
| simc.f:montecarlo | `driver::MonteCarlo::run` (MonteCarlo.cpp); π0 photons l.1489-1556 `driver::PizeroCalorimeter` |
| simc.f:calculate_central | `driver::calculateCentral` (Central.cpp) |
| simc.f:report, .geni l.462-532, .gen l.534-612, `inc` | `driver::FormatHist` (RunReportHist.cpp), `FormatGeni` (RunReportGeni.cpp), `FortranHistSink::formatGen`, `FortranHistSink::Inc` (ReportSink.cpp) |
| simc.f l.360-432 normalisation | `driver::Normalise` (RunReport.cpp) |
| simc.f:226-238 pass_cuts | `driver::PassCuts` (SimcRun.cpp) |
| simc.f:120-155 field angles | `driver::TargetFieldArmAngles` (TargetFieldAngles.cpp) |
| dbase.f:dbase_read (CTP read) | `InpReader::loadRun` (src/input/InpReader.cpp), `SimcParameters` (variable table in the same order as regallvars) |
| dbase.f post-read logic l.124-547 | `BuildReactionAndKinematics` (src/setup/SetupReaction.cpp) |
| init.f:limits_init | `BuildLimits` (SetupLimits.cpp) |
| init.f:theory_init, momentum distributions | SetupTheory.cpp |
| init.f:min_max_init | `RunSetup` default member values (`kMinMaxInit`) |
| init.f:target_init, target.f | `material::HallCCryoTarget::targetInit`, `trip` (trip_thru_target), `extremeTripThruTarget`, `targetMusc`, `extremeTargetMusc`, `generateZ` (event.f:185) |
| enerloss_new.f, shared/musc.f, musc_ext.f | `material::SimcMaterial` |
| init.f:radc_init, radc_init_ev, basicrad_init_ev | `radiation::radcInit`, `radcInitEv`, `basicradInitEv` (src/radiation/RadcInit.cpp) |
| radc.f | `radiation::generateRad`, `basicrad`, `peakedRadWeight`, `extradFriedrich`, `extradPhi`, `schwinger`, `spen`, `lambdaDave` (Radc.cpp) |
| brem.f | src/radiation/Brem.cpp (`brem`, `bremos`, `inter`, `inter_prime`, `srad`, `spence`) |
| event.f:generate (beam part l.150-212) | `reactions::StandardBeam::generate` |
| event.f:generate / complete_ev / complete_recon_ev / complete_main, (e,e'p) | `reactions::EepReactionBase` (`completeEv`, `completeRecon`, `completeMain`, `genVolume`) + `HydrogenElastic::generate`, `NuclearEep::generate` |
| event.f, doing_pion / doing_kaon / doing_delta / doing_semi / doing_rho branches | `PionProduction`, `KaonProduction`, `DeltaProduction`, `SidisProduction`, `RhoProduction` (src/reactions/*/) |
| event.f:complete_main l.1432-1563 | `reactions::assembleWeight` (EventChain.cpp) |
| event.f:limits_update, update_range, update_slop_range | `reactions::limitsUpdate`, `updateRange`, `updateSlopRange` |
| event.f:physics_angles, spectrometer_angles; loren.f | `PhysicsAngles`, `SpectrometerAngles`, `loren` (src/setup/KinematicsUtils.cpp) |
| jacobians.f:transform_to_cm | `reactions::transformToCm` |
| physics_proton.f | `proton::sigMott`, `fofa_best_fit`, `sigep`, `deForest` (ProtonCrossSections.cpp); plug-ins in ProtonModels.cpp |
| sf_lookup.f | `nuclear::SfTable` (`lookup`, `lookupDiff`, `generateEm`), theory-file weight (SpectralFunctionTables.cpp) |
| physics_pion.f | `pion::peepi`, `sigBlok`, `sigParam04`, `sigParam3000`, `sigParam2021`, `exclfit`, `sigmaid` (+ `MaidTables`), `SigMultipole` |
| pizero_decay.f | `pion::pizeroDecay` |
| physics_kaon.f, cern/fint.f | `kaon::peeK`, `eekeek`, `eekeeks`, `phspwght`, `sigFactorized`, `fint`, `SaghaiTables` |
| physics_delta.f | `delta::peedelta` (DeltaModel.cpp) |
| semi_physics.f, semi_dilution.f | `sidis::peepiX`, `slacemcfit`, `rhadGlobal`, `semiDilution` |
| cteq5/Ctq5Pdf.f, fdss/fdss.f, F1F2IN21_v1.0.f | `sidis::Cteq5Pdf`, `FdssGrid`, `f1f2in21` / `sf` / `rescsp` / `rescsn` / `resmodp` / `resmodn` |
| rho_physics.f, generate_rho.f, rho_decay.f | `rho::peerho`, `generateRho`, `rhoDecay` (RhoPhysics.cpp) |
| shared/transp.f, project.f, rotate_[hv]axis.f, locforunt.f, mc_*_recon.f | `optics::CosyEngine` (`transp`, `loadForward`, `loadRecon`, recon), `optics::project`, `rotateHaxis`, `rotateVaxis`, `powi` |
| cern/lfit.f | `optics::lfit` |
| hms/mc_hms.f, mc_hms_hut.f, mc_hms_coll.f, pion_coll_absorb.f | `hms::McHms::mcHms`, `hitDipole`, McHmsHut.cpp, McHmsColl.cpp; plug-in `HmsSpectrometer` |
| shms/*.f | `shms::McShms::mcShms` (+ hut, coll); plug-in `ShmsSpectrometer` |
| sos/*.f | `sos::McSos::mcSos` (+ hut); plug-in `SosSpectrometer` |
| hrsl/*.f, hrsr/*.f | `hrs::McHrs::mcHrs` (one class, side parameter; + hut); plug-in `HrsSpectrometer` |
| calo/mc_calo.f, mc_calo_recon.f | `calo::McCalo::mcCalo`, `mcCaloRecon`; plug-in `CaloSpectrometer` |
| trg_track.f | `optics::TargetField` (`init`, `field`, `deriv`, `rk4`, `trackToPlane`, `trackFromTgt`, `trackToTgt`) |
| call_ranlux.f, cern/ranlux.f | `RanluxRandom` (`seed` = RLUXGO, `uniform` = grnd with the skip of every 1000th number) |
| gauss1.f | `gauss1Polar` / `Gauss1` (core/random/Gauss1.h) |
| results_write.f, NtupleInit.f, util/root_tree | `NtupleLayout` (`FillNtuple`, `FillEepNtuple`, `FillPionNtuple`, `FillSemiNtuple`, `FillRhoNtuple`), `RootTreeSink` |
| histograms.inc + simc.f:250-287 | `FortranHistSink` (text `.gen`), `GenHistSink` (`.gen.root`) |

## 3.3 COMMON blocks and SAVEd locals

| Fortran | C++ owner |
|---|---|
| `/gnrl/` (simulate.inc:91) | `RunSetup` (Mh, Ebeam, ...); per try for doing_rho: `MainEvent::gnrl` with `gnrl_per_event` |
| `/target_info/ targ` | `RunSetup::targ` (`TargetInfo`) |
| `/ntuple_variables/ ntup, decdist` | `MainEvent::ntup` (`NtupVars`) |
| `/pfermi/` | `MainEvent::pfermi` (`FermiMomentum`); the tables `eval, pval, mprob` in `RunSetup` |
| `/theory/` | `RunSetup` (theory tables read by the set-up) |
| `/decd/` | `RunSetup` (ctau, doing_decay) and `SpecInput::decay_flag` |
| `/sigkaon/` | `kaon::SaghaiTables`, the kaon model |
| `/radccom/`, `/radctest/`, brem.inc `/bremcom/` | `radiation::RadcState`, `RadcRun` (kernel/RadcState.h) |
| `/sftable/` + SAVEd sf_lookup locals | `nuclear::SfTable` |
| `/track/`, `/driftdistances/` | the spectrometer kernels (`optics::Track`, `CosyEngine`) |
| `/spectrometer1/`, `/spectrometer_shms/`, SOS/HRS/calo counters | the kernels' counters, exposed by `ISpectrometer::stopCounters()` |
| histograms.inc `H` | `FortranHistSink::HistSet` |
| simc.f main-program locals `main, vertex, orig, recon` | `driver::EventState` (one per run) |
| simc.f:1334 SAVEd `tmpfact` | `MonteCarlo::tmpfact_` |
| other SAVEd locals | members of the kernel object that represents the Fortran routine (`EepReactionBase::krel_`, `HallCCryoTarget::thCan_`, ...) |

Zero-initialised COMMON semantics are reproduced by value-initialisation of the members (`{0}`, `{}`).

## 3.4 Naming and comment conventions

* **Origin comment** on every ported function or block, at the definition:
  `// Fortran: event.f:complete_ev (l.432-1054)`. For a block inside a routine: `// l.226-238 pass_cuts` or
  `// Fortran: simc.f:1334`. File headers list the routines of the file. The generator of §3.1 relies on the
  `file.f:routine` form: keep it.
* **Markers** (chapter 5): `FORTRAN-BUG(D-xxx)` (corrected in the default mode; the original code and the reason),
  `FORTRAN-DIFF-CANDIDATE(D-xxx)` (ambiguous, kept, with the question), `FORTRAN-OPEN(D-031 flag 1)` (acceptance mode),
  `FORTRAN-QUIRK` (faithful odd behaviour, often with the finding id `D-xxx Fn`). At P15.D1 the code has 136
  FORTRAN-BUG, 87 FORTRAN-DIFF-CANDIDATE, 72 FORTRAN-QUIRK and 3 FORTRAN-OPEN markers.
* **Names**: C++ names are camelCase versions of the Fortran ones (`peaked_rad_weight` → `peakedRadWeight`,
  `trip_thru_target` → `tripThruTarget`), in a namespace per module (`simc::radiation`, `simc::optics`, `simc::pion`,
  ...). Kernel argument mirrors keep the Fortran names (`McHmsArgs`). Fortran 1-based arrays are either shifted to
  0-based with a comment (`teff(1:3) → [0..2]`) or kept 1-based with an unused element 0 in kernels where the index
  arithmetic matters (`FortranHistSink::Hist::buf[1..50]`).
* Fortran variables of structure fields keep their names in C++ (`main.SP.p.delta`, `targ.Eloss[i].ave`), so a reader
  can grep both sources for the same word.
