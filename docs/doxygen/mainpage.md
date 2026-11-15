# SIMC++ API reference {#mainpage}

SIMC++ is a C++20 / ROOT 6 Monte Carlo for coincidence (e,e'X) experiments in Hall C (and Hall A or calorimeter
configurations). For one reaction it generates events, transports the particles through the target and the
spectrometers, reconstructs them, weights them with a cross-section model and writes ROOT trees, histograms and text
reports. Its physics and numerical results reproduce the Fortran SIMC from which it was translated.

This reference documents every public class, function and variable of the headers in `include/simc`. It is generated
from the header comments by Doxygen (`cmake --build <build> --target doc`).

## The other manuals

| document | for whom |
|---|---|
| [User manual](../manual/user/index.md) | installing and running SIMC++: inputs, reactions, spectrometers, outputs, reproducibility |
| [Developer manual](../manual/developer/index.md) | architecture, data flow, build and tests, writing plug-ins, parallelism, coding standards |
| [Translation report](../translation/index.md) | how SIMC++ was translated from the Fortran SIMC and validated, routine map, Fortran bugs and quirks |

## The framework in one page

* **A run** is a simc::driver::SimcRun. It reads the input (`.inp` or JSON) into simc::SimcParameters, derives the run
  constants (simc::RunSetup, built by simc::BuildRunSetup), creates the **stages** of the run by name, computes the
  central event, loops over tries and finally normalises and writes the outputs.
* **Stages are plug-ins.** Every replaceable piece of the event chain is an abstract interface (module
  @ref stages): random source simc::IRandom, target simc::ITarget, material simc::IMaterialModel, beam
  simc::IBeam, reaction simc::IReaction, cross-section model simc::IPhysicsModel, spectral function
  simc::ISpectralFunction, radiative corrections simc::IRadiator, spectrometer arm simc::ISpectrometer and output sinks
  simc::IOutputSink. An implementation registers a name with #SIMC_REGISTER (simc::Registry) and is selected by
  that name in the JSON `pipeline` section.
* **Dependencies are injected** at `init()` through a simc::StageContext: the run set-up and the collaborators of the
  stage.
* **Events are value types** (simc::Event, simc::MainEvent) held in one simc::driver::EventState for the whole run;
  fields a try does not rewrite keep their previous values (this is part of the defined behaviour).
* **Kernels and stages.** The numerical algorithms are *kernels* (free functions and small state classes, some in
  private `kernel/` headers) that keep the exact operation order and random-number order of the original program;
  *stages* are the classes around them that implement the interfaces.
* **One unit convention** (simc::units): MeV, MeV/c, MeV², cm, rad, g/cm², delta in %, track slopes. The input
  parameters keep their `.inp` units and are converted once by the set-up; the ntuple has its own output units.
* **Two behaviour modes.** The default mode corrects documented defects of the original program; the run option
  `fortran_compat` (simc::RunSetup::fortran_compat, `--fortran-compat`) reproduces the original exactly.
* **Threads.** With `--threads N` every worker owns its own stages and random stream; the results are merged in a fixed
  order. Unless a class says otherwise, one object must be used by one thread at a time; the registries and the
  read-only data tables are safe to share once built.

## Modules

| module | content |
|---|---|
| @ref core | parameters, run set-up, events, run summary, registry, errors, units, constants, random sources |
| @ref input | `.inp` reader, effective values, the run set-up builder |
| @ref stages | the stage interfaces and simc::StageContext |
| @ref reactions | beam and reaction stages, the common event completion |
| @ref physics | cross-section models, spectral functions, PDFs and fragmentation functions, decays; @ref material |
| @ref radiation | radiative corrections (simc::SimcRadiator) and their kernels |
| @ref spectrometers | HMS, SHMS, SOS, HRS and calorimeter arms and their Monte Carlo kernels |
| @ref optics | COSY transport engine, forward chains, target magnetic field |
| @ref output | ROOT tree, summary JSON, histograms, ntuple layout |
| @ref driver | simc::driver::SimcRun, the Monte Carlo of one try, reports, pipeline, command line |
| @ref validation | scripted random numbers and the kernel entry points used by the routine-level tests |

## Conventions of this reference

* `@note Origin:` lines name the Fortran SIMC routine (`file.f:routine`, line numbers where useful) that a function or
  field reproduces, for the reader who wants to compare; the full mapping is in the translation report.
* Units are given in brackets, e.g. `[MeV]`, `[rad]`, `[cm]`; "slope" means dx/dz or dy/dz.
* `FORTRAN-BUG(D-xxx)` / `FORTRAN-QUIRK(D-xxx)` mark behaviour that differs between the default mode and
  `fortran_compat`, catalogued in `docs/FORTRAN_DIFFERENCES.md`.
