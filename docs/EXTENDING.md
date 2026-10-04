# Extending SIMC++ with plug-ins

Status: P7.4 (2026-10-01). This page shows how to add a spectrometer, reaction, physics model, target or output
sink without touching the driver. The worked example is `examples/ideal_box/`, a complete ISpectrometer plug-in
with its test. Read [ARCHITECTURE.md](ARCHITECTURE.md) first for the layers and the per-event data flow.

## 1. The mechanism in three steps

1. **Implement an interface** from `include/simc/stages/` (one class).
2. **Register it under a name**, in its `.cpp`:
   ```cpp
   SIMC_REGISTER(simc::ISpectrometer, "ideal_box", simc::examples::IdealBoxSpectrometer)
   ```
3. **Select it in the JSON pipeline** of a run (the `pipeline` section of a run JSON, or a file passed with
   `simc++ ... --pipeline FILE`). A stage is a name or `{"type": name, <config keys>...}`:
   ```json
   {"pipeline": {"spectrometers": {"hadron": {"type": "ideal_box", "delta_min": -8, "delta_max": 8}}}}
   ```

The driver (`driver::SimcRun`) creates every stage by name through `Registry<I>::create`, calls
`configure(config, params)` with the stage's JSON object (without `"type"`) and the input parameters, then
`init(StageContext)`. It never casts a stage to a concrete class, so a plug-in has exactly the same possibilities
as the built-in stages.

### Build: make the registration reach the executable

`SIMC_REGISTER` registers during static initialisation. The linker drops an object file that nothing references,
so a plug-in must be compiled into an **OBJECT library** that the executable links directly. This is how
`examples/CMakeLists.txt` does it:

```cmake
add_library(simcfw_example_ideal_box OBJECT ideal_box/IdealBoxSpectrometer.cpp)
target_link_libraries(simcfw_example_ideal_box PUBLIC simcfw)
list(APPEND SIMC_EXAMPLE_PLUGINS simcfw_example_ideal_box)
set(SIMC_EXAMPLE_PLUGINS ${SIMC_EXAMPLE_PLUGINS} PARENT_SCOPE)
```

The top-level `CMakeLists.txt` links `${SIMC_DRIVER_LIBS} ${SIMC_EXAMPLE_PLUGINS}` into `simc++`
(option `SIMC_BUILD_EXAMPLES`, default ON). For a new plug-in, either add it to an existing module library
(for example `src/spectrometers/<name>/` with its own OBJECT target added to `SIMC_DRIVER_LIBS` in
`src/driver/CMakeLists.txt`), or follow the `examples/` pattern. A test or a program of your own links
`${SIMC_DRIVER_LIBS}` plus your OBJECT library.

## 2. What every stage has (`IStage`)

```cpp
class IStage {
    virtual std::string name() const = 0;                                   // the registry name
    virtual void configure(const Json& config, const SimcParameters& params) = 0;   // once
    virtual void init(const StageContext& ctx) = 0;                         // once, after configure
};
```

* `config` is the stage's JSON object from the pipeline, merged over the driver's default config for that stage.
  Unknown keys are ignored. Every built-in stage understands `"verbose": false` (the driver adds it for
  `simc++ --quiet`); a plug-in that prints messages should honour it.
* `params` (`SimcParameters`) are the input variables in **.inp units** (degrees, mrad, mg/cm^2 ...).
* `ctx.setup` (`RunSetup`) holds the derived run quantities in **internal units** (MeV, MeV/c, cm, rad, delta in %;
  see `core/Units.h`). Prefer it to `params`.
* `ctx.target` and `ctx.radiator` are the run's target and radiator, through their interfaces. They are null when
  the target or the radiator itself is initialised. Use `ctx.requireTarget(name())` to get a clear error if missing.
* Keep the state of one run inside the object. There are no global variables. Fortran COMMON blocks and SAVEd
  locals become members (see the kernels).

Errors: throw `simc::FortranStop` (a Fortran STOP reached), `simc::NotPorted` (a supported Fortran feature that is
not there yet), `std::invalid_argument` (bad configuration) or `std::logic_error` (misuse). All are in
`core/Errors.h` or the standard library. `simc++` reports them and exits with status 3.

## 3. Spectrometer (`ISpectrometer`) — the worked example

```cpp
struct SpecInput  { double p_spec, th_spec, delta, x, y, z, dxdz, dydz, m2; bool ms_flag, wcs_flag, decay_flag;
                    double fry; };
struct SpecResult { bool ok; std::string stop_reason; double x_fp, dx_fp, y_fp, dy_fp;
                    double delta_recon, dxdz_recon, dydz_recon, y_recon, pathlen;
                    std::optional<double> resolution_factor; double m2_final, decay_distance, Mh2_final; };
class ISpectrometer : public IStage {
    virtual SpecResult transport(const SpecInput& in, IRandom& rng) = 0;
    virtual const std::map<std::string, long>& stopCounters() const = 0;
};
```

The input is the track in the spectrometer TRANSPORT frame (x down, y left, z along the central ray), already
drifted to the plane z = 0 through the target centre: positions in cm, slopes dx/dz and dy/dz, delta in %, the
central momentum in MeV/c and the mass squared in MeV^2. `driver::MonteCarlo` builds it from the event
(simc.f:montecarlo), so a plug-in does not need to know about the target or the beam. Fill `SpecResult`:

* `ok`. If false, set `stop_reason` to the name of the counter that stopped the track.
* The reconstructed target quantities (`*_recon`) and the focal plane (`*_fp`). The driver reads them only if `ok`.
* `resolution_factor`: set it if your detector model defines one (it goes to the ntuple `resfac`). Leave it empty
  otherwise.
* `m2_final`, `Mh2_final`: the mass squared after a possible decay (= `in.m2` without decay). `decay_distance` is in
  cm (0 without decay).
* `stopCounters()`: name → count. It should contain `trials` and `successes`. The counters appear in the
  `.summary.json` (`stops.<name>`) and, for non-Fortran spectrometers, in a generic block of the `.geni`.

`examples/ideal_box/IdealBoxSpectrometer.cpp` (about 90 lines) has perfect resolution and a rectangular
acceptance in delta, xptar and yptar (JSON `delta_min/max` in %, `xptar_min/max` and `yptar_min/max` in rad,
`path_length` in cm). It draws no random numbers. Run it with:

```bash
simc++ /home/ayerbe/NewDisk/simc_gfortran/infiles/eep_h_elastic_rad_A.inp \
       --pipeline examples/ideal_box/hadron_arm.json --ngen 2000 -o out --print-pipeline
```

`examples/ideal_box/both_arms.json` uses it for both arms, with different boxes. Two arms with the same type
**and** the same config share one object (one Fortran COMMON block of counters, as HMS used for both arms).
Different configs give two objects, reported as `ideal_box_hadron` and `ideal_box_electron`.

Test: `examples/ideal_box/test_ideal_box.cpp` (CTest `example_ideal_box_{hadron_arm,both_arms}`, label
`examples`) runs the full driver and checks the pipeline resolution, the counters (trials = successes + cuts),
perfect resolution (RECON = SP − offset, bitwise) inside the box, and the output files. `example_ideal_box_cli`
does the same through `simc++ --pipeline`.

Notes:

* The arm **codes** of the input (`electron_arm`, `hadron_arm`) still select the default plug-in (1 → `hms`,
  5/6 → `shms`). They also select the Fortran target-window materials in `hallc_cryo` (trip_thru_target) and the
  Fortran `.geni` blocks. A plug-in replaces the single-arm Monte Carlo only.
* To reuse the HMS/SHMS optics in a new spectrometer, use the COSY engine `optics::CosyEngine`
  (`include/simc/optics/`; transp, project, recon with the Fortran map files).

## 4. Reaction (`IReaction`)

```cpp
class IReaction : public IStage {
    virtual bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) = 0;   // event.f:generate
    virtual bool completeEv(MainEvent& main, Event& vertex, IRandom& rng) = 0;             // complete_ev
    virtual bool completeRecon(Event& recon) = 0;                                          // complete_recon_ev
    virtual double genVolume() const = 0;                                                  // simc.f:383-396
};
```

* `generate` runs after `IBeam::generate` has set the vertex position and `vertex.Ein`. It generates the reaction
  variables, completes the vertex (`completeEv`), and then either radiates through
  `ctx.radiator->radiate(main, vertex, orig, rng, completeEvFn)` or copies `orig = vertex`. Keep the radiator and
  target pointers from `init(ctx)`.
* `completeEv` must compute the target traversal of the outgoing particles with `ctx.target->tripThruTarget(...)`
  (`main.target.Eloss/teff[1..2]`) and call `radiator->initEvent(main, vertex)` at the end, like event.f:1045.
* `genVolume()` is the generated phase-space volume used in the normalisation.
* Select it with `"reaction": "my_reaction"`. The default comes from the run flags (`hydrogen_elastic` or
  `phase_space`).

Current limitation (to be lifted in P8): the event weight is assembled by `reactions::completeMain` (event.f:
complete_main), which implements only the H(e,e'p)/phase-space spectral-function branch and throws `NotPorted`
for the other Fortran reaction flags. A new reaction therefore runs on top of the H(e,e'p)/phsp flags of the input.
The weight assembly will move behind the reaction interface when the other Fortran reactions are ported.

Template: `include/simc/reactions/HydrogenElastic.h` / `src/reactions/HydrogenElastic.cpp`.

## 5. Physics model (`IPhysicsModel`)

```cpp
class IPhysicsModel : public IStage {
    virtual double sigcc(const Event& vertex, const MainEvent& main) = 0;       // main%sigcc
    virtual double sigccRecon(const Event& recon, const MainEvent& main) = 0;   // main%sigcc_recon
    virtual void beginEvent(long nevent) {}                                     // optional (messages)
};
```

A cross-section model takes the completed vertex (MeV, rad) and returns the cross section in the units of the
Fortran routine it replaces (for H(e,e'p): those of physics_proton.f:sigep). Select it with `"physics": "name"`.
Minimal example:

```cpp
class FlatModel : public simc::IPhysicsModel {
public:
    std::string name() const override { return "flat"; }
    void configure(const simc::Json& c, const simc::SimcParameters&) override { value_ = c.value("value", 1.0); }
    void init(const simc::StageContext&) override {}
    double sigcc(const simc::Event&, const simc::MainEvent&) override { return value_; }
    double sigccRecon(const simc::Event&, const simc::MainEvent&) override { return value_; }
private:
    double value_{1.0};
};
SIMC_REGISTER(simc::IPhysicsModel, "flat", FlatModel)
```
`{"pipeline": {"physics": {"type": "flat", "value": 2.5}}}`. Templates: `src/physics/proton/ProtonModels.cpp`.

## 6. Target (`ITarget`)

```cpp
class ITarget : public IStage, public ITargetLossProvider {
    virtual double generateZ(IRandom& rng) = 0;                                 // vertex z (cm)
    virtual Traversal tripThruTarget(Arm arm, double zpos, double energy, double theta, double mass,
                                     LossMode mode, IRandom& rng) = 0;          // {eloss MeV, radlen}
    virtual ScatteringAngles multipleScattering(double p, double beta, double teff, IRandom& rng) = 0;
    // ITargetLossProvider: the start-up losses of the set-up
    virtual void targetInit(const TargetInitArgs&, TargetInfo&, std::vector<std::string>& messages) = 0;
    virtual void extremeLosses(const ExtremeTripArgs&, TargetInfo&) = 0;
};
```

A target model provides the per-event traversals and the **start-up** quantities that the set-up uses for the
limits and the `.hist` (init.f:target_init: radiation length, average losses, Coulomb energy;
extreme_trip_thru_target: min/max losses and the maximum multiple-scattering angles). `SimcRun` creates a
dedicated instance of the pipeline's target for the set-up, so `targetInit`/`extremeLosses` get everything through
their arguments and must not rely on `configure()`. Select it with `"target": "name"`. Template:
`material::HallCCryoTarget` (`hallc_cryo`).

## 7. Output sink (`IOutputSink`)

```cpp
class IOutputSink : public IStage {
    virtual void beginRun(const RunSetup& setup) = 0;
    virtual void fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                      bool success) = 0;     // every try; success = the event contributes
    virtual void endRun(const RunSummary& summary) = 0;
};
```

`fill` sees every try, including failed ones (the Fortran fills its `geni` histograms with them). Use `success`
to select the contributing events and `main.weight` as their weight. A sink that needs the complete driver report
(central event, limits, resolutions) can also implement `driver::IReportSink::writeReport(setup, report)`
(`include/simc/driver/ReportSink.h`). The pipeline key `sinks` **replaces** the default list:

```json
{"pipeline": {"sinks": ["root_tree", "run_summary", {"type": "my_sink", "file": "out.csv"}]}}
```

Templates: `src/output/RunSummarySink.cpp` (simplest), `src/output/RootTreeSink.cpp`. A sink created in code
(tests, a user program) can be added with `SimcRun::addSink(std::make_unique<MySink>())`, no registration needed.

## 8. Other stages

`IBeam` (`standard`: beam spread, raster, vertex z; template `StandardBeam`), `IRadiator` (`simc_rad`; it must also
report `currentValues()`, `photonEnergies()` and `hadronRadiates()` if the driver should report them, the
defaults are zero/false), `IMaterialModel` (`simc_material`, used inside `hallc_cryo` via its `"material"` key) and
`IRandom` (`mt64`). They follow the same three steps.

## 9. Checklist for a new plug-in

* [ ] One class, one `SIMC_REGISTER` line, in an OBJECT library linked by `simc++` and your tests.
* [ ] Internal units (`core/Units.h`). Convert input units in `configure`/`init` only.
* [ ] No globals. Run state lives in the object, collaborators come from the `StageContext`.
* [ ] If it ports a Fortran routine: cite it (`// Fortran: file.f:routine (l.a-b)`), keep the operation order, put
      the Fortran-mirror structures in a `kernel/` header of its module, and add an L1 test against a Fortran
      reference table.
* [ ] A CTest that runs the driver with the plug-in selected through the JSON (see `examples/ideal_box`).
