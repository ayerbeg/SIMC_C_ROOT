# 4. Extending SIMC++

The mechanics of every stage interface (signatures, configuration keys, the in-tree and out-of-tree build) are in
[EXTENDING.md](../../EXTENDING.md). This chapter gives the step-by-step recipes with code that is compiled and tested:
`examples/dev_manual/` (CTest `example_dev_manual_plugins`), `examples/ideal_box/` (`example_ideal_box_*`) and
`examples/out_of_tree_plugin/` (`package_out_of_tree_plugin`).

## 4.1 The three steps (every plug-in)

1. Implement one interface from `include/simc/stages/` (§3.5 of chapter 3 gives the contracts).
2. Register it under a name in its `.cpp`: `SIMC_REGISTER(simc::<Interface>, "<name>", <Class>)`, and compile it into an
   **OBJECT** library that the executable links directly (chapter 2 §2.3).
3. Select it in the JSON `pipeline` of a run (or a file given to `simc++ --pipeline FILE`). A stage is `"name"` or
   `{"type": "name", <config keys>}`; the keys reach `configure()`.

Where to put the code:

* **in tree** (a new built-in stage): `src/<module>/` with its own `CMakeLists.txt` (OBJECT library), added with
  `add_subdirectory` in the top-level CMakeLists.txt and appended to `SIMC_DRIVER_LIBS` in src/driver/CMakeLists.txt; the
  public header in `include/simc/<module>/`, private kernel state in `include/simc/<module>/kernel/`; then add the header
  to the reference table of chapter 3 §3.9 (the `docs_dev_manual` test checks it).
* **examples**: an OBJECT library in examples/CMakeLists.txt; append it to `SIMC_EXAMPLE_PLUGINS` to make it available in
  `simc++` (like `ideal_box`), or link it only into its test (like `dev_manual`).
* **out of tree**: compile the plug-in source into your own executable that links `SIMC::simc` from the installed
  package and calls `simc::driver::CommandLineMain` (§4.6).

## 4.2 A new physics model (`IPhysicsModel`)

The model receives the completed vertex and returns the cross section in the units the reaction expects
(for H(e,e'p): those of the `sigep` model). Complete example:

<!-- from: examples/dev_manual/FlatModel.cpp -->
```cpp
class FlatModel : public simc::IPhysicsModel {
public:
    std::string name() const override { return "example_flat"; }
    void configure(const simc::Json& config, const simc::SimcParameters&) override {
        if (config.is_object()) value_ = config.value("value", value_);
    }
    void init(const simc::StageContext&) override {}
    double sigcc(const simc::Event&, const simc::MainEvent&) override { return value_; }
    double sigccRecon(const simc::Event&, const simc::MainEvent&) override { return value_; }

private:
    double value_{1.0};
};

}  // namespace mylab

SIMC_REGISTER(simc::IPhysicsModel, "example_flat", mylab::FlatModel)
```

Select it with `{"pipeline": {"physics": {"type": "example_flat", "value": 2.5}}}`. The test runs H(e,e'p) with it and
checks that every success has `sigcc = 2.5`. A model that must also write event fields (CM angles, jacobians,
`ntup.sigcm`, as the `pion` model does) overrides `sigccMain(vertex, main&)`; the reactions call that one. Models
that read tables take a `data_dir` key and default to `simc::DataDir("<sub>")`, and read the tables in `init()` (a missing
file is a start-up error). Templates: src/physics/proton/ProtonModels.cpp, src/physics/pion/PionModel.cpp.

Put the numerical part in a kernel function (`double peeX(const Event&, MainEvent&, const XConfig&)`) with a config
struct for the run constants (and `fortran_compat` if the model has a compat behaviour), and test it on fixed inputs. A
model ported from the Fortran also needs a reference table (see the
[translation report](../../translation/02_validation.md#27-adding-a-new-fortran-reference-end-to-end)).

## 4.3 A new output sink (`IOutputSink`, `IParallelSink`)

`fill` sees **every** try (use `success` to select contributing events, `main.weight` as their weight); `endRun` gets
the `RunSummary`. To work with `--threads N`, also implement `IParallelSink`: `makeWorkerSink` returns a new, configured
sink (it gets `beginRun` and `fill`, never `endRun`; return `nullptr` if `fill` does nothing) and `mergeWorkerSink` adds a
worker's content, called in the worker order 0..N−1 before `endRun`:

<!-- from: examples/dev_manual/CountSink.cpp -->
```cpp
    void fill(const simc::MainEvent& main, const simc::Event&, const simc::Event&, const simc::Event&,
              bool success) override {
        ++tries_;
        if (success) {
            ++successes_;
            sumWeight_ += main.weight;
        }
    }
    void endRun(const simc::RunSummary& summary) override {
        std::ofstream f(file_);
        if (!f) throw std::runtime_error("example_count: cannot write " + file_);
        f << "tries " << tries_ << "\nsuccesses " << successes_ << "\nsum_weight " << sumWeight_ << "\nntried "
          << summary.ntried << "\n";
    }

    // --threads N: a fresh sink per worker (fill only), merged in the worker order 0..N-1 before endRun().
    std::unique_ptr<simc::IOutputSink> makeWorkerSink(int, int) const override {
        auto w = std::make_unique<CountSink>();
        w->file_ = file_;   // configured like this one; a worker sink never writes (no endRun)
        return w;
    }
    void mergeWorkerSink(simc::IOutputSink& sink) override {
        const auto& w = dynamic_cast<const CountSink&>(sink);
        tries_ += w.tries_;
        successes_ += w.successes_;
        sumWeight_ += w.sumWeight_;
    }
```

The pipeline key `sinks` **replaces** the default list, so name the defaults you still want:
`{"pipeline": {"sinks": ["run_summary", "example_count"], "threads": 2}}`. The test checks that the one-thread counts equal
the run totals and that the two-thread merged counts equal `ngen` and the merged `ntried`. A sink created in code is added
with `SimcRun::addSink(std::make_unique<MySink>())` (no registration). A sink that needs the full report (central event,
limits) also implements `driver::IReportSink::writeReport`. Templates: src/output/RunSummarySink.cpp (simplest),
RootTreeSink.cpp (worker temporary files), GenHistSink.cpp (histograms summed).

## 4.4 A new spectrometer (`ISpectrometer`)

The input is the track in the TRANSPORT frame at z = 0 (cm, slopes, delta %, central momentum MeV/c, mass² MeV²);
return `SpecResult` with `ok`/`stop_reason`, the reconstructed target quantities, the focal plane, the path length,
optionally `resolution_factor`, and `m2_final`/`Mh2_final` (= `in.m2` without decay). `stopCounters()` must contain
`trials` and `successes`. The minimal out-of-tree example:

<!-- from: examples/out_of_tree_plugin/ToyWindowSpectrometer.cpp -->
```cpp
    simc::SpecResult transport(const simc::SpecInput& in, simc::IRandom&) override {
        simc::SpecResult r;
        r.m2_final = in.m2;
        r.Mh2_final = in.m2;
        ++counts_["trials"];
        if (std::fabs(in.delta) > deltaMax_) {
            ++counts_["delta"];
            r.stop_reason = "delta";
            return r;
        }
        ++counts_["successes"];
        r.ok = true;
        r.delta_recon = in.delta;
        r.dxdz_recon = in.dxdz;
        r.dydz_recon = in.dydz;
        r.y_recon = in.y;
        r.x_fp = in.x;
        r.dx_fp = in.dxdz;
        r.y_fp = in.y;
        r.dy_fp = in.dydz;
        r.resolution_factor = 1.0;
        return r;
    }
```

Select it per arm: `{"pipeline": {"spectrometers": {"hadron": {"type": "toy_window", "delta_max": 8.0}}}}`. Arms with the
same type and the same config share one object (one counter block). The arm *codes* of the input still select the
default plug-in, the target-window materials in `hallc_cryo` and the built-in `.geni` blocks; a plug-in replaces the
single-arm Monte Carlo only, and its counters appear in a generic `.geni` block and in the summary JSON. For real optics
reuse `optics::CosyEngine` (COSY maps, `transp`, `project`, recon); put the kernel state (track, counters,
argument mirrors) in `include/simc/spectrometers/<name>/kernel/` behind a pimpl, as `HmsSpectrometer` does. A spectrometer ported
from the Fortran also needs reference ray tables and a full-run identity test (translation report). For the target-field tracking implement `reconstructLastTrack`.

## 4.5 A new reaction (`IReaction`)

A reaction is the largest plug-in: it owns the generation, the vertex completion and the weight. Steps:

1. **Choose a base.** An (e,e'p)-like reaction derives from `reactions::EepReactionBase` and implements `generate()` and
   `spectralWeight()`; `completeEv`, `completeRecon`, `completeMain` and `genVolume` are shared by the family.
   The structure, from `include/simc/reactions/HydrogenElastic.h`:

   <!-- from: include/simc/reactions/HydrogenElastic.h -->
   ```cpp
   class HydrogenElastic : public EepReactionBase {
   public:
       explicit HydrogenElastic(std::string name = "hydrogen_elastic");

       void init(const StageContext& ctx) override;
       bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) override;
       double spectralWeight(const Event& vertex) override;
   };
   ```

   A reaction with its own kinematics derives from `IReaction` directly or from `PionProduction` (its own `completeEv`
   and `completeMain`), like the kaon, Delta and SIDIS reactions; rho is the example of a reaction that changes the hadron
   mass per try (`MainEvent::gnrl`).
2. **`init(ctx)`**: refuse set-ups it cannot handle (check `ctx.setup.flags`, throw `std::invalid_argument` or
   `NotPorted`), keep `ctx.target`, `ctx.radiator`, `ctx.physics` (`requireTarget(name())` ...).
3. **`generate`** (called after `IBeam::generate`): set `main.gen_weight`, draw the generated variables from the
   generation limits in `RunSetup::gen` (every draw through `rng`, in a fixed order), call `completeEv`, then
   `radiator->radiate(main, vertex, orig, rng, completeEvFn)` when `using_rad`, else `orig = vertex`.
4. **`completeEv`**: the vertex kinematics, `target->tripThruTarget` for arms 2 and 3 (`main.target.Eloss/teff`), and
   finally `radiator->initEvent(main, vertex)`.
5. **`completeRecon`**: the reconstructed kinematics from the spectrometer quantities in `recon`.
6. **`completeMain`**: `main.SF_weight`, then `sigcc`/`sigcc_recon`/`weight` with the model (for the (e,e'p) family
   `reactions::assembleWeight(force_sigcc, main, vertex, vertex0, recon, model, MainContext::FromSetup(setup))`); return
   false when the event cannot contribute unless `force_sigcc` (the central event).
7. **`genVolume`**: the generated phase-space volume.
8. Optional: `prepareCentral`, `finalTargetState`, `beginTry`.
9. **Register and select**: `SIMC_REGISTER(IReaction, "my_reaction", MyReaction)`, `{"pipeline": {"reaction":
   "my_reaction", "physics": "my_model"}}`. `PipelineConfig::Defaults` picks the built-in reactions from the `doing_*`
   flags; a new reaction is selected through the JSON (or add its flag rule to `Defaults` in src/driver/Pipeline.cpp).
   The ntuple layout follows the reaction flags (output/NtupleLayout.h); a reaction with new columns needs a layout there.
10. **Tests**: a run through the driver (like `example_dev_manual_plugins`), tests of the kernels on fixed inputs (with a
    `ScriptedRandom`), a regression baseline in both modes (`tests/regression`) and a statistical comparison with
    `simc_compare`. A reaction ported from the Fortran needs the reference tests of the translation report.

## 4.6 Out of tree

<!-- from: examples/out_of_tree_plugin/CMakeLists.txt -->
```cmake
find_package(SIMC 2.0 REQUIRED)

add_executable(my_simc main.cpp ToyWindowSpectrometer.cpp)
target_link_libraries(my_simc PRIVATE SIMC::simc)
```

<!-- from: examples/out_of_tree_plugin/main.cpp -->
```cpp
#include "simc/driver/CommandLine.h"

int main(int argc, char** argv) { return simc::driver::CommandLineMain(argc, argv); }
```

`cmake -S my_plugin -B build -DCMAKE_PREFIX_PATH="<prefix>;<ROOT install>"`, then `build/my_simc run.inp --pipeline
toy_window.json -o out`. The installed data are found without `--data-dir` (chapter 2 §2.7). Only the installed public
headers are available (no `kernel/` headers). Plug-ins are linked, not loaded at run time (no dlopen).

## 4.7 Rules: state and thread safety

* **No global or static mutable state.** Values kept between calls become members of the object; caches
  are members too. With `--threads N` every worker builds its own instance of every stage, so a `static` cache would be
  shared between threads (a data race) and between runs. The only process-wide mutable state is the Registry (filled at
  static initialisation, read-only afterwards) and the column lists of `NtupleLayout` (under a mutex).
* A stage takes its collaborators from the `StageContext` and never casts them to concrete classes.
* All random numbers come from the `IRandom` passed to the call, never from another generator.
* A worker sink sees only the `fill` arguments; it must not write the run's files or use the master's collaborators.
* Internal units everywhere; convert input units in `configure`/`init` only.
* Honour `"verbose": false` (`--quiet`) for messages; write messages to the log stream the driver gives, not to stderr.
* ROOT objects: each worker owns its `TFile`/`TTree`/`TH1` (`SetDirectory(nullptr)`); call
  `ROOT::EnableThreadSafety()` in `makeWorkerSink`.
