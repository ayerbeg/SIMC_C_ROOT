// tests/unit/test_stage_headers.cpp — every stage interface compiles, can be implemented,
// and self-registers in the name-keyed Registry (docs/ARCHITECTURE.md §3).
#include <stdexcept>

#include "TestUtil.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcParameters.h"
#include "simc/stages/IBeam.h"
#include "simc/stages/IMaterialModel.h"
#include "simc/stages/IOutputSink.h"
#include "simc/stages/IPhysicsModel.h"
#include "simc/stages/IRadiator.h"
#include "simc/stages/IRandom.h"
#include "simc/stages/IReaction.h"
#include "simc/stages/ISpectrometer.h"
#include "simc/stages/ITarget.h"

using namespace simc;
using simctest::check;

namespace {

#define DUMMY_STAGE_BASICS(N)                                            \
    std::string name() const override { return N; }                      \
    void configure(const Json&, const SimcParameters&) override {}       \
    void init(RunSetup&) override {}

struct DummyRandom : IRandom {
    DUMMY_STAGE_BASICS("dummy")
    double uniform() override { return 0.5; }
    double gauss1(double) override { return 0.0; }
    void seed(std::uint64_t) override {}
};
struct DummyBeam : IBeam {
    DUMMY_STAGE_BASICS("dummy")
    void generate(MainEvent&, Event&, IRandom&) override {}
};
struct DummyTarget : ITarget {
    DUMMY_STAGE_BASICS("dummy")
    double generateZ(IRandom& r) override { return r.uniform(); }
    void tripThruTarget(int, double, double, double, double, LossMode, IRandom&, double& e, double& t) override { e = t = 0; }
    void multipleScattering(double, double, double, double d[2], IRandom&) override { d[0] = d[1] = 0; }
};
struct DummyReaction : IReaction {
    DUMMY_STAGE_BASICS("dummy")
    bool generate(MainEvent&, Event&, Event&, IRandom&) override { return true; }
    bool completeEv(MainEvent&, Event&, IRandom&) override { return true; }
    bool completeRecon(Event&) override { return true; }
    double genVolume() const override { return 1.0; }
};
struct DummyModel : IPhysicsModel {
    DUMMY_STAGE_BASICS("dummy")
    double sigcc(const Event&, const MainEvent&) override { return 1.0; }
    double sigccRecon(const Event&, const MainEvent&) override { return 1.0; }
};
struct DummyRadiator : IRadiator {
    DUMMY_STAGE_BASICS("dummy")
    void initEvent(MainEvent&, const Event&) override {}
    bool radiate(MainEvent&, Event&, Event&, IRandom&, const CompleteEvFn&) override { return true; }
};
struct DummyMaterial : IMaterialModel {
    DUMMY_STAGE_BASICS("dummy")
    double energyLoss(double, double, double, double, double, double, int, IRandom&) override { return 0; }
    void musc(double, double, double, double&, double&, IRandom&) override {}
    void muscExt(double, double, double, double, double&, double&, double&, double&, IRandom&) override {}
};
struct DummySpec : ISpectrometer {
    DUMMY_STAGE_BASICS("dummy")
    SpecResult transport(const SpecInput&, IRandom&) override {
        ++counters_["trials"];
        SpecResult r;
        r.stop_reason = "slit_hor";
        ++counters_[r.stop_reason];
        return r;
    }
    const std::map<std::string, long>& stopCounters() const override { return counters_; }
    std::map<std::string, long> counters_;
};
struct DummySink : IOutputSink {
    DUMMY_STAGE_BASICS("dummy")
    void beginRun(const RunSetup&) override {}
    void fill(const MainEvent&, const Event&, const Event&, const Event&, bool) override {}
    void endRun(const RunSummary&) override {}
};

}  // namespace

SIMC_REGISTER(IRandom, "dummy", DummyRandom)
SIMC_REGISTER(IBeam, "dummy", DummyBeam)
SIMC_REGISTER(ITarget, "dummy", DummyTarget)
SIMC_REGISTER(IReaction, "dummy", DummyReaction)
SIMC_REGISTER(IPhysicsModel, "dummy", DummyModel)
SIMC_REGISTER(IRadiator, "dummy", DummyRadiator)
SIMC_REGISTER(IMaterialModel, "dummy", DummyMaterial)
SIMC_REGISTER(ISpectrometer, "dummy", DummySpec)
SIMC_REGISTER(IOutputSink, "dummy", DummySink)

int main() {
    auto rng = Registry<IRandom>::create("dummy");
    check(rng && rng->name() == "dummy", "IRandom created by name");
    check(Registry<IBeam>::contains("dummy"), "IBeam registered");
    check(Registry<ITarget>::create("dummy")->generateZ(*rng) == 0.5, "ITarget call");
    check(Registry<IReaction>::create("dummy")->genVolume() == 1.0, "IReaction call");
    check(Registry<IPhysicsModel>::contains("dummy"), "IPhysicsModel registered");
    check(Registry<IRadiator>::contains("dummy"), "IRadiator registered");
    check(Registry<IMaterialModel>::contains("dummy"), "IMaterialModel registered");
    check(Registry<IOutputSink>::contains("dummy"), "IOutputSink registered");
    auto spec = Registry<ISpectrometer>::create("dummy");
    SpecInput in;
    auto res = spec->transport(in, *rng);
    check(!res.ok && res.stop_reason == "slit_hor" && spec->stopCounters().at("trials") == 1, "ISpectrometer stop counters");
    bool threw = false;
    try { Registry<ISpectrometer>::create("nope"); } catch (const std::out_of_range&) { threw = true; }
    check(threw, "unknown plug-in name throws");
    threw = false;
    try { Registry<IBeam>::add("dummy", [] { return std::unique_ptr<IBeam>(new DummyBeam()); }); }
    catch (const std::logic_error&) { threw = true; }
    check(threw, "duplicate plug-in name throws");
    return simctest::finish("stage_headers_and_registry");
}
