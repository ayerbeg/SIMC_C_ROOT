// examples/dev_manual/FlatModel.cpp — the IPhysicsModel example of the developer manual (task P15.D1,
// docs/manual/markdown/developer/04_extending.md). "example_flat": a constant cross section, JSON {"type": "example_flat",
// "value": 2.5}. The manual quotes this file verbatim; tests/docs checks that the quotes match.
#include <string>

#include "simc/core/Registry.h"
#include "simc/stages/IPhysicsModel.h"

namespace mylab {

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
