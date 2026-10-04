// src/driver/Pipeline.cpp — see include/simc/driver/Pipeline.h.
#include "simc/driver/Pipeline.h"

#include <stdexcept>

#include "simc/driver/MonteCarlo.h"

namespace simc {
namespace driver {

void StageSpec::apply(const Json& j, const std::string& what) {
    if (j.is_string()) {
        type = j.get<std::string>();
        return;
    }
    if (!j.is_object()) throw std::invalid_argument("pipeline." + what + ": expected a name or an object");
    for (auto it = j.begin(); it != j.end(); ++it) {
        if (it.key() == "type")
            type = it.value().get<std::string>();
        else
            config[it.key()] = it.value();
    }
}

Json StageSpec::toJson() const {
    Json j = config.is_object() ? config : Json::object();
    j["type"] = type;
    return j;
}

namespace {
constexpr const char* kDefaultTarget = "hallc_cryo";   // the Fortran's target model (target.f)
}  // namespace

std::string PipelineConfig::TargetType(const Json& section) {
    StageSpec t;
    t.type = kDefaultTarget;
    if (section.is_object() && section.contains("target")) t.apply(section.at("target"), "target");
    return t.type;
}

PipelineConfig PipelineConfig::Defaults(const RunSetup& s) {
    PipelineConfig p;
    const ReactionFlags& f = s.flags;
    p.rng.type = "mt64";
    if (s.params.random_seed > 0) p.rng.config["seed"] = static_cast<std::uint64_t>(s.params.random_seed);
    p.target.type = kDefaultTarget;
    p.target.config["Mh"] = s.Mh;
    p.beam.type = "standard";
    // event.f dispatch: the reaction is chosen by the flags. FORTRAN-OPEN(D-031 flag 1, D-009): H + doing_phsp keeps
    // doing_hyd_elast, so complete_ev runs the elastic kinematics (event.f:546 before 709).
    if (f.doing_phsp)
        p.reaction.type = "phase_space";
    else if (f.doing_hyd_elast)
        p.reaction.type = "hydrogen_elastic";
    else
        p.reaction.type = "";   // not ported: SimcRun reports NotPorted with the reaction title
    // event.f:complete_main l.1444-1456
    if (f.doing_deuterium || f.doing_heavy) {
        p.physics.type = "deforest";
        p.physics.config["Mh"] = s.Mh;
    } else {
        p.physics.type = "sigep";
    }
    p.radiator.type = "simc_rad";
    p.hadron_arm.type = "";
    p.electron_arm.type = "";
    if (s.params.Nntu > 0) p.sinks.push_back({"root_tree", Json::object()});
    p.sinks.push_back({"run_summary", Json::object()});
    p.sinks.push_back({"gen_hist", Json::object()});
    p.sinks.push_back({"fortran_hist", Json::object()});
    return p;
}

void PipelineConfig::apply(const Json& j) {
    if (j.is_null()) return;
    if (!j.is_object()) throw std::invalid_argument("pipeline: expected an object");
    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string& k = it.key();
        if (!k.empty() && k[0] == '_') continue;   // comments
        if (k == "rng") rng.apply(it.value(), k);
        else if (k == "target") target.apply(it.value(), k);
        else if (k == "beam") beam.apply(it.value(), k);
        else if (k == "reaction") reaction.apply(it.value(), k);
        else if (k == "physics") physics.apply(it.value(), k);
        else if (k == "radiator") radiator.apply(it.value(), k);
        else if (k == "spectrometers") {
            if (!it.value().is_object()) throw std::invalid_argument("pipeline.spectrometers: expected an object");
            for (auto a = it.value().begin(); a != it.value().end(); ++a) {
                if (a.key() == "hadron") hadron_arm.apply(a.value(), "spectrometers.hadron");
                else if (a.key() == "electron") electron_arm.apply(a.value(), "spectrometers.electron");
                else throw std::invalid_argument("pipeline.spectrometers: unknown key '" + a.key() + "'");
            }
        } else if (k == "sinks") {
            if (!it.value().is_array()) throw std::invalid_argument("pipeline.sinks: expected an array");
            sinks.clear();
            for (const auto& e : it.value()) {
                StageSpec s;
                s.apply(e, "sinks[]");
                sinks.push_back(s);
            }
        } else {
            throw std::invalid_argument("pipeline: unknown stage '" + k + "'");
        }
    }
}

Json PipelineConfig::toJson() const {
    Json j;
    j["rng"] = rng.toJson();
    j["target"] = target.toJson();
    j["beam"] = beam.toJson();
    j["reaction"] = reaction.toJson();
    j["physics"] = physics.toJson();
    j["radiator"] = radiator.toJson();
    j["spectrometers"] = {{"hadron", hadron_arm.toJson()}, {"electron", electron_arm.toJson()}};
    j["sinks"] = Json::array();
    for (const auto& s : sinks) j["sinks"].push_back(s.toJson());
    return j;
}

}  // namespace driver
}  // namespace simc
