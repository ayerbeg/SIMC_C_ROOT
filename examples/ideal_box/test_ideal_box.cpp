// examples/ideal_box/test_ideal_box.cpp — the complete driver with the example plug-in "ideal_box" (task P7.4).
//
//   test_example_ideal_box <file.inp> <pipeline.json> <outdir>
//
// Runs SimcRun (2000 successes, seed 4242) with the "pipeline" section of <pipeline.json>, which selects ideal_box for
// one or both arms, and checks that the output is sane:
//   * the resolved pipeline uses ideal_box for the arm(s) the JSON names;
//   * the run completes: nevent = ngen, successes > 0, every weight finite;
//   * the stop counters add up: trials = successes + delta + xptar + yptar, successes >= ncontribute;
//   * perfect resolution: for every success RECON = SP (minus the spectrometer offsets) BITWISE, inside the box;
//   * the output files exist, the .geni has the plug-in's counter block and the summary JSON its counters.
// Exit 0 = all checks pass, 1 = a check failed, 77 = input missing (CTest SKIP).
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "simc/core/Json.h"
#include "simc/core/SimcParameters.h"
#include "simc/driver/SimcRun.h"
#include "simc/input/InpReader.h"
#include "simc/stages/IOutputSink.h"

using namespace simc;

namespace {

int failures = 0;
void expect(bool ok, const std::string& what) {
    std::printf("  [%s] %s\n", ok ? " ok " : "FAIL", what.c_str());
    if (!ok) ++failures;
}

/// Records SP, RECON and the weight of every success.
struct Capture : IOutputSink {
    struct Row {
        ArmSP SPp, SPe, REp, REe;
        double weight;
    };
    std::vector<Row> rows;
    std::string name() const override { return "capture"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}
    void beginRun(const RunSetup&) override {}
    void fill(const MainEvent& m, const Event&, const Event&, const Event&, bool success) override {
        if (success) rows.push_back({m.SP.p, m.SP.e, m.RECON.p, m.RECON.e, m.weight});
    }
    void endRun(const RunSummary&) override {}
};

std::string slurp(const std::string& path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

/// Perfect-resolution and box checks for one arm.
void checkArm(const char* arm, const Json& box, const SpecArm& spec, const std::vector<Capture::Row>& rows,
              bool hadron) {
    long bad = 0, outside = 0;
    for (const auto& r : rows) {
        const ArmSP& sp = hadron ? r.SPp : r.SPe;
        const ArmSP& re = hadron ? r.REp : r.REe;
        const double dx = sp.xptar - spec.offset.xptar, dy = sp.yptar - spec.offset.yptar;
        if (re.delta != sp.delta || re.xptar != dx || re.yptar != dy) ++bad;
        if (re.delta < box.value("delta_min", -15.) || re.delta > box.value("delta_max", 15.) ||
            re.xptar < box.value("xptar_min", -0.1) || re.xptar > box.value("xptar_max", 0.1) ||
            re.yptar < box.value("yptar_min", -0.1) || re.yptar > box.value("yptar_max", 0.1))
            ++outside;
    }
    expect(bad == 0, std::string(arm) + " arm: RECON == SP - offset bitwise for every success (perfect resolution), " +
                         std::to_string(bad) + " mismatches");
    expect(outside == 0, std::string(arm) + " arm: every success inside the box (" + std::to_string(outside) + " outside)");
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: test_example_ideal_box <file.inp> <pipeline.json> <outdir>\n");
        return 2;
    }
    const std::string inp = argv[1], pipeFile = argv[2], outdir = argv[3];
    SimcParameters params;
    InpReader rd(params);
    if (!rd.loadRun(inp)) {
        std::printf("SKIP: cannot read %s\n", inp.c_str());
        return 77;
    }
    Json j;
    {
        std::ifstream f(pipeFile);
        if (!f) {
            std::printf("cannot read %s\n", pipeFile.c_str());
            return 2;
        }
        f >> j;
    }
    const Json pipeline = j.at("pipeline");
    const Json& specs = pipeline.at("spectrometers");
    const bool hadronBox = specs.contains("hadron"), electronBox = specs.contains("electron");

    std::filesystem::create_directories(outdir);
    driver::RunOptions opt;
    opt.run_name = "ideal_box";
    opt.output_dir = outdir;
    opt.override_ngen = true;
    opt.ngen = 2000;
    opt.override_seed = true;
    opt.seed = 4242;
    opt.quiet_stages = true;
    std::ostringstream log;
    opt.log = &log;

    driver::SimcRun run(params, pipeline, opt);
    auto capOwned = std::make_unique<Capture>();
    Capture* cap = capOwned.get();
    run.addSink(std::move(capOwned));
    try {
        run.run();
    } catch (const std::exception& e) {
        std::printf("run error: %s\n", e.what());
        return 1;
    }
    std::printf("ideal_box: %s\n", pipeFile.c_str());

    // pipeline resolution
    const Json resolved = run.pipeline().toJson();
    if (hadronBox) expect(resolved["spectrometers"]["hadron"]["type"] == "ideal_box", "hadron arm resolved to ideal_box");
    if (electronBox)
        expect(resolved["spectrometers"]["electron"]["type"] == "ideal_box", "electron arm resolved to ideal_box");
    expect(run.montecarlo().hadronArm().name() == (hadronBox ? "ideal_box" : "hms"), "hadron arm object");

    // run totals
    const auto& T = run.totals();
    expect(T.nevent == 2000 && T.ncontribute == 2000, "2000 successes generated (nevent " + std::to_string(T.nevent) + ")");
    bool finite = true;
    for (const auto& r : cap->rows) finite = finite && std::isfinite(r.weight);
    expect(!cap->rows.empty() && finite, "every success has a finite weight");
    expect(run.report().norm.normfac > 0 && std::isfinite(run.report().norm.normfac), "normfac finite and positive");

    // stop counters of the plug-in(s)
    for (const auto& sc : run.report().stops) {
        if (sc.spectrometer.rfind("ideal_box", 0) != 0) continue;
        const auto& c = sc.counts;
        const long sum = c.at("successes") + c.at("delta") + c.at("xptar") + c.at("yptar");
        expect(c.at("trials") == sum, sc.spectrometer + ": trials (" + std::to_string(c.at("trials")) +
                                          ") = successes + delta + xptar + yptar");
        expect(c.at("successes") >= T.ncontribute && c.at("delta") + c.at("xptar") + c.at("yptar") > 0,
               sc.spectrometer + ": successes >= ncontribute and the box cuts something");
    }

    // perfect resolution
    if (hadronBox) checkArm("hadron", specs["hadron"], run.setup().spec.p, cap->rows, true);
    if (electronBox) checkArm("electron", specs["electron"], run.setup().spec.e, cap->rows, false);

    // files
    const std::string base = outdir + "/ideal_box";
    const std::string geni = slurp(base + ".geni"), summary = slurp(base + ".summary.json");
    expect(!slurp(base + ".hist").empty() && !slurp(base + ".gen").empty(), ".hist and .gen written");
    expect(geni.find("ideal_box") != std::string::npos && geni.find("stop counters") != std::string::npos,
           ".geni has the plug-in's counter block");
    if (params.Nntu > 0) expect(!slurp(base + ".root").empty(), ".root written (Nntu > 0)");
    Json sj = Json::parse(summary);
    bool inSummary = false;
    for (auto it = sj["stops"].begin(); it != sj["stops"].end(); ++it)
        inSummary = inSummary || it.key().rfind("ideal_box", 0) == 0;
    expect(inSummary, "summary JSON has the ideal_box counters");

    std::printf("%s (%d failure(s))\n", failures ? "FAIL" : "PASS", failures);
    return failures ? 1 : 0;
}
