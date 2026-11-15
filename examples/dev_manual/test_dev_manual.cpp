// examples/dev_manual/test_dev_manual.cpp — runs the developer-manual example plug-ins through the complete driver
// (task P15.D1, docs/manual/developer/08_extending.md):
//
//   test_dev_manual <file.inp> <outdir>
//
// 1. one thread: pipeline {"physics": {"type": "example_flat", "value": 2.5}, "sinks": ["run_summary",
//    "example_count"]}, 300 successes, seed 4242: every success has sigcc = 2.5 (checked by a sink added in code with
//    SimcRun::addSink), the example_count file agrees with the run totals (tries = ntried, successes = ncontribute);
// 2. two threads (pipeline "threads": 2): the example_count file of the merged run has successes = ngen and
//    tries = the merged ntried.
// Exit 0 = pass, 1 = a check failed, 77 = input missing (CTest SKIP).
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

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

struct SigccCheck : IOutputSink {
    long successes{0}, bad{0};
    std::string name() const override { return "sigcc_check"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}
    void beginRun(const RunSetup&) override {}
    void fill(const MainEvent& m, const Event&, const Event&, const Event&, bool success) override {
        if (!success) return;
        ++successes;
        if (m.sigcc != 2.5) ++bad;
    }
    void endRun(const RunSummary&) override {}
};

std::map<std::string, double> readCount(const std::string& path) {
    std::map<std::string, double> m;
    std::ifstream f(path);
    std::string k;
    double v;
    while (f >> k >> v) m[k] = v;
    return m;
}

driver::RunOptions options(const std::string& outdir, const std::string& name) {
    driver::RunOptions opt;
    opt.run_name = name;
    opt.output_dir = outdir;
    opt.override_ngen = true;
    opt.ngen = 300;
    opt.override_seed = true;
    opt.seed = 4242;
    opt.quiet_stages = true;
    return opt;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: test_dev_manual <file.inp> <outdir>\n");
        return 2;
    }
    const std::string inp = argv[1], outdir = argv[2];
    SimcParameters params;
    InpReader rd(params);
    if (!rd.loadRun(inp)) {
        std::printf("SKIP: cannot read %s\n", inp.c_str());
        return 77;
    }
    std::filesystem::create_directories(outdir);
    const Json pipeline = Json::parse(R"({"physics": {"type": "example_flat", "value": 2.5},
                                          "sinks": ["run_summary", "example_count"]})");
    std::ostringstream log;

    std::printf("1. one thread\n");
    {
        driver::RunOptions opt = options(outdir, "flat_t1");
        opt.log = &log;
        driver::SimcRun run(params, pipeline, opt);
        auto& chk = static_cast<SigccCheck&>(run.addSink(std::make_unique<SigccCheck>()));
        run.run();
        expect(run.pipeline().physics.type == "example_flat", "the pipeline uses the physics model example_flat");
        expect(chk.successes == run.totals().ncontribute && chk.successes > 0,
               "sigcc_check saw every success (" + std::to_string(chk.successes) + ")");
        expect(chk.bad == 0, "every success has sigcc = 2.5 (" + std::to_string(chk.bad) + " others)");
        const auto c = readCount(outdir + "/flat_t1.count.txt");
        expect(c.count("tries") && static_cast<long>(c.at("tries")) == run.totals().ntried,
               "example_count tries = ntried (" + std::to_string(run.totals().ntried) + ")");
        expect(c.count("successes") && static_cast<long>(c.at("successes")) == run.totals().ncontribute,
               "example_count successes = ncontribute");
    }

    std::printf("2. two threads\n");
    {
        Json p2 = pipeline;
        p2["threads"] = 2;
        driver::RunOptions opt = options(outdir, "flat_t2");
        opt.log = &log;
        driver::SimcRun run(params, p2, opt);
        run.run();
        expect(run.threadCount() == 2, "the run used 2 threads");
        const auto c = readCount(outdir + "/flat_t2.count.txt");
        expect(c.count("successes") && static_cast<long>(c.at("successes")) == 300,
               "merged example_count successes = ngen = 300");
        expect(c.count("tries") && static_cast<long>(c.at("tries")) == run.totals().ntried &&
                   static_cast<long>(c.at("ntried")) == run.totals().ntried,
               "merged example_count tries = merged ntried (" + std::to_string(run.totals().ntried) + ")");
    }
    std::printf("%s (%d failures)\n", failures ? "FAIL" : "PASS", failures);
    return failures ? 1 : 0;
}
