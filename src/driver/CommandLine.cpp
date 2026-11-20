// src/driver/CommandLine.cpp — the `simc++` command line (a library function, CommandLineMain, so that a program that
// links extra plug-ins can reuse it: tools/simc_main.cpp, examples/out_of_tree_plugin).
//
//   simc++ <file.inp | file.json> [-o OUTDIR] [--seed N] [--ngen N] [--pipeline overrides.json]
//          [--run-name NAME] [--data-dir DIR] [--infiles-dir DIR] [--quiet] [--print-pipeline] [--fortran-compat] [--compat-target C]
//          [--threads N] [--print-data-dir] [--tree-layout upstream|legacy]
//
// Input:
//   * a Fortran .inp (CTP syntax, + its extra_dbase_file), read exactly like dbase.f (InpReader), or
//   * a SIMC++ JSON: {"parameters": {...} (inp2json output), "pipeline": {...} (optional, docs/manual/markdown/reference/DRIVER.md),
//                     "run": {"name": ..., "output_dir": ..., "fortran_compat": false} (optional)}.
// --fortran-compat (or "run": {"fortran_compat": true}) keeps the exact Fortran behaviour, including the clear Fortran
// bugs that SIMC++ corrects by default (D-031, docs/translation/markdown/FORTRAN_DIFFERENCES.md). It exists for the bitwise validation.
// --threads N (docs/manual/markdown/reference/PERFORMANCE.md): the event loop on N threads with independent mt64 streams (reproducible for a
// given seed and N; statistically equivalent to, not bitwise equal to, the one-thread run). Default 1 (or the pipeline
// "threads"). Refused (exit 2) with --fortran-compat and with an rng other than mt64.
// --tree-layout upstream|legacy (or JSON "output": {"layout": "upstream"|"legacy"}; task U.2, D-082): the column layout
// of the tree h10. upstream (default) = upstream SIMC 59ccdf6 (NtupleInit.f/results_write.f after e4dd175), legacy = SIMC
// 294f662 (docs/manual/markdown/reference/OUTPUT_FORMAT.md). The CLI option wins over the JSON.
// --print-data-dir prints the default data directory (include/simc/core/DataPath.h: env SIMC_DATA_DIR, the installed
// package's, <exe>/../share/simc/data, or the source tree) and exits.
// --pipeline FILE merges a pipeline JSON (either the section itself or {"pipeline": {...}}) over the input's.
// Outputs in OUTDIR (default "."): <run>.root (tree h10, if Nntu > 0), <run>.summary.json/.txt, <run>.gen.root,
// <run>.hist and <run>.geni (simc.f report layout). <run> defaults to the input file stem (Fortran `base`).
// Exit status: 0 ok, 1 usage, 2 input error, 3 run error (Fortran STOP, not-ported feature, I/O).
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "simc/core/DataPath.h"
#include "simc/core/Json.h"
#include "simc/core/SimcParameters.h"
#include "simc/driver/SimcRun.h"
#include "simc/input/InpReader.h"

#include "simc/driver/CommandLine.h"

namespace simc {
namespace driver {

namespace {

void usage() {
    std::cerr << "usage: simc++ <file.inp|file.json> [-o OUTDIR] [--seed N] [--ngen N] [--pipeline FILE]\n"
                 "               [--run-name NAME] [--data-dir DIR] [--infiles-dir DIR] [--quiet] [--print-pipeline]\n"
                 "               [--fortran-compat] [--compat-target 59ccdf6|294f662] [--threads N] [--print-data-dir]\n"
                 "               [--tree-layout upstream|legacy]\n"
                 "  --print-pipeline  print the resolved pipeline (JSON), then run\n"
                 "  --fortran-compat  exact Fortran behaviour (known Fortran bugs kept; bitwise validation)\n"
                 "  --compat-target C upstream Fortran commit reproduced by --fortran-compat: 59ccdf6 (default) or\n"
                 "                    294f662 (frozen tree: kaon decay daughter mass Mk, H which_pion=3 STOP)\n"
                 "  --threads N       event loop on N threads (independent mt64 streams per thread; a given seed and N\n"
                 "                    always give the same output; not with --fortran-compat or rng ranlux). Default 1.\n"
                 "  --tree-layout L   tree h10 column layout: upstream (default, upstream SIMC 59ccdf6) or legacy (SIMC 294f662)\n"
                 "  --print-data-dir  print the default data directory (env SIMC_DATA_DIR overrides it) and exit\n";
}

std::string stem(const std::string& path) {
    std::string b = path.substr(path.find_last_of('/') == std::string::npos ? 0 : path.find_last_of('/') + 1);
    const auto dot = b.find_last_of('.');
    return dot == std::string::npos ? b : b.substr(0, dot);
}

bool endsWith(const std::string& s, const std::string& e) {
    return s.size() >= e.size() && s.compare(s.size() - e.size(), e.size(), e) == 0;
}

// task U.2: "upstream" / "legacy" -> the layout; false for anything else
bool parseLayout(const std::string& v, TreeLayout& out) {
    if (v == "upstream") { out = TreeLayout::Upstream; return true; }
    if (v == "legacy") { out = TreeLayout::Legacy; return true; }
    return false;
}

// task U.R: "59ccdf6" / "294f662" -> the compat target; false for anything else
bool parseCompatTarget(const std::string& v, CompatTarget& out) {
    if (v == "59ccdf6") { out = CompatTarget::U59ccdf6; return true; }
    if (v == "294f662") { out = CompatTarget::U294f662; return true; }
    return false;
}

bool readJson(const std::string& path, Json& j) {
    std::ifstream f(path);
    if (!f) return false;
    try {
        f >> j;
    } catch (const std::exception& e) {
        std::cerr << "simc++: " << path << ": " << e.what() << "\n";
        return false;
    }
    return true;
}

}  // namespace

// The simc++ program: parse the options, read the input (.inp or JSON), run, write the outputs. Exit status 0 on
// success, 1 for a command-line error, 2 for an unreadable or rejected input/configuration, 3 when the run stops
// (a Fortran STOP, a refused feature or any other error during the run).
int CommandLineMain(int argc, char** argv) {
    std::string input, outdir = ".", pipelineFile, runName, dataDir, infilesDir;
    bool quiet = false, printPipeline = false, compatFlag = false, layoutFlag = false, targetFlag = false;
    CompatTarget targetCli = CompatTarget::U59ccdf6;
    TreeLayout layoutCli = TreeLayout::Upstream;
    driver::RunOptions opt;
    // --- 1. the options
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](std::string& v) {
            if (i + 1 >= argc) { usage(); std::exit(1); }
            v = argv[++i];
        };
        std::string v;
        if (a == "-o" || a == "--output-dir") next(outdir);
        else if (a == "--seed") { next(v); opt.override_seed = true; opt.seed = std::stol(v); }
        else if (a == "--ngen") { next(v); opt.override_ngen = true; opt.ngen = std::stol(v); }
        else if (a == "--pipeline") next(pipelineFile);
        else if (a == "--run-name") next(runName);
        else if (a == "--data-dir") next(dataDir);
        else if (a == "--infiles-dir") next(infilesDir);
        else if (a == "--quiet") quiet = true;
        else if (a == "--print-pipeline") printPipeline = true;
        else if (a == "--fortran-compat") compatFlag = true;
        else if (a == "--compat-target") {   // task U.R
            next(v);
            if (!parseCompatTarget(v, targetCli)) {
                std::cerr << "simc++: --compat-target: expected '59ccdf6' or '294f662', got '" << v << "'\n";
                return 1;
            }
            targetFlag = true;
        }
        else if (a == "--tree-layout") {
            next(v);
            if (!parseLayout(v, layoutCli)) {
                std::cerr << "simc++: --tree-layout: expected 'upstream' or 'legacy', got '" << v << "'\n";
                return 1;
            }
            layoutFlag = true;
        }
        else if (a == "--threads") {
            next(v);
            char* end = nullptr;
            const long n = std::strtol(v.c_str(), &end, 10);
            if (end == v.c_str() || *end != '\0' || n < 1 || n > 1024) {
                std::cerr << "simc++: --threads: expected an integer 1..1024, got '" << v << "'\n";
                return 1;
            }
            opt.threads = static_cast<int>(n);
        }
        else if (a == "--print-data-dir") { std::cout << DataRoot() << "\n"; return 0; }
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (!a.empty() && a[0] == '-') { usage(); return 1; }
        else if (input.empty()) input = a;
        else { usage(); return 1; }
    }
    if (input.empty()) { usage(); return 1; }

    // --- 2. the input: a SIMC++ JSON (parameters + optional pipeline/run sections) or a Fortran .inp
    SimcParameters params;

    std::vector<std::string> ignoredNames;   // InpReader::notRegistered()
    Json pipeline = Json::object();
    if (endsWith(input, ".json")) {
        Json j;
        if (!readJson(input, j)) { std::cerr << "simc++: cannot read " << input << "\n"; return 2; }
        std::vector<std::string> warnings;
        try {
            FromJson(j.contains("parameters") ? j["parameters"] : Json::object(), params, &warnings);
        } catch (const std::invalid_argument& e) {   // rejected input (D-032: arm code 6)
            std::cerr << "simc++: input error: " << input << ": " << e.what() << "\n";
            return 2;
        }
        for (const auto& w : warnings) std::cerr << "simc++: warning: " << w << "\n";
        if (j.contains("pipeline")) pipeline = j["pipeline"];
        if (j.contains("run")) {
            if (runName.empty()) runName = j["run"].value("name", std::string());
            if (outdir == ".") outdir = j["run"].value("output_dir", outdir);
            opt.fortran_compat = j["run"].value("fortran_compat", false);
            if (j["run"].contains("compat_target")) {   // task U.R
                const Json& t = j["run"]["compat_target"];
                if (!t.is_string() || !parseCompatTarget(t.get<std::string>(), opt.compat_target)) {
                    std::cerr << "simc++: input error: " << input
                              << ": \"run\": {\"compat_target\"} must be \"59ccdf6\" or \"294f662\"\n";
                    return 2;
                }
            }
        }
        if (j.contains("output")) {   // task U.2: "output": {"layout": "upstream" | "legacy"}
            const Json& o = j["output"];
            if (!o.is_object() || (o.contains("layout") && (!o["layout"].is_string() ||
                                                            !parseLayout(o["layout"].get<std::string>(), opt.tree_layout)))) {
                std::cerr << "simc++: input error: " << input
                          << ": \"output\": {\"layout\"} must be \"upstream\" or \"legacy\"\n";
                return 2;
            }
        }
    } else {
        InpReader::Options ro;
        ro.infiles_dir = infilesDir;
        InpReader rd(params, ro);
        if (!rd.loadRun(input)) {
            for (const auto& e : rd.errors()) std::cerr << "simc++: input error: " << input << ": " << e << "\n";
            if (rd.errors().empty()) {
                // Say which file failed: the main input or the extra_dbase_file it names (looked up next to the
                // main input, or in --infiles-dir).
                std::cerr << "simc++: cannot read " << input << "\n";
                for (const auto& m : rd.messages()) std::cerr << "simc++:   " << m.text << "\n";
                if (!params.extra_dbase_file.empty())
                    std::cerr << "simc++:   the input names extra_dbase_file='" << params.extra_dbase_file
                              << "': put that file next to the input, or give its folder with --infiles-dir DIR\n";
            }
            return 2;
        }
        for (const auto& m : rd.messages()) std::cout << m.text << "\n";   // CTP messages, like the Fortran
        ignoredNames = rd.notRegistered();
        if (!ignoredNames.empty()) {
            // D-065 W1: the Fortran ignores unregistered names silently (CTP registers them as dummy ints); a run
            // with legacy '.' names (targ.A, spec.e.P, ...) keeps the defaults of those variables. Same behaviour here,
            // but said loudly.
            std::cerr << "simc++: WARNING: " << ignoredNames.size()
                      << " input variable(s) are not SIMC variables and are IGNORED (legacy SIMC behaviour; their defaults are "
                         "used):\n  ";
            for (std::size_t i = 0; i < ignoredNames.size(); ++i)
                std::cerr << ignoredNames[i] << (i + 1 < ignoredNames.size() ? (i % 6 == 5 ? "\n  " : ", ") : "\n");
            std::cerr << "  Legacy '.' names (targ.A, spec.e.P, ...) must be converted to '%' names with "
                         "simc_gfortran/infiles/convert_inputfile.sh.\n";
        }
    }
    // --- 3. --pipeline overrides, merged key by key over the input's pipeline section
    if (!pipelineFile.empty()) {
        Json j;
        if (!readJson(pipelineFile, j)) { std::cerr << "simc++: cannot read " << pipelineFile << "\n"; return 2; }
        const Json& sec = j.contains("pipeline") ? j["pipeline"] : j;
        for (auto it = sec.begin(); it != sec.end(); ++it) pipeline[it.key()] = it.value();
    }
    // --- 4. run options (the run name defaults to the input file name without extension, as the Fortran)
    opt.run_name = runName.empty() ? stem(input) : runName;
    opt.output_dir = outdir;
    opt.data_dir = dataDir;
    opt.quiet_stages = quiet;
    if (compatFlag) opt.fortran_compat = true;   // the CLI flag wins over the JSON
    if (targetFlag) opt.compat_target = targetCli;   // task U.R
    if (layoutFlag) opt.tree_layout = layoutCli;   // idem (task U.2)
    if (opt.threads > 1 && opt.fortran_compat) {  // refused before the run starts (also SimcRun::threadsRefusal)
        std::cerr << "simc++: --threads " << opt.threads
                  << " (more than one event-loop thread) cannot be used with --fortran-compat: the exact Fortran "
                     "reproduction needs the single random stream in the Fortran order (run with --threads 1)\n";
        return 2;
    }

    // --- 5. the run
    try {
        std::filesystem::create_directories(outdir);   // like the Fortran outfiles/, but created on demand
        driver::SimcRun run(params, pipeline, opt);
        run.initialise();
        if (printPipeline) std::cout << run.pipeline().toJson().dump(2) << "\n";
        const int nthreads = run.threadCount();
        if (nthreads > 1) {   // the parallel event loop (refused with a message for compat / ranlux)
            if (const std::string why = run.threadsRefusal(nthreads); !why.empty()) {
                std::cerr << "simc++: " << why << "\n";
                return 2;
            }
            run.runThreaded(nthreads);
        } else {
            while (!run.done()) run.step();
            run.finish();
        }
    } catch (const std::invalid_argument& e) {   // a configuration rejected at start-up (e.g. D-037 F9)
        std::cerr << "simc++: configuration error: " << e.what() << "\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "simc++: " << e.what() << "\n";
        if (!ignoredNames.empty())
            std::cerr << "simc++: note: the input has " << ignoredNames.size()
                      << " ignored variable(s) (listed above), so their defaults were used — the likely cause of this "
                         "STOP. Convert legacy '.' names with simc_gfortran/infiles/convert_inputfile.sh.\n";
        return 3;
    }
    return 0;
}

}  // namespace driver
}  // namespace simc
