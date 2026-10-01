// src/output/RunSummarySink.cpp
// Fortran: simc.f:report (l.644-1142) and .geni (l.462-532); see docs/OUTPUT_FORMAT.md.
#include "simc/output/RunSummarySink.h"

#include <fstream>
#include <stdexcept>

#include "simc/core/Registry.h"
#include "simc/core/RunSummary.h"
#include "simc/output/SummaryJson.h"

namespace simc {

void RunSummarySink::configure(const Json& config, const SimcParameters&) {
    if (config.is_object()) file_ = config.value("file", file_);
}

void RunSummarySink::init(const StageContext&) {}

void RunSummarySink::beginRun(const RunSetup& setup) {
    setup_ = setup;
    if (file_.empty()) file_ = setup.output_dir + "/" + setup.run_name + ".summary.json";
}

void RunSummarySink::endRun(const RunSummary& summary) {
    last_ = BuildSummaryJson(setup_, summary);
    std::ofstream f(file_);
    if (!f) throw std::runtime_error("RunSummarySink: cannot write " + file_);
    f << last_.dump(2) << "\n";
    std::string txt = file_;
    const std::string ext = ".json";
    if (txt.size() > ext.size() && txt.compare(txt.size() - ext.size(), ext.size(), ext) == 0)
        txt.replace(txt.size() - ext.size(), ext.size(), ".txt");
    else
        txt += ".txt";
    std::ofstream t(txt);
    t << FormatSummaryText(last_);
}

}  // namespace simc

SIMC_REGISTER(simc::IOutputSink, "run_summary", simc::RunSummarySink)
