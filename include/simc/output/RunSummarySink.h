// include/simc/output/RunSummarySink.h
//
// Writes <run>.summary.json (schema docs/summary.schema.json) and a human-readable
// <run>.summary.txt at the end of the run. The numbers are those of simc.f:report (.hist)
// and the .geni stop counters.
//
// JSON config (optional): {"file": "<path>.summary.json"}; default <output_dir>/<run_name>.summary.json
#ifndef SIMC_OUTPUT_RUNSUMMARYSINK_H
#define SIMC_OUTPUT_RUNSUMMARYSINK_H

#include <string>

#include "simc/core/RunSetup.h"
#include "simc/stages/IOutputSink.h"

namespace simc {

class RunSummarySink : public IOutputSink {
public:
    std::string name() const override { return "run_summary"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    void beginRun(const RunSetup& setup) override;
    void fill(const MainEvent&, const Event&, const Event&, const Event&, bool) override {}
    void endRun(const RunSummary& summary) override;

    const std::string& fileName() const { return file_; }
    const Json& written() const { return last_; }

private:
    std::string file_;
    RunSetup setup_;
    Json last_;
};

}  // namespace simc

#endif
