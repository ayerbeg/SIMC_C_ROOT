// include/simc/output/SummaryJson.h
//
// The run-summary JSON ("simc_run_summary", docs/summary.schema.json, docs/OUTPUT_FORMAT.md):
// one format for both codes so the comparison tool reads a single layout.
//   * C++:     BuildSummaryJson(setup, summary)       (RunSummarySink)
//   * Fortran: ParseFortranHist(.hist, .geni)          (tools/hist2json)
// ROOT-free.
#ifndef SIMC_OUTPUT_SUMMARYJSON_H
#define SIMC_OUTPUT_SUMMARYJSON_H

#include <map>
#include <string>
#include <vector>

#include "simc/core/Json.h"
#include "simc/core/RunSetup.h"
#include "simc/core/RunSummary.h"

namespace simc {

constexpr int kSummarySchemaVersion = 1;
constexpr const char* kSummarySchemaName = "simc_run_summary";

/// Canonical stop-counter keys per spectrometer, in .geni print order.
/// Fortran: simc.f l.462-532 (.geni) and the *STOP_* counters of mc_hms.f, mc_shms.f, mc_sos.f,
/// mc_hrsr.f, mc_hrsl.f, mc_calo.f.
const std::vector<std::string>& StopCounterKeys(const std::string& spectrometer);

/// Spectrometer name for a Fortran arm code (1 hms, 2 sos, 3 hrsr, 4 hrsl, 5/6 shms, 7/8 calo).
std::string SpectrometerName(int arm_code);

/// C++ side: build the summary from the run setup and end-of-run numbers.
Json BuildSummaryJson(const RunSetup& setup, const RunSummary& summary);

/// Fortran side: parse a .hist (simc.f:report) and a .geni (simc.f l.462-532) into the same layout.
/// `geni_text` may be empty. Unparsed/unknown lines are reported in `warnings`.
Json ParseFortranHist(const std::string& hist_text, const std::string& geni_text, const std::string& run_name,
                      std::vector<std::string>* warnings = nullptr);

/// Parse only a .geni text into {"<spectrometer>": {"<key>": n, ...}} (keys: StopCounterKeys).
Json ParseFortranGeni(const std::string& geni_text, std::vector<std::string>* warnings = nullptr);

/// Human-readable .hist-like rendering of a summary JSON.
std::string FormatSummaryText(const Json& summary);

/// Minimal JSON-Schema (draft-07 subset: type, required, properties, additionalProperties,
/// items, enum, minimum) validator used by the tests. Returns the list of violations.
std::vector<std::string> ValidateAgainstSchema(const Json& instance, const Json& schema);

}  // namespace simc

#endif  // SIMC_OUTPUT_SUMMARYJSON_H
