/**
 * @file
 * @ingroup output
 * @brief The run-summary JSON ("simc_run_summary"): builder for the C++ run, parser for the Fortran .hist/.geni.
 *
 * The run-summary JSON ("simc_run_summary", docs/manual/markdown/reference/summary.schema.json, docs/manual/markdown/reference/OUTPUT_FORMAT.md):
 * one format for both codes so the comparison tool reads a single layout.
 *   * C++:     BuildSummaryJson(setup, summary)       (RunSummarySink)
 *   * Fortran: ParseFortranHist(.hist, .geni)          (tools/hist2json)
 * ROOT-free.
 */
#ifndef SIMC_OUTPUT_SUMMARYJSON_H
#define SIMC_OUTPUT_SUMMARYJSON_H

#include <map>
#include <string>
#include <vector>

#include "simc/core/Json.h"
#include "simc/core/RunSetup.h"
#include "simc/core/RunSummary.h"

namespace simc {

/// @addtogroup output
/// @{

/// Version of the summary schema written in the JSON.
constexpr int kSummarySchemaVersion = 1;
/// Name of the summary schema.
constexpr const char* kSummarySchemaName = "simc_run_summary";

/// Canonical stop-counter keys per spectrometer, in .geni print order.
/// Fortran: simc.f l.462-532 (.geni) and the *STOP_* counters of mc_hms.f, mc_shms.f, mc_sos.f,
/// mc_hrsr.f, mc_hrsl.f, mc_calo.f.
/// @param spectrometer spectrometer key ("hms", "sos", "hrsr", "hrsl", "shms", "calo").
/// @return the keys (empty for an unknown spectrometer).
const std::vector<std::string>& StopCounterKeys(const std::string& spectrometer);

/// Stop-counter key for an arm code (core/SpectrometerId.h StopKey: hms, sos, hrsr, hrsl, shms, calo; else "none").
/// @param arm_code arm code of the input.
/// @return the key.
std::string SpectrometerName(int arm_code);

/// C++ side: build the summary from the run setup and end-of-run numbers.
/// @param setup the run set-up.
/// @param summary the end-of-run numbers.
/// @return the summary JSON.
Json BuildSummaryJson(const RunSetup& setup, const RunSummary& summary);

/// Fortran side: parse a .hist (simc.f:report) and a .geni (simc.f l.462-532) into the same layout.
/// `geni_text` may be empty. Unparsed/unknown lines are reported in `warnings`.
/// @param hist_text the .hist text.
/// @param geni_text the .geni text (may be empty).
/// @param run_name the run name to store.
/// @param[out] warnings unparsed lines (optional, may be nullptr).
/// @return the summary JSON.
Json ParseFortranHist(const std::string& hist_text, const std::string& geni_text, const std::string& run_name,
                      std::vector<std::string>* warnings = nullptr);

/// Parse only a .geni text into {"<spectrometer>": {"<key>": n, ...}} (keys: StopCounterKeys).
/// @param geni_text the .geni text.
/// @param[out] warnings unparsed lines (optional, may be nullptr).
/// @return the stop counters per spectrometer.
Json ParseFortranGeni(const std::string& geni_text, std::vector<std::string>* warnings = nullptr);

/// Human-readable .hist-like rendering of a summary JSON.
/// @param summary the summary JSON.
/// @return the text.
std::string FormatSummaryText(const Json& summary);

/// Minimal JSON-Schema (draft-07 subset: type, required, properties, additionalProperties,
/// items, enum, minimum) validator used by the tests. Returns the list of violations.
/// @param instance the JSON to check.
/// @param schema the schema.
/// @return the violations (empty if valid).
std::vector<std::string> ValidateAgainstSchema(const Json& instance, const Json& schema);

/// @}

}  // namespace simc

#endif  // SIMC_OUTPUT_SUMMARYJSON_H
