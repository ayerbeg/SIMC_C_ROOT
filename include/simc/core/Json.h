/**
 * @file
 * @ingroup core
 * @brief The JSON type of the framework (nlohmann/json), selected in one place.
 */
#ifndef SIMC_CORE_JSON_H
#define SIMC_CORE_JSON_H

#include <nlohmann/json.hpp>

namespace simc {
/**
 * @ingroup core
 * @brief The JSON value type used for configurations (stage `configure()`), the pipeline and the run summary.
 */
using Json = nlohmann::json;
}  // namespace simc

#endif  // SIMC_CORE_JSON_H
