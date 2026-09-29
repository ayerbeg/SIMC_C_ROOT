// include/simc/core/Json.h
// Single place that selects the JSON library (nlohmann/json) for the framework.
#ifndef SIMC_CORE_JSON_H
#define SIMC_CORE_JSON_H

#include <nlohmann/json.hpp>

namespace simc {
using Json = nlohmann::json;
}  // namespace simc

#endif  // SIMC_CORE_JSON_H
