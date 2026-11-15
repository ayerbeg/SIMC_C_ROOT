// cmake/SimcInstalledDataRoot.cpp — installed as <prefix>/share/simc/cmake/; an INTERFACE source of SIMC::simc set by
// SIMCConfig.cmake: compiled into every program built with find_package(SIMC), it registers the package's data directory
// (SIMC_INSTALLED_DATA_DIR = <prefix>/share/simc/data, defined by SIMCConfig.cmake) as the default data root
// (include/simc/core/DataPath.h, task P15.4).
#include "simc/core/DataPath.h"

namespace {
const bool simc_installed_data_root_registered = simc::detail::SetInstalledDataRoot(SIMC_INSTALLED_DATA_DIR);
}  // namespace
