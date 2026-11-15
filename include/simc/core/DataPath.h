/**
 * @file
 * @ingroup core
 * @brief Where SIMC++ finds its data files by default (COSY maps, model tables, spectral functions, field map).
 *
 * The data directory holds `matrices/` (COSY forward/reconstruction maps), `pion/`, `kaon/`, `sidis/` (model tables),
 * `spectral/` (nuclear spectral functions) and `field/` (the polarised-target field map). simc::DataRoot() returns the
 * first of these that applies:
 *   1. the environment variable `SIMC_DATA_DIR` (if set and not empty);
 *   2. for a program built against the INSTALLED package (`find_package(SIMC)`): `<prefix>/share/simc/data` of that
 *      package, registered by the package's interface source `SimcInstalledDataRoot.cpp` through
 *      simc::detail::SetInstalledDataRoot() (relocatable: the prefix is the one the package was found in);
 *   3. for an installed program (`<prefix>/bin/simc++`): `<exe dir>/../share/simc/data`, if that directory exists;
 *   4. the source tree the libraries were built from (`<source>/data`): in-tree builds and tests.
 *
 * Each stage can still be pointed elsewhere (`--data-dir` for the spectral/field files, the `"data_dir"` /
 * `"matrix_dir"` JSON keys of the stages). Task P15.4 (decision D-069).
 */
#ifndef SIMC_CORE_DATAPATH_H
#define SIMC_CORE_DATAPATH_H

#include <string>

namespace simc {

/** @addtogroup core
 *  @{ */

/**
 * @brief The default data directory of SIMC++ (resolution order in the file description).
 * @return Absolute or relative path of the data directory, without a trailing slash. The directory is not checked for
 *         existence in steps 1, 2 and 4; a missing file is reported by the stage that opens it.
 * @note Thread-safe after static initialisation (it only reads the environment and a value set during static
 *       initialisation).
 */
std::string DataRoot();

/**
 * @brief A sub-directory of the data directory.
 * @param sub Relative path below the data directory, e.g. `"pion"` or `"matrices/simc_fortran"`.
 * @return `DataRoot() + "/" + sub`.
 */
std::string DataDir(const std::string& sub);

/// @brief Implementation details of the data-directory lookup (not for use by applications).
namespace detail {
/**
 * @brief Registers the data directory of an installed SIMC package (step 2 of the lookup).
 *
 * Called during static initialisation by the interface source `SimcInstalledDataRoot.cpp` that `find_package(SIMC)`
 * adds to the consumer's executable. Applications do not call it.
 * @param dir Absolute path of `<prefix>/share/simc/data` of the installed package (copied).
 * @return Always true (so that it can initialise a static variable).
 */
bool SetInstalledDataRoot(const char* dir);
}  // namespace detail

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_DATAPATH_H
