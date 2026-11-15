// include/simc/core/DataPath.h
//
// Where SIMC++ finds its data files (COSY maps data/matrices, model tables data/pion|kaon|sidis, spectral functions
// data/spectral, the target field map data/field) by default — task P15.4 (D-069). The first that applies:
//   1. the environment variable SIMC_DATA_DIR (if set and not empty);
//   2. a program built against the INSTALLED package (find_package(SIMC)): <prefix>/share/simc/data of that package
//      (registered by the package's interface source SimcInstalledDataRoot.cpp; relocatable: the prefix is the one the
//      package was found in);
//   3. an installed program (<prefix>/bin/simc++): <exe dir>/../share/simc/data, if that directory exists;
//   4. the source tree the libraries were built from (<source>/data): in-tree builds and tests.
// Each stage can still be pointed elsewhere (--data-dir for the spectral/field files, the "data_dir"/"matrix_dir" JSON keys).
#ifndef SIMC_CORE_DATAPATH_H
#define SIMC_CORE_DATAPATH_H

#include <string>

namespace simc {

/// The data directory (see the header comment).
std::string DataRoot();
/// DataRoot() + "/" + sub (e.g. DataDir("pion"), DataDir("matrices/simc_fortran")).
std::string DataDir(const std::string& sub);

namespace detail {
/// Called by the installed package's interface source (step 2 above); returns true.
bool SetInstalledDataRoot(const char* dir);
}  // namespace detail

}  // namespace simc

#endif  // SIMC_CORE_DATAPATH_H
