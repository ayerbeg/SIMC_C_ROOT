// src/core/DataPath.cpp — see include/simc/core/DataPath.h (task P15.4, D-069).
#include "simc/core/DataPath.h"

#include <cstdlib>
#include <filesystem>
#include <system_error>

#ifndef SIMC_SOURCE_DATA_ROOT
#define SIMC_SOURCE_DATA_ROOT "data"
#endif

namespace simc {

namespace {
// Set during static initialisation (before main) by the installed package's interface source; read-only afterwards.
std::string& installedRoot() {
    static std::string r;
    return r;
}
bool isDir(const std::filesystem::path& p) {
    std::error_code ec;
    return std::filesystem::is_directory(p, ec);
}
}  // namespace

namespace detail {
bool SetInstalledDataRoot(const char* dir) {
    installedRoot() = dir ? dir : "";
    return true;
}
}  // namespace detail

std::string DataRoot() {
    if (const char* e = std::getenv("SIMC_DATA_DIR"); e && *e) return e;                    // 1
    if (!installedRoot().empty() && isDir(installedRoot())) return installedRoot();           // 2
    std::error_code ec;                                                                       // 3
    const std::filesystem::path exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec) {
        const std::filesystem::path d = exe.parent_path().parent_path() / "share" / "simc" / "data";
        if (isDir(d)) return d.lexically_normal().string();
    }
    return SIMC_SOURCE_DATA_ROOT;                                                             // 4
}

std::string DataDir(const std::string& sub) { return DataRoot() + "/" + sub; }

}  // namespace simc
