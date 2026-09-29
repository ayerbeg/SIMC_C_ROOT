// tests/unit/TestUtil.h — minimal check helpers for the framework unit tests (no external deps).
#ifndef SIMC_TESTS_TESTUTIL_H
#define SIMC_TESTS_TESTUTIL_H

#include <cstdio>
#include <fstream>
#include <string>

namespace simctest {

inline int& failures() { static int n = 0; return n; }
inline int& checks() { static int n = 0; return n; }

inline void check(bool ok, const std::string& what) {
    ++checks();
    if (!ok) {
        ++failures();
        std::printf("FAIL: %s\n", what.c_str());
    }
}

inline bool exists(const std::string& p) { return std::ifstream(p).good(); }

/// CTest skip code (SKIP_RETURN_CODE 77).
constexpr int kSkip = 77;

inline int finish(const char* name) {
    std::printf("%s: %d checks, %d failures\n", name, checks(), failures());
    return failures() ? 1 : 0;
}

}  // namespace simctest

#endif
