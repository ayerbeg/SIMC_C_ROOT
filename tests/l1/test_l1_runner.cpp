// tests/l1/test_l1_runner.cpp — self-test of the generic L1 runner (table parsing,
// bitwise/rel/abs accounting, NaN mismatches, worst row, echoed-input check, skip code).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>

#include "L1Runner.h"

using namespace simcl1;

static int fails = 0;
static void check(bool ok, const char* what) {
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        ++fails;
    }
}

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : ".";
    {
        std::ofstream f(dir + "/l1selftest_input.txt");
        f << "x y\n# comment line\n  1.0000000000000000E+00  2.0E+00\n  3.0 4.0\n -5.0 0.0\n 7.0 1.0\n";
        std::ofstream g(dir + "/l1selftest_output.txt");
        g << "  x  y  s  n\n# c\n 1.0000000000000000E+000 2.0 3.0000000000000000E+000 1.0\n"
             " 3.0 4.0 7.0000000100000000E+000 NaN\n -5.0 0.0 -5.0 0.0\n 7.0 1.0 0.0 1.0\n";
    }
    Table in, out;
    std::string err;
    check(in.load(dir + "/l1selftest_input.txt", err), "load input");
    check(out.load(dir + "/l1selftest_output.txt", err), "load output");
    check(in.rows.size() == 4 && in.names.size() == 2, "input shape");
    check(out.col("s") == 2 && out.col("zz") == -1, "col lookup");
    check(std::isnan(out.rows[1][3]), "NaN parsed");

    std::vector<ColumnSpec> cols = {{"s", 1e-12, 0.0}, {"n", 1e-12, 0.0}};
    auto res = compare(in, out, cols, [](const Row& r, Outputs& o) {
        o["s"] = r["x"] + r["y"];        // row 2: 7 vs 7.00000001 (rel 1.4e-9), row 4: 8 vs 0
        o["n"] = r["y"] > 0 ? 1.0 : 0.0; // row 2: 1 vs NaN
    }, err);
    check(err.empty(), "compare ok");
    check(res.size() == 2, "two results");
    if (res.size() == 2) {
        check(res[0].rows == 4 && res[0].bitwise_equal == 2, "s bitwise count");
        check(res[0].failed == 2, "s failures (row 2 above 1e-12, row 4 ref 0)");
        check(res[0].worst_row == 4 && std::isinf(res[0].max_rel), "s worst row 4 (ref 0)");
        check(res[0].max_abs == 8.0, "s max abs");
        check(res[1].nan_mismatch == 1 && res[1].worst_row == 2, "n NaN mismatch at row 2");
    }
    check(!report("selftest (expected FAIL)", res, cols), "report returns false");

    // Tolerances: rel 1e-8 passes row 2, abs 10 passes row 4 (ref 0).
    std::vector<ColumnSpec> loose = {{"s", 1e-8, 10.0}};
    res = compare(in, out, loose, [](const Row& r, Outputs& o) { o["s"] = r["x"] + r["y"]; }, err);
    check(res.size() == 1 && res[0].failed == 0, "abs_tol / rel_tol respected");

    // Misaligned echo column is detected.
    Table bad = out;
    bad.rows[0][0] = 1.5;
    err.clear();
    compare(in, bad, cols, [](const Row&, Outputs& o) { o["s"] = 0; o["n"] = 0; }, err);
    check(!err.empty(), "echoed input mismatch detected");

    check(run("missing", dir + "/no_such_dir", cols, [](const Row&, Outputs&) {}) == 77, "missing table -> 77");
    check(relDiff(1.0, 1.0) == 0.0 && relDiff(0.0, 0.0) == 0.0 && std::isinf(relDiff(1.0, 0.0)), "relDiff");
    std::printf("test_l1_runner: %s\n", fails ? "FAIL" : "PASS");
    return fails ? 1 : 0;
}
