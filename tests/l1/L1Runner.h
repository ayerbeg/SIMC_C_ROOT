// tests/l1/L1Runner.h
//
// Generic L1 (function-level) test runner for the Fortran harness tables in
// reference/functions/<routine>/ (format: reference/functions/README.md):
//   * first non-comment line = header with column names; '#' lines are comments;
//   * values in %24.16E / ES24.16E3 (17 significant digits: text -> double is exact);
//   * output.txt echoes the input columns and appends the output columns, same row order.
//
// A routine needs only a short adapter: given one input row (columns by name) it evaluates
// the C++ and fills the output columns. The runner compares every requested output column
// with a per-column tolerance and prints max|rel|, max|abs|, the worst row, the number of
// bitwise-equal rows and NaN mismatches. Exit codes: 0 pass, 1 fail, 77 table missing (skip).
#ifndef SIMC_TESTS_L1_L1RUNNER_H
#define SIMC_TESTS_L1_L1RUNNER_H

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "Upstream294.h"

namespace simcl1 {

/// A whitespace-separated harness table.
struct Table {
    std::vector<std::string> names;
    std::vector<std::vector<double>> rows;
    std::string path;

    /// Index of column `name`, or -1.
    int col(const std::string& name) const;
    /// Load `path`. Returns false (and sets `err`) if missing or malformed.
    bool load(const std::string& path, std::string& err);
};

/// Read access to one row by column name (throws std::out_of_range on unknown names).
class Row {
public:
    Row(const Table& t, std::size_t i) : t_(t), i_(i) {}
    double operator[](const std::string& name) const;
    bool has(const std::string& name) const { return t_.col(name) >= 0; }
    std::size_t index() const { return i_; }        ///< 0-based
    std::size_t number() const { return i_ + 1; }   ///< 1-based (the Fortran drivers' row counter)
private:
    const Table& t_;
    std::size_t i_;
};

/// Output column and its tolerance: pass if |rel| <= rel_tol or |abs| <= abs_tol.
struct ColumnSpec {
    std::string name;
    double rel_tol{1e-12};
    double abs_tol{0.0};
};

struct ColumnResult {
    std::string name;
    std::size_t rows{0}, bitwise_equal{0}, failed{0}, nan_mismatch{0};
    double max_rel{0}, max_abs{0};
    std::size_t worst_row{0};           ///< 1-based row of max_rel (0 if all equal)
    double worst_ref{0}, worst_cpp{0};
};

/// Output values produced by an adapter for one row.
using Outputs = std::map<std::string, double>;
using Adapter = std::function<void(const Row& in, Outputs& out)>;

/// Relative difference used by the runner: |c-f|/|f|; 0 if c==f; inf if f==0 != c.
double relDiff(double cpp, double ref);

/// Compare the adapter's outputs with the reference output table.
/// `in` and `ref` must have the same number of rows; echoed input columns of `ref` are
/// checked to be bit-identical to `in` (guards against a stale/misaligned table).
std::vector<ColumnResult> compare(const Table& in, const Table& ref, const std::vector<ColumnSpec>& cols,
                                  const Adapter& adapter, std::string& err);

/// Print the per-column report; returns true if all columns pass.
bool report(const std::string& title, const std::vector<ColumnResult>& res, const std::vector<ColumnSpec>& cols);

/// Full test: load <dir>/input.txt and <dir>/output.txt, run, report.
/// Returns 0 (pass), 1 (fail / error), 77 (tables missing -> CTest SKIP).
int run(const std::string& title, const std::string& dir, const std::vector<ColumnSpec>& cols, const Adapter& adapter);

}  // namespace simcl1

#endif  // SIMC_TESTS_L1_L1RUNNER_H
