// tests/l1/L1Runner.cpp — see L1Runner.h.
#include "L1Runner.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace simcl1 {

int Table::col(const std::string& name) const {
    for (std::size_t i = 0; i < names.size(); ++i)
        if (names[i] == name) return static_cast<int>(i);
    return -1;
}

static bool blankOrComment(const std::string& line) {
    std::size_t p = line.find_first_not_of(" \t\r");
    return p == std::string::npos || line[p] == '#';
}

bool Table::load(const std::string& p, std::string& err) {
    path = p;
    names.clear();
    rows.clear();
    std::ifstream f(p);
    if (!f) {
        err = "cannot open " + p;
        return false;
    }
    std::string line;
    std::size_t lineno = 0;
    while (std::getline(f, line)) {
        ++lineno;
        if (blankOrComment(line)) continue;
        std::istringstream ss(line);
        if (names.empty()) {
            std::string w;
            while (ss >> w) names.push_back(w);
            continue;
        }
        std::vector<double> r;
        std::string w;
        while (ss >> w) {
            char* end = nullptr;
            // strtod accepts 3-digit exponents, NaN and Infinity (any case).
            double v = std::strtod(w.c_str(), &end);
            if (end == w.c_str() || *end != '\0') {
                err = p + ":" + std::to_string(lineno) + ": cannot parse '" + w + "'";
                return false;
            }
            r.push_back(v);
        }
        if (r.size() != names.size()) {
            err = p + ":" + std::to_string(lineno) + ": " + std::to_string(r.size()) + " values, header has " +
                  std::to_string(names.size());
            return false;
        }
        rows.push_back(std::move(r));
    }
    if (names.empty()) {
        err = p + ": no header";
        return false;
    }
    return true;
}

double Row::operator[](const std::string& name) const {
    int c = t_.col(name);
    if (c < 0) throw std::out_of_range("L1: no column '" + name + "' in " + t_.path);
    return t_.rows[i_][static_cast<std::size_t>(c)];
}

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

double relDiff(double cpp, double ref) {
    if (cpp == ref) return 0.0;
    double d = std::fabs(cpp - ref);
    if (ref == 0.0) return std::numeric_limits<double>::infinity();
    return d / std::fabs(ref);
}

std::vector<ColumnResult> compare(const Table& in, const Table& ref, const std::vector<ColumnSpec>& cols,
                                  const Adapter& adapter, std::string& err) {
    std::vector<ColumnResult> res;
    if (in.rows.size() != ref.rows.size()) {
        err = "row count mismatch: " + in.path + " has " + std::to_string(in.rows.size()) + ", " + ref.path + " has " +
              std::to_string(ref.rows.size());
        return res;
    }
    // Echoed input columns must be bit-identical.
    for (std::size_t c = 0; c < in.names.size(); ++c) {
        int rc = ref.col(in.names[c]);
        if (rc < 0) continue;
        for (std::size_t i = 0; i < in.rows.size(); ++i) {
            double a = in.rows[i][c], b = ref.rows[i][static_cast<std::size_t>(rc)];
            if (!sameBits(a, b) && !(std::isnan(a) && std::isnan(b))) {
                err = "echoed input column '" + in.names[c] + "' differs at row " + std::to_string(i + 1);
                return res;
            }
        }
    }
    std::vector<int> refcol;
    for (const auto& s : cols) {
        int rc = ref.col(s.name);
        if (rc < 0) {
            err = "reference has no output column '" + s.name + "'";
            return res;
        }
        refcol.push_back(rc);
        ColumnResult r;
        r.name = s.name;
        res.push_back(r);
    }
    Outputs out;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        out.clear();
        adapter(Row(in, i), out);
        for (std::size_t k = 0; k < cols.size(); ++k) {
            auto it = out.find(cols[k].name);
            if (it == out.end()) {
                err = "adapter did not set column '" + cols[k].name + "'";
                return {};
            }
            double c = it->second, f = ref.rows[i][static_cast<std::size_t>(refcol[k])];
            ColumnResult& r = res[k];
            ++r.rows;
            if (sameBits(c, f) || (std::isnan(c) && std::isnan(f))) {
                ++r.bitwise_equal;
                continue;
            }
            if (std::isnan(c) != std::isnan(f)) {
                ++r.nan_mismatch;
                ++r.failed;
                if (!std::isinf(r.max_rel)) {   // first NaN mismatch becomes the worst row
                    r.max_rel = std::numeric_limits<double>::infinity();
                    r.worst_row = i + 1;
                    r.worst_ref = f;
                    r.worst_cpp = c;
                }
                continue;
            }
            double a = std::fabs(c - f), rel = relDiff(c, f);
            if (a > r.max_abs) r.max_abs = a;
            if (r.worst_row == 0 || rel > r.max_rel) {
                r.max_rel = rel;
                r.worst_row = i + 1;
                r.worst_ref = f;
                r.worst_cpp = c;
            }
            if (!(rel <= cols[k].rel_tol || a <= cols[k].abs_tol)) ++r.failed;
        }
    }
    return res;
}

bool report(const std::string& title, const std::vector<ColumnResult>& res, const std::vector<ColumnSpec>& cols) {
    bool ok = !res.empty();
    std::printf("L1 %s\n", title.c_str());
    std::printf("  %-12s %6s %8s %6s %5s %12s %12s %8s  %-24s %-24s  tol(rel/abs)\n", "column", "rows", "bitwise",
                "fail", "NaN!", "max|rel|", "max|abs|", "worstrow", "ref(worst)", "cpp(worst)");
    for (std::size_t k = 0; k < res.size(); ++k) {
        const auto& r = res[k];
        std::printf("  %-12s %6zu %8zu %6zu %5zu %12.4e %12.4e %8zu  %24.16e %24.16e  %.1e/%.1e\n", r.name.c_str(), r.rows,
                    r.bitwise_equal, r.failed, r.nan_mismatch, r.max_rel, r.max_abs, r.worst_row, r.worst_ref,
                    r.worst_cpp, cols[k].rel_tol, cols[k].abs_tol);
        if (r.failed) ok = false;
    }
    std::printf("  => %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

int run(const std::string& title, const std::string& dir, const std::vector<ColumnSpec>& cols, const Adapter& adapter) {
    Table in, ref;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !ref.load(dir + "/output.txt", err)) {
        std::printf("L1 %s: SKIP (%s)\n", title.c_str(), err.c_str());
        return 77;
    }
    auto res = compare(in, ref, cols, adapter, err);
    if (!err.empty()) {
        std::printf("L1 %s: ERROR %s\n", title.c_str(), err.c_str());
        return 1;
    }
    return report(title, res, cols) ? 0 : 1;
}

}  // namespace simcl1
