// tests/unit/test_parameters.cpp
//   completeness <inp_variables.csv> : every registered .inp variable has a SimcParameters
//       field with the right CTP kind and the true Fortran default; no extra fields; the
//       46 obsolete '.' names are NOT registered.
//   json : ToJson/FromJson round trip, unknown-key warnings.
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TestUtil.h"
#include "simc/core/SimcParameters.h"

using namespace simc;
using simctest::check;

namespace {

// Minimal RFC-4180 CSV reader (quoted fields may contain commas and doubled quotes).
std::vector<std::vector<std::string>> readCsv(const std::string& path) {
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string s = ss.str();
    std::vector<std::vector<std::string>> rows;
    std::vector<std::string> row;
    std::string field;
    bool q = false;
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (q) {
            if (c == '"') {
                if (i + 1 < s.size() && s[i + 1] == '"') { field += '"'; ++i; }
                else q = false;
            } else field += c;
        } else if (c == '"') q = true;
        else if (c == ',') { row.push_back(field); field.clear(); }
        else if (c == '\n') { row.push_back(field); field.clear(); rows.push_back(row); row.clear(); }
        else if (c != '\r') field += c;
    }
    if (!field.empty() || !row.empty()) { row.push_back(field); rows.push_back(row); }
    return rows;
}

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// First number in the CSV default text ("0 (.false.)", "-1.0d10", "+1.0d10", "0.0 x3", "all 0 ...").
bool firstNumber(std::string t, double& v) {
    for (auto& c : t) if (c == 'd' || c == 'D') c = 'e';
    for (size_t i = 0; i < t.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(t[i])) ||
            ((t[i] == '-' || t[i] == '+') && i + 1 < t.size() && std::isdigit(static_cast<unsigned char>(t[i + 1])))) {
            v = std::strtod(t.c_str() + i, nullptr);
            return true;
        }
    }
    return false;
}

ParamKind expectedKind(const std::string& reg, const std::string& ftype) {
    if (reg.find("regparmintarray") != std::string::npos) return ParamKind::FlagArray;
    if (reg.find("regparmdoublearray") != std::string::npos) return ParamKind::DoubleArray;
    if (reg.find("regparmstring") != std::string::npos) return ParamKind::String;
    if (reg.find("regparmdouble") != std::string::npos) return ParamKind::Double;
    // regparmint: stored as bool when convert_to_logical turns it into a logical
    if (ftype.find("logical") != std::string::npos) return ParamKind::Flag;
    return ParamKind::Int;
}

int completeness(const std::string& csv) {
    if (!simctest::exists(csv)) {
        std::printf("SKIP: %s not found\n", csv.c_str());
        return simctest::kSkip;
    }
    auto rows = readCsv(csv);
    check(!rows.empty(), "csv has rows");
    std::map<std::string, int> col;
    for (size_t i = 0; i < rows[0].size(); ++i) col[rows[0][i]] = static_cast<int>(i);
    for (const char* need : {"name", "kind", "ctp_registration", "fortran_type", "default_if_absent"})
        check(col.count(need) == 1, std::string("csv column ") + need);

    const SimcParameters def{};
    std::set<std::string> csv_registered;
    int n_registered = 0, n_legacy = 0;
    for (size_t r = 1; r < rows.size(); ++r) {
        const auto& row = rows[r];
        if (row.size() < rows[0].size()) continue;
        const std::string name = row[col["name"]];
        const std::string kind = row[col["kind"]];
        if (kind == "registered") {
            ++n_registered;
            csv_registered.insert(lower(name));
            const ParamDescriptor* d = FindParameter(name);
            check(d != nullptr, "SimcParameters field for registered variable " + name);
            if (!d) continue;
            check(std::string(d->name) == name, "exact CTP spelling of " + name + " (have " + d->name + ")");
            const ParamKind ek = expectedKind(row[col["ctp_registration"]], row[col["fortran_type"]]);
            check(d->kind == ek, "storage kind of " + name);
            // true default
            const std::string dtext = row[col["default_if_absent"]];
            auto& p = const_cast<SimcParameters&>(def);
            if (d->kind == ParamKind::String) {
                check(d->asString(p)->empty(), "default of " + name + " is blank (" + dtext + ")");
            } else {
                double v = 0;
                check(firstNumber(dtext, v), "numeric default parsed for " + name + ": '" + dtext + "'");
                const int n = (d->kind == ParamKind::FlagArray || d->kind == ParamKind::DoubleArray) ? d->size : 1;
                for (int i = 0; i < n; ++i) {
                    double have = 0;
                    switch (d->kind) {
                        case ParamKind::Int: have = *d->asInt(p); break;
                        case ParamKind::Flag:
                        case ParamKind::FlagArray: have = d->asBool(p)[i] ? 1 : 0; break;
                        case ParamKind::Double:
                        case ParamKind::DoubleArray: have = d->asDouble(p)[i]; break;
                        default: break;
                    }
                    const bool ok = (d->kind == ParamKind::Flag || d->kind == ParamKind::FlagArray)
                                        ? ((have > 0) == (v > 0))
                                        : std::abs(have - v) <= 1e-12 * std::max(1.0, std::abs(v));
                    check(ok, "default of " + name + " = " + std::to_string(have) + " expected '" + dtext + "'");
                }
            }
        } else {
            ++n_legacy;
            check(FindParameter(name) == nullptr, "obsolete name " + name + " must NOT be registered");
        }
    }
    check(n_registered == 106, "csv lists 106 registered variables (got " + std::to_string(n_registered) + ")");
    check(n_legacy == 46, "csv lists 46 obsolete names (got " + std::to_string(n_legacy) + ")");
    // U.1: variables registered only by the upstream SIMC (59ccdf6), absent from the 294f662 inventory csv:
    // doing_nuc_elast (dbase.f:1007, commit 4df7351; integer flag, default 0).
    const std::set<std::string> upstream_only = {"doing_nuc_elast"};
    for (const auto& u : upstream_only) {
        const ParamDescriptor* d = FindParameter(u);
        check(d != nullptr && d->kind == ParamKind::Flag && !*d->asBool(const_cast<SimcParameters&>(def)),
              "upstream-only flag " + u + " registered, default F");
    }
    check(ParameterTable().size() == 106 + upstream_only.size(), "ParameterTable has 106 + " +
              std::to_string(upstream_only.size()) + " upstream entries (got " + std::to_string(ParameterTable().size()) + ")");
    for (const auto& d : ParameterTable())
        check(csv_registered.count(lower(d.name)) == 1 || upstream_only.count(d.name) == 1,
              std::string("table entry ") + d.name + " is in the csv (or an upstream-only variable)");
    std::printf("registered=%d obsolete=%d table=%zu\n", n_registered, n_legacy, ParameterTable().size());
    return simctest::finish("params_completeness");
}

int json() {
    SimcParameters a;
    a.Ebeam = 10600.0;
    a.spec.e.theta = -12.5;
    a.targ.fr_pattern = 3;
    a.debug[1] = true;
    a.lambda = {0.1, 0.2, 0.3};
    a.SPedge.p.xptar.max = 55.5;
    a.extra_dbase_file = "extra_hydrogen";
    a.using_rad = true;
    a.cuts.Em.min = 0.0;

    Json j = ToJson(a);
    check(j["targ"]["fr_pattern"] == 3, "json path targ%fr_pattern -> targ.fr_pattern");
    check(j["SPedge"]["p"]["xptar"]["max"] == 55.5, "json path SPedge%p%xptar%max");
    check(j["EXPER"]["charge"] == 0.0, "json path EXPER%charge");
    check(j["debug"].size() == 6 && j["lambda"].size() == 3, "arrays have Fortran sizes");

    SimcParameters b;
    std::vector<std::string> w;
    FromJson(j, b, &w);
    check(w.empty(), "round trip has no warnings");
    check(ToJson(b) == j, "round trip preserves every field");

    // Integers given as numbers for flags, case-insensitive keys, unknown keys.
    Json k = Json::parse(R"({"ebeam": 5000, "using_rad": 1, "mc_smear": -1, "targ": {"FR_PATTERN": 1.6, "bogus": 3},
                             "not_a_param": 1, "_comment": "ignored", "debug": [0, 2]})");
    SimcParameters c;
    w.clear();
    FromJson(k, c, &w);
    check(c.Ebeam == 5000.0, "case-insensitive key ebeam");
    check(c.using_rad, "flag from integer 1");
    check(!c.mc_smear, "flag from integer -1 is .false. (convert_to_logical: >0)");
    check(c.targ.fr_pattern == 2, "real into int is rounded like CTP floatToLong");
    check(c.debug[1] && !c.debug[0], "debug array");
    int unknown = 0;
    for (const auto& s : w) if (s.find("unknown parameter") != std::string::npos) ++unknown;
    check(unknown == 2, "two unknown-key warnings (targ%bogus, not_a_param), got " + std::to_string(unknown));
    return simctest::finish("params_json_roundtrip");
}

}  // namespace

int main(int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "";
    if (mode == "completeness") return completeness(argc > 2 ? argv[2] : "");
    if (mode == "json") return json();
    std::printf("usage: test_parameters completeness <csv> | json\n");
    return 2;
}
