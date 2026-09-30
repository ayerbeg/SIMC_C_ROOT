// tests/unit/test_hist2json.cpp
//   test_hist2json <run.hist> <docs/summary.schema.json>
// Parses a Fortran .hist (+ .geni next to it) and checks that
//  * every line is understood (no parser warnings),
//  * every standalone number printed in the .hist appears in the JSON with the same value,
//  * every integer printed in the .geni appears in stops (the SHMS Q3 "mid" column is stored
//    as q2_mid_dup, since the Fortran prints shmsSTOP_q2_mid there),
//  * the result validates against the schema.
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>

#include "TestUtil.h"
#include "simc/output/SummaryJson.h"

using namespace simc;
using simctest::check;

namespace {

std::string slurp(const std::string& p) {
    std::ifstream f(p);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

void collect(const Json& n, std::multiset<double>& out) {
    if (n.is_number()) out.insert(n.get<double>());
    else if (n.is_boolean()) out.insert(n.get<bool>() ? 1.0 : 0.0);
    else if (n.is_structured()) for (const auto& c : n) collect(c, out);
}

bool isNumberToken(const std::string& t, double& v) {
    if (t.find_first_of("0123456789") == std::string::npos) return false;
    char* end = nullptr;
    v = std::strtod(t.c_str(), &end);
    return end && *end == '\0';
}

bool present(const std::multiset<double>& s, double v) {
    auto it = s.lower_bound(v - std::abs(v) * 1e-12 - 1e-300);
    return it != s.end() && std::abs(*it - v) <= std::abs(v) * 1e-12 + 1e-300;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string hist = argc > 1 ? argv[1] : "";
    const std::string schema_path = argc > 2 ? argv[2] : "";
    if (!simctest::exists(hist)) { std::printf("SKIP: %s missing\n", hist.c_str()); return simctest::kSkip; }
    std::string geni = hist.substr(0, hist.size() - 5) + ".geni";
    const std::string htxt = slurp(hist), gtxt = simctest::exists(geni) ? slurp(geni) : "";

    std::vector<std::string> warnings;
    Json j = ParseFortranHist(htxt, gtxt, "test", &warnings);
    for (const auto& w : warnings) std::printf("  warning: %s\n", w.c_str());
    check(warnings.empty(), "no unparsed lines");

    std::multiset<double> values;
    collect(j, values);

    // every standalone number of the .hist
    int n_numbers = 0, n_missing = 0;
    {
        std::istringstream is(htxt);
        for (std::string line; std::getline(is, line);) {
            if (line.find("Time:") != std::string::npos) continue;          // dates
            std::istringstream ls(line);
            for (std::string t; ls >> t;) {
                double v;
                if (!t.empty() && t.back() == ',') t.pop_back();
                if (!isNumberToken(t, v)) continue;
                ++n_numbers;
                if (!present(values, v)) {
                    ++n_missing;
                    std::printf("  missing %s from line: %s\n", t.c_str(), line.c_str());
                }
            }
        }
    }
    std::printf("hist numbers checked: %d, geni: see below\n", n_numbers);
    check(n_missing == 0, "all " + std::to_string(n_numbers) + " printed .hist numbers are in the JSON");

    // every integer of the .geni
    int g_numbers = 0, g_missing = 0;
    std::multiset<double> stops;
    collect(j["stops"], stops);
    {
        std::istringstream is(gtxt);
        for (std::string line; std::getline(is, line);) {
            std::istringstream ls(line);
            int k = 0;
            const bool q3 = line.find("Q3 phys") != std::string::npos;
            for (std::string t; ls >> t;) {
                double v;
                if (t.find_first_not_of("-0123456789") != std::string::npos || !isNumberToken(t, v)) continue;
                ++k;
                (void)q3;   // the Q3 "mid" column is stored as q2_mid_dup
                ++g_numbers;
                if (!present(stops, v)) { ++g_missing; std::printf("  geni missing %s: %s\n", t.c_str(), line.c_str()); }
            }
        }
    }
    std::printf("geni counters checked: %d\n", g_numbers);
    check(g_missing == 0, "all " + std::to_string(g_numbers) + " .geni counters are in stops");

    // Placement of the key normalisation numbers (not just presence).
    auto printed = [&](const std::string& re) -> double {
        std::smatch m;
        std::regex r(re);
        if (std::regex_search(htxt, m, r)) return std::strtod(m[1].str().c_str(), nullptr);
        return std::nan("");
    };
    auto same = [&](const std::string& what, double a, const Json& b) {
        check(b.is_number() && std::abs(a - b.get<double>()) <= 1e-12 * std::abs(a), what + " placed correctly");
    };
    same("ngen", printed(R"(Ngen \(request\) =\s*(\S+))"), j["counters"]["ngen"]);
    same("ntried", printed(R"(Ntried\s*=\s*(\S+))"), j["counters"]["ntried"]);
    same("ncontribute", printed(R"(Ncontribute\s*=\s*(\S+))"), j["counters"]["ncontribute"]);
    same("wtcontr", printed(R"(wtcontr=\s*(\S+))"), j["normalization"]["wtcontr"]);
    same("central_sigcc", printed(R"(CENTRAL\.sigcc =\s*(\S+))"), j["normalization"]["sigcc_central"]);
    same("average_sigcc", printed(R"(AVERAGE\.sigcc =\s*(\S+))"), j["normalization"]["sigcc_average"]);
    same("charge", printed(R"(charge =\s*(\S+))"), j["normalization"]["charge_mC"]);
    same("targetfac", printed(R"(targetfac =\s*(\S+))"), j["normalization"]["targetfac"]);
    same("luminosity_ub", printed(R"(luminosity =\s*(\S+)\s+ub)"), j["normalization"]["luminosity_ub"]);
    same("luminosity_GeV2", printed(R"(luminosity =\s*(\S+)\s+GeV)"), j["normalization"]["luminosity_GeV2"]);
    same("genvol", printed(R"(genvol =\s*(\S+))"), j["normalization"]["genvol"]);
    same("normfac", printed(R"(normfac =\s*(\S+))"), j["normalization"]["normfac"]);

    Json schema = Json::parse(slurp(schema_path));
    auto errs = ValidateAgainstSchema(j, schema);
    for (const auto& e : errs) std::printf("  schema: %s\n", e.c_str());
    check(errs.empty(), "hist2json output validates against the schema");

    std::printf("key: ngen=%s ntried=%s ncontribute=%s wtcontr=%s central_sigcc=%s average_sigcc=%s "
                "luminosity_ub=%s genvol=%s normfac=%s\n",
                j["counters"]["ngen"].dump().c_str(), j["counters"]["ntried"].dump().c_str(),
                j["counters"]["ncontribute"].dump().c_str(), j["normalization"]["wtcontr"].dump().c_str(),
                j["normalization"]["sigcc_central"].dump().c_str(), j["normalization"]["sigcc_average"].dump().c_str(),
                j["normalization"]["luminosity_ub"].dump().c_str(), j["normalization"]["genvol"].dump().c_str(),
                j["normalization"]["normfac"].dump().c_str());
    for (auto it = j["stops"].begin(); it != j["stops"].end(); ++it)
        std::printf("stops %s (%s): trials=%s successes=%s\n", it.key().c_str(),
                    j["stop_roles"].value(it.key(), std::string("?")).c_str(),
                    it.value().value("trials", Json()).dump().c_str(), it.value().value("successes", Json()).dump().c_str());
    return simctest::finish("hist2json");
}
