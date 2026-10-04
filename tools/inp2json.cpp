// tools/inp2json.cpp
//
// Convert a Fortran SIMC input file (.inp, CTP syntax, plus its extra_dbase_file)
// into the SIMC++ JSON configuration (D-003).
//
// Usage:
//   inp2json [options] file.inp
//     -o FILE            write the JSON to FILE (default: stdout)
//     --infiles-dir DIR  directory for extra_dbase_file (Fortran: 'infiles/'); default: dir of file.inp
//     --no-extra         do not load extra_dbase_file
//     --table[=FILE]     also print the Fortran-equivalent effective-value table (stdout or FILE)
//                        (with --table and no -o, only the table is written to stdout)
//     --mapping          print the .inp -> JSON mapping table (markdown) and exit
//     --strict           exit with status 3 if any name is "not registered"
//     -q, --quiet        do not echo CTP messages to stderr
//
// CTP messages ("<name> not registered, registering as int", ...) are echoed to stderr
// in the same words as the Fortran prints them. Exit status: 0 ok, 1 usage, 2 file error.
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "simc/core/SimcParameters.h"
#include "simc/input/EffectiveValues.h"
#include "simc/input/InpReader.h"

using namespace simc;

namespace {

const char* kindName(ParamKind k) {
    switch (k) {
        case ParamKind::Int: return "int";
        case ParamKind::Flag: return "bool (int>0)";
        case ParamKind::Double: return "double";
        case ParamKind::String: return "string";
        case ParamKind::FlagArray: return "bool[]";
        case ParamKind::DoubleArray: return "double[]";
    }
    return "?";
}

void printMapping() {
    const SimcParameters def;
    std::cout << "| # | Fortran (.inp) name | JSON path under \"parameters\" | type | .inp units | default | Fortran variable (registration) |\n"
                 "|---|---|---|---|---|---|---|\n";
    int i = 0;
    for (const auto& d : ParameterTable()) {
        std::string path;
        for (const auto& c : d.jsonPath()) path += (path.empty() ? "" : ".") + c;
        std::string kind = kindName(d.kind);
        if (d.kind == ParamKind::FlagArray || d.kind == ParamKind::DoubleArray || d.kind == ParamKind::String)
            kind += "[" + std::to_string(d.size) + "]";
        std::cout << "| " << ++i << " | `" << d.name << "` | `" << path << "` | " << kind << " | " << d.units
                  << " | " << FormatValue(d, def) << " | " << d.fortran << " |\n";
    }
}

void usage() {
    std::cerr << "usage: inp2json [-o out.json] [--infiles-dir DIR] [--no-extra] [--table[=FILE]] "
                 "[--strict] [-q] file.inp\n";
}

bool fileExists(const std::string& p) { return std::ifstream(p).good(); }

}  // namespace

int main(int argc, char** argv) {
    std::string input, output, table_file, infiles_dir;
    bool want_table = false, strict = false, quiet = false, load_extra = true;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-o" && i + 1 < argc) output = argv[++i];
        else if (a == "--infiles-dir" && i + 1 < argc) infiles_dir = argv[++i];
        else if (a == "--no-extra") load_extra = false;
        else if (a == "--table") want_table = true;
        else if (a.rfind("--table=", 0) == 0) { want_table = true; table_file = a.substr(8); }
        else if (a == "--strict") strict = true;
        else if (a == "--mapping") { printMapping(); return 0; }
        else if (a == "-q" || a == "--quiet") quiet = true;
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (!a.empty() && a[0] == '-') { usage(); return 1; }
        else if (input.empty()) input = a;
        else { usage(); return 1; }
    }
    if (input.empty()) { usage(); return 1; }

    // dbase.f:84-88: '.inp' is appended when the file name has no '.'.
    if (!fileExists(input) && input.find('.', input.find_last_of('/') == std::string::npos
                                                   ? 0 : input.find_last_of('/')) == std::string::npos)
        input += ".inp";

    SimcParameters params;
    InpReader::Options opt;
    opt.infiles_dir = infiles_dir;
    opt.load_extra = load_extra;
    InpReader reader(params, opt);
    const bool ok = reader.loadRun(input);

    if (!quiet)
        for (const auto& m : reader.messages()) std::cerr << m.text << "\n";
    if (!ok) {
        for (const auto& e : reader.errors()) std::cerr << "inp2json: input error: " << e << "\n";
        std::cerr << "inp2json: could not read " << input << "\n";
        return 2;
    }

    Json out;
    out["_source"]["generator"] = "inp2json (SIMC++ P1.3)";
    out["_source"]["files"] = reader.filesLoaded();
    out["_source"]["note"] =
        "parameters = every dbase.f:regallvars variable, in .inp units; extra_dbase_file is already merged";
    out["_source"]["not_registered"] = Json::object();
    for (const auto& [name, vals] : reader.unregisteredValues()) out["_source"]["not_registered"][name] = vals;
    out["_source"]["set_by"] = reader.provenance();
    out["parameters"] = ToJson(params);

    if (output.empty()) {
        if (!(want_table && table_file.empty())) std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream f(output);
        if (!f) { std::cerr << "inp2json: cannot write " << output << "\n"; return 2; }
        f << out.dump(2) << "\n";
    }

    if (want_table) {
        const std::string t = FormatEffectiveTable(EffectiveValues(params));
        if (table_file.empty()) std::cout << t;
        else {
            std::ofstream f(table_file);
            if (!f) { std::cerr << "inp2json: cannot write " << table_file << "\n"; return 2; }
            f << t;
        }
    }

    if (strict && !reader.notRegistered().empty()) return 3;
    return 0;
}
