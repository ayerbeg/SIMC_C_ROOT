// src/input/InpReader.cpp
// Fortran: dbase.f:dbase_read (l.66-122); CTP thLoad.c:getblock, thParm.c:thParmLineSet.
#include "simc/input/InpReader.h"

#include "simc/core/SpectrometerId.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>

namespace simc {

namespace {

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

bool istartsWith(const std::string& s, size_t pos, const std::string& prefix) {
    if (s.size() < pos + prefix.size()) return false;
    for (size_t i = 0; i < prefix.size(); ++i)
        if (std::tolower(static_cast<unsigned char>(s[pos + i])) !=
            std::tolower(static_cast<unsigned char>(prefix[i])))
            return false;
    return true;
}

std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

// thSpaceStrip: remove all white space.
std::string stripAllSpace(const std::string& s) {
    std::string out;
    for (char c : s)
        if (!std::isspace(static_cast<unsigned char>(c))) out += c;
    return out;
}

// CTP/daVar.h:131 floatToLong
int floatToLong(double x) { return static_cast<int>(x > 0.0 ? x + 0.5 : x - 0.5); }

// thLoad.c: types[] and qualifiers[] (prefix, case-insensitive matches)
const char* kTypes[] = {"parm", "gethit", "test", "hist", "uhist", "tree", "report", nullptr};
const char* kQualifiers[] = {"obey", "read", "write", nullptr};

// ---------------------------------------------------------------- expression evaluator
struct Expr {
    const std::string& s;
    size_t i{0};
    bool ok{true};
    explicit Expr(const std::string& str) : s(str) {}
    void ws() { while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i; }
    double primary() {
        ws();
        if (i < s.size() && s[i] == '(') {
            ++i;
            double v = sum();
            ws();
            if (i < s.size() && s[i] == ')') ++i; else ok = false;
            return v;
        }
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
            char c = s[i++];
            double v = primary();
            return c == '-' ? -v : v;
        }
        const char* start = s.c_str() + i;
        char* end = nullptr;
        double v;
        if (s.compare(i, 2, "0x") == 0 || s.compare(i, 2, "0X") == 0) {
            v = static_cast<double>(std::strtol(start, &end, 16));
        } else {
            v = std::strtod(start, &end);
            // CTP float literals (D-018): thTestParse.c:391-394 stores a TOKFLOAT literal through
            // `DAFLOAT f = atof(string)` (DAFLOAT = float, daVar.h:114). Integer literals (only digits)
            // stay exact. The arithmetic of thTestExecute.c:487-560 then promotes the operands to
            // double (dl = fl) and computes + - * / in double, and unary minus is exact, so rounding
            // each non-integer literal to float and evaluating in double reproduces CTP.
            const std::string lit(start, static_cast<size_t>(end - start));
            if (end != start && lit.find_first_of(".eE") != std::string::npos)
                v = static_cast<double>(static_cast<float>(v));
        }
        if (end == start) { ok = false; return 0.0; }
        i += static_cast<size_t>(end - start);
        return v;
    }
    double product() {
        double v = primary();
        for (;;) {
            ws();
            if (i < s.size() && (s[i] == '*' || s[i] == '/')) {
                char c = s[i++];
                double r = primary();
                v = (c == '*') ? v * r : v / r;
            } else return v;
        }
    }
    double sum() {
        double v = product();
        for (;;) {
            ws();
            if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
                char c = s[i++];
                double r = product();
                v = (c == '+') ? v + r : v - r;
            } else return v;
        }
    }
};

bool isIntToken(const std::string& t) {
    size_t i = 0;
    if (i < t.size() && (t[i] == '+' || t[i] == '-')) ++i;
    if (i >= t.size()) return false;
    if (t.compare(i, 2, "0x") == 0 || t.compare(i, 2, "0X") == 0) {
        i += 2;
        if (i >= t.size()) return false;
        for (; i < t.size(); ++i) if (!std::isxdigit(static_cast<unsigned char>(t[i]))) return false;
        return true;
    }
    for (; i < t.size(); ++i) if (!std::isdigit(static_cast<unsigned char>(t[i]))) return false;
    return true;
}

// thCommas: split on commas outside quotes.
std::vector<std::string> splitCommas(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    char q = 0;
    for (char c : s) {
        if (q) {
            if (c == q) q = 0;
            cur += c;
        } else if (c == '\'' || c == '"') {
            q = c;
            cur += c;
        } else if (c == ',') {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    out.push_back(cur);
    return out;
}

std::string dirOf(const std::string& path) {
    auto pos = path.find_last_of('/');
    return pos == std::string::npos ? std::string(".") : path.substr(0, pos);
}

}  // namespace

bool EvalCtpNumber(const std::string& text, double& value, bool& is_int) {
    const std::string t = stripAllSpace(text);
    is_int = isIntToken(t);
    if (t.empty()) return false;
    Expr e(t);
    value = e.sum();
    return e.ok && e.i == t.size();
}

void InpReader::note(const std::string& file, int line, const std::string& text) {
    messages_.push_back({file, line, text});
}

// Fortran: dbase.f:78-117
bool InpReader::loadRun(const std::string& main_file) {
    p_.extra_dbase_file.clear();   // dbase.f:76 extra_dbase_file=' '
    if (!loadFile(main_file)) return false;
    if (opt_.load_extra && !trim(p_.extra_dbase_file).empty()) {
        std::string name = trim(p_.extra_dbase_file);
        // dbase.f:106-110: 'infiles/'//name, '.inp' appended when the name has no '.'
        if (name.find('.') == std::string::npos) name += ".inp";
        const std::string dir = opt_.infiles_dir.empty() ? dirOf(main_file) : opt_.infiles_dir;
        const std::string path = dir + "/" + name;
        if (!loadFile(path)) return false;
    }
    // D-032: the retired arm code 6 (SHMS LSA tune) is rejected when the input is read (simc/core/SpectrometerId.h).
    errors_ = ValidateArmCodes(p_);
    for (const auto& e : errors_) note(main_file, 0, "ERROR: " + e);
    return errors_.empty();
}

// CTP thload (getblock) + thbook (thLoadParameters) for one file.
bool InpReader::loadFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        note(path, 0, "(thLoad) Failed to open " + path);
        return false;
    }
    files_.push_back(path);

    std::vector<std::string> raw;
    for (std::string l; std::getline(in, l);) {
        if (!l.empty() && l.back() == '\r') l.pop_back();
        raw.push_back(l);
    }

    // ---- thLoad.c:getblock — collect blocks
    std::vector<Block> blocks;
    size_t n = 0;
    while (n < raw.size()) {
        // scan for "begin"
        std::string line = raw[n];
        const int lineno = static_cast<int>(n) + 1;
        ++n;
        std::string code = line.substr(0, line.find(';'));
        std::istringstream ts(code);
        std::vector<std::string> toks;
        for (std::string t; ts >> t;) toks.push_back(t);
        if (toks.empty()) continue;
        if (lower(toks[0]) == "#include") {
            note(path, lineno, "#include is not supported by inp2json (no Fortran input uses it); line ignored");
            continue;
        }
        if (lower(toks[0]) != "begin") continue;

        size_t it = 1;
        if (it < toks.size())
            for (int q = 0; kQualifiers[q]; ++q)
                if (istartsWith(toks[it], 0, kQualifiers[q])) { ++it; break; }
        Block b;
        b.first_line = lineno;
        if (it < toks.size()) {
            int ti = 0;
            for (; kTypes[ti]; ++ti)
                if (istartsWith(toks[it], 0, kTypes[ti])) break;
            if (kTypes[ti]) b.type = kTypes[ti];
            else {
                b.type = toks[it];
                note(path, lineno, "Unknown block type " + toks[it]);
            }
            ++it;
        }
        if (it < toks.size()) b.name = toks[it];
        b.lines.push_back({lineno, line});

        bool ended = false;
        while (n < raw.size()) {
            const std::string& l = raw[n];
            const int ln = static_cast<int>(n) + 1;
            ++n;
            size_t s = 0;
            while (s < l.size() && std::isspace(static_cast<unsigned char>(l[s]))) ++s;
            if (istartsWith(l, s, "#include")) {
                note(path, ln, "#include is not supported by inp2json; line ignored");
                continue;
            }
            b.lines.push_back({ln, l});
            if (!istartsWith(l, s, "end")) continue;
            s += 4;                                        // strlen(ENDSTR)+1
            while (s < l.size() && std::isspace(static_cast<unsigned char>(l[s]))) ++s;
            if (!istartsWith(l, s, b.type)) continue;
            ended = true;
            break;
        }
        if (!ended) {
            // getblock at EOF appends a synthetic "end <type> block" line.
            note(path, b.first_line, "block '" + b.name + "' has no end line before EOF (CTP closes it)");
            b.lines.push_back({static_cast<int>(raw.size()) + 1, "end " + b.type + " block"});
        }

        // Same block name in the same file: CTP re-registers the block variable with the new text.
        auto same = std::find_if(blocks.begin(), blocks.end(), [&](const Block& x) {
            return lower(x.type) == lower(b.type) && lower(x.name) == lower(b.name);
        });
        if (same != blocks.end()) {
            note(path, b.first_line, "block '" + b.name + "' appears twice in this file; the later one replaces the earlier (CTP)");
            *same = b;
        } else {
            blocks.push_back(b);
        }
    }

    // ---- thbook: book parm blocks
    for (const auto& b : blocks) {
        if (b.type != "parm") {
            note(path, b.first_line, "block type '" + b.type + "' is not a parameter block; ignored by inp2json");
            continue;
        }
        bookParm(b, path);
    }
    return true;
}

// CTP thParm.c:thLoadParameters: first line must be the BEGIN, last line the END.
void InpReader::bookParm(const Block& b, const std::string& file) {
    cur_ = Cursor{};
    const size_t nl = b.lines.size();
    for (size_t k = 0; k < nl; ++k) {
        const auto& [ln, text] = b.lines[k];
        if (k == 0) continue;                       // the BEGIN line
        if (k + 1 == nl) {                          // the END line (or last line at EOF)
            if (lower(text).find("end") == std::string::npos)
                note(file, ln, "Last line of parm block " + b.name + " is not an END");
            break;
        }
        std::string t = trim(text);
        if (!t.empty() && t[0] == '#') {            // thSpecial directives (#real, #int, ...)
            note(file, ln, "CTP directive '" + t + "' ignored by inp2json");
            continue;
        }
        parmLine(text, file, ln);
    }
}

// CTP thParm.c:thParmLineSet
void InpReader::parmLine(std::string line, const std::string& file, int lineno) {
    // Remove comment (outside strings), turn white space into blanks.
    {
        char q = 0;
        bool blank = true;
        for (size_t i = 0; i < line.size(); ++i) {
            char& c = line[i];
            if (q) {
                if (c == q) {
                    if (i + 1 < line.size() && line[i + 1] == q) ++i;
                    else q = 0;
                }
            } else if (c == '\'' || c == '"') {
                q = c;
                blank = false;
            } else if (std::isspace(static_cast<unsigned char>(c))) {
                c = ' ';
            } else if (c == ';') {
                line.resize(i);
                break;
            } else {
                blank = false;
            }
        }
        if (blank) return;
    }

    std::string varnam, orgargs, rhs;
    const size_t eq = line.find('=');
    const bool has_name = eq != std::string::npos;
    if (has_name) {
        varnam = stripAllSpace(line.substr(0, eq));
        orgargs = line.substr(eq + 1);
        rhs = orgargs;
    } else {
        rhs = line;
    }
    std::vector<std::string> args = splitCommas(rhs);
    for (auto& a : args) a = stripAllSpace(a);
    if (!args.empty() && args.back().empty()) args.pop_back();

    if (has_name) {
        // Array element?  "(i)" is 1-based, "[i]" is 0-based (thUtils.c:thTokenArray).
        int index = 0;
        std::string base = varnam;
        size_t lp = varnam.find('(');
        size_t lb = varnam.find('[');
        if (lp != std::string::npos || lb != std::string::npos) {
            const bool cstyle = lp != std::string::npos;
            const size_t left = cstyle ? lp : lb;
            const size_t right = varnam.find(cstyle ? ')' : ']', left);
            if (right != std::string::npos && right > left) {
                index = std::atoi(varnam.substr(left + 1, right - left - 1).c_str());
                if (cstyle) --index;
                base = varnam.substr(0, left);
            }
        }
        double dummy;
        bool isint;
        if (EvalCtpNumber(base, dummy, isint)) {
            note(file, lineno, "Variable name " + base + " can't be a number");
            cur_ = Cursor{};
            return;
        }
        std::string lookup = base;
        if (istartsWith(lookup, 0, "parm.")) lookup = lookup.substr(5);
        cur_ = Cursor{};
        cur_.index = index;
        cur_.valid = true;
        cur_.desc = FindParameter(lookup);
        if (!cur_.desc) {
            const std::string key = lower(lookup);
            auto found = unreg_lower_.find(key);
            if (found == unreg_lower_.end()) {
                unreg_lower_[key] = base;
                not_registered_.push_back(base);
                note(file, lineno, base + " not registered, registering as int");
                cur_.unreg = base;
            } else {
                cur_.unreg = found->second;
            }
        }
    } else if (!cur_.valid) {
        note(file, lineno, "line without '=' and no previous variable; ignored");
        return;
    }

    // ---- Unregistered (auto-allocated dynamic) variable: keep the values for reference only.
    if (!cur_.desc) {
        auto& vals = unreg_values_[cur_.unreg];
        for (const auto& a : args) {
            double v = 0;
            bool isint;
            if (!EvalCtpNumber(a, v, isint)) note(file, lineno, "Parm: Error interpreting " + a);
            if (static_cast<int>(vals.size()) <= cur_.index) vals.resize(cur_.index + 1, 0.0);
            vals[cur_.index++] = v;
        }
        return;
    }

    const ParamDescriptor& d = *cur_.desc;
    if (d.kind == ParamKind::String) {
        if (!has_name) return;
        // thParm.c: first non-blank char after '='; if quoted, up to the matching quote.
        size_t s = 0;
        while (s < orgargs.size() && std::isspace(static_cast<unsigned char>(orgargs[s]))) ++s;
        std::string val = orgargs.substr(s);
        if (!val.empty() && (val[0] == '\'' || val[0] == '"')) {
            const char q = val[0];
            size_t e = val.find(q, 1);
            val = val.substr(1, e == std::string::npos ? std::string::npos : e - 1);
        }
        if (static_cast<int>(val.size()) > d.size) {
            note(file, lineno, std::string(d.name) + ": string truncated to character*" + std::to_string(d.size));
            val.resize(d.size);
        }
        while (!val.empty() && val.back() == ' ') val.pop_back();   // Fortran blank padding
        *d.asString(p_) = val;
        provenance_[d.name] = file + ":" + std::to_string(lineno);
        return;
    }

    const int size = (d.kind == ParamKind::FlagArray || d.kind == ParamKind::DoubleArray) ? d.size : 1;
    for (const auto& a : args) {
        if (cur_.index < 0 || cur_.index >= size) {
            note(file, lineno, std::string("Tried to fill past end of array ") + d.name);
            return;
        }
        double dval = 0.0;
        bool isint = false;
        if (!EvalCtpNumber(a, dval, isint)) {
            note(file, lineno, "Parm: Error interpreting " + a);
            ++cur_.index;
            continue;
        }
        (void)isint;
        const int lval = floatToLong(dval);   // CTP: integer variables get floatToLong(value)
        switch (d.kind) {
            case ParamKind::Int: *d.asInt(p_) = lval; break;
            case ParamKind::Flag: *d.asBool(p_) = lval > 0; break;               // convert_to_logical
            case ParamKind::Double: *d.asDouble(p_) = dval; break;
            case ParamKind::FlagArray: d.asBool(p_)[cur_.index] = lval > 0; break;
            case ParamKind::DoubleArray: d.asDouble(p_)[cur_.index] = dval; break;
            case ParamKind::String: break;
        }
        ++cur_.index;
        provenance_[d.name] = file + ":" + std::to_string(lineno);
    }
}

}  // namespace simc
