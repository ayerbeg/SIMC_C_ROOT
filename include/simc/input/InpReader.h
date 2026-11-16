/**
 * @file
 * @ingroup input
 * @brief simc::InpReader: reader of the Fortran SIMC input files (.inp, CTP "parm" syntax) into SimcParameters.
 *
 * Reader for Fortran SIMC input files (.inp, CTP "parm" syntax) into SimcParameters.
 *
 * Fortran: dbase.f:dbase_read (l.66-122) — thload/thbook of the main file, then of
 *          'infiles/'//extra_dbase_file (+'.inp' if the name has no '.').
 * CTP:     thLoad.c:getblock (block scanning), thParm.c:thLoadParameters/thParmLineSet
 *          (line parsing), thUtils.c:thTokenArray (array indices), daVar.h:floatToLong.
 *
 * Reproduced CTP behaviour:
 *  * ';' starts a comment (except inside '...' or "..." strings); tabs are spaces.
 *  * Only lines between "begin <type> <name>" and "end <type>" are read; the block
 *    name is irrelevant (the end line only has to name the type).
 *  * Variable names are case-insensitive. Unknown names are accepted and printed as
 *    "<name> not registered, registering as int" (once per name), exactly like CTP;
 *    their values do not reach any registered variable. This is how the obsolete
 *    '.' names (targ.A, spec.e.P, ...) are silently ignored by the current Fortran.
 *  * "name(i)" is 1-based, "name[i]" is 0-based; "name = a, b, c" fills consecutive
 *    elements; a line without '=' continues the previous array.
 *  * Values are evaluated (numbers and + - * / ( ) expressions); a real value stored
 *    in an integer variable is rounded with floatToLong ((int)(x+-0.5)).
 *  * Strings: first quoted token (or the raw remainder), truncated to character*N.
 *  * Fortran integer flags (*_int) become bool with dbase.f:convert_to_logical (>0).
 *  * Blocks of one file are all loaded before booking; within one file a repeated
 *    block name replaces the earlier block (same variable in CTP).
 */
#ifndef SIMC_INPUT_INPREADER_H
#define SIMC_INPUT_INPREADER_H

#include <map>
#include <string>
#include <vector>

#include "simc/core/SimcParameters.h"

namespace simc {

/// @addtogroup input
/// @{

/// One message of the reader (a CTP message or a reader note), with its source position.
struct InpMessage {
    std::string file;  ///< input file.
    int line{0};  ///< 1-based line number.
    std::string text;   ///< CTP-style message text (e.g. "targ.A not registered, registering as int")
};

/// @brief Reads .inp files with the CTP rules listed in the file description and fills a SimcParameters.
///
/// Holds a reference to the parameters (must outlive the reader). Not thread-safe.
class InpReader {
public:
    /// Reader options.
    struct Options {
        /// Directory used for extra_dbase_file ('infiles/' in the Fortran). Empty = directory of the main file.
        std::string infiles_dir;
        /// Follow extra_dbase_file after the main file (dbase.f:105-117).
        bool load_extra{true};
    };

    /// @param params parameters to fill (not owned).
    explicit InpReader(SimcParameters& params) : p_(params) {}
    /// @param params parameters to fill (not owned).
    /// @param opt reader options.
    InpReader(SimcParameters& params, Options opt) : p_(params), opt_(std::move(opt)) {}

    /// dbase.f:dbase_read: main file then the extra file. Returns false if a file cannot be opened or if the input
    /// is rejected (errors(): electron_arm/hadron_arm = 6, the retired SHMS LSA tune, D-032).
    /// @param main_file path of the main .inp file.
    /// @return true if all files were read and the input is accepted.
    bool loadRun(const std::string& main_file);
    /// Input errors found by loadRun (D-032 arm-code validation). Empty when the input is accepted.
    /// @return the error messages.
    const std::vector<std::string>& errors() const { return errors_; }

    /// CTP thload + thbook of a single file.
    /// @param path the .inp file.
    /// @return false if the file cannot be opened.
    bool loadFile(const std::string& path);

    /// Names reported as "not registered", in the order CTP prints them.
    /// @return the unregistered names.
    const std::vector<std::string>& notRegistered() const { return not_registered_; }
    /// Values assigned to unregistered names (CTP keeps them in auto-registered variables).
    /// @return name -> values.
    const std::map<std::string, std::vector<double>>& unregisteredValues() const { return unreg_values_; }
    /// Every message CTP would print, plus reader notes.
    /// @return the messages, in order.
    const std::vector<InpMessage>& messages() const { return messages_; }
    /// Files loaded, in order.
    /// @return the file paths.
    const std::vector<std::string>& filesLoaded() const { return files_; }
    /// For each registered parameter set by the input: "file:line" of the last assignment.
    /// @return parameter name -> "file:line".
    const std::map<std::string, std::string>& provenance() const { return provenance_; }

private:
    struct Block {
        std::string type;                 // "parm", "test", ...
        std::string name;
        int first_line{0};
        std::vector<std::pair<int, std::string>> lines;  // (line number, raw text) incl. begin/end
    };

    void bookParm(const Block& b, const std::string& file);
    void parmLine(std::string line, const std::string& file, int lineno);
    void note(const std::string& file, int line, const std::string& text);

    SimcParameters& p_;
    Options opt_;

    // CTP "current variable" state (thParm.c globals) used by continuation lines.
    struct Cursor {
        const ParamDescriptor* desc{nullptr};
        std::string unreg;   // name of the auto-registered variable, if not registered
        int index{0};
        bool valid{false};
    } cur_;

    std::vector<std::string> not_registered_;
    std::map<std::string, std::string> unreg_lower_;   // lower-case name -> spelling
    std::map<std::string, std::vector<double>> unreg_values_;
    std::vector<InpMessage> messages_;
    std::vector<std::string> errors_;
    std::vector<std::string> files_;
    std::map<std::string, std::string> provenance_;
};

/// Evaluate a CTP numeric argument (number or simple arithmetic expression).
/// Returns false if it cannot be evaluated. `is_int` is true for a plain integer token.
/// @param text the argument text.
/// @param[out] value the value.
/// @param[out] is_int true for a plain integer token.
/// @return false if the text cannot be evaluated.
bool EvalCtpNumber(const std::string& text, double& value, bool& is_int);

/// @}

}  // namespace simc

#endif  // SIMC_INPUT_INPREADER_H
