// tests/unit/test_inp_reader.cpp — CTP parsing semantics on a synthetic .inp
// (thParm.c:thParmLineSet, thLoad.c:getblock, dbase.f:dbase_read extra file).
#include <cstdio>
#include <fstream>
#include <string>

#include "TestUtil.h"
#include "simc/input/EffectiveValues.h"
#include "simc/input/InpReader.h"

using namespace simc;
using simctest::check;

int main() {
    const std::string dir = "inp_reader_tmp";
    std::system(("mkdir -p " + dir).c_str());
    {
        std::ofstream f(dir + "/main.inp");
        f << "; comment line\n"
             "Ebeam = 1.0   ; outside any block: ignored by CTP\n"
             "begin parm experiment\n"
             "  EBEAM = 3056.\t\t; case-insensitive, tab separated\n"
             "  ngen = -1000\n"
             "  EXPER%charge = 2*0.5        ; expression\n"
             "  extra_dbase_file='extra_x' ;  quoted string\n"
             "  random_state_file = 'a;b'  ; ';' inside a string is not a comment\n"
             "  doing_pion = 2            ; flag: >0 => .true.\n"
             "  mc_smear = -1             ; flag: <=0 => .false.\n"
             "end parm whatever_name_is_irrelevant\n"
             "begin parm target\n"
             "  targ%fr_pattern = 1.6     ; real into int: floatToLong -> 2\n"
             "  targ%can = -1.5           ; -> -2\n"
             "  targ.A = 12               ; obsolete name\n"
             "  TARG.a = 13               ; same unknown name, different case: reported once\n"
             "  foo = 1, 2.5\n"
             "end parm target\n"
             "begin parm debug\n"
             "  debug(2) = 1              ; 1-based\n"
             "  debug[4] = 1              ; 0-based -> debug(5)\n"
             "  lambda = 0.1,\n"
             "   0.2, 0.3                 ; continuation line fills lambda(2),lambda(3)\n"
             "  lambda(4) = 9             ; past the end: error, ignored\n"
             "end parm debug\n";
    }
    {
        std::ofstream f(dir + "/extra_x.inp");
        f << "begin parm simulate\n"
             "  Ebeam = 4000.   ; the extra file overrides the main one\n"
             "  using_rad = 1\n"
             "  cuts%Em%min = 0.\n"
             "  cuts%Em%max = 0.\n"
             "  extra_dbase_file = 'never_loaded'\n"
             "end parm\n";
    }

    SimcParameters p;
    InpReader r(p);
    check(r.loadRun(dir + "/main.inp"), "loadRun succeeds");
    check(r.filesLoaded().size() == 2, "main + extra file loaded (extra not followed twice)");
    check(p.Ebeam == 4000.0, "extra file overrides Ebeam");
    check(p.ngen == -1000, "ngen");
    check(p.EXPER.charge == 1.0, "expression 2*0.5");
    check(p.extra_dbase_file == "never_loaded", "extra_dbase_file keeps the last assignment");
    check(p.random_state_file == "a;b", "string with ';'");
    check(p.doing_pion, "flag 2 -> true");
    check(!p.mc_smear, "flag -1 -> false");
    check(p.targ.fr_pattern == 2, "floatToLong(1.6) = 2");
    check(p.targ.can == -2, "floatToLong(-1.5) = -2");
    check(p.targ.A == 0.0, "obsolete name does not reach targ%A");
    check(p.debug[1] && p.debug[4] && !p.debug[0], "debug(2) and debug[4]");
    check(p.lambda[0] == 0.1 && p.lambda[1] == 0.2 && p.lambda[2] == 0.3, "lambda with continuation");
    check(p.using_rad, "using_rad from extra");
    check(p.cuts.Em.min == 0.0 && p.cuts.Em.max == 0.0, "cuts%Em from extra");

    check(r.notRegistered().size() == 2, "two unknown names reported (targ.A, foo)");
    if (r.notRegistered().size() == 2) {
        check(r.notRegistered()[0] == "targ.A", "first unknown spelled as in the file");
        check(r.notRegistered()[1] == "foo", "second unknown");
    }
    bool past_end = false, msg_ok = false;
    for (const auto& m : r.messages()) {
        if (m.text.find("Tried to fill past end of array lambda") != std::string::npos) past_end = true;
        if (m.text == "targ.A not registered, registering as int") msg_ok = true;
    }
    check(past_end, "lambda(4) reported as past the end");
    check(msg_ok, "CTP wording of the not-registered message");
    check(r.unregisteredValues().at("foo").size() == 2, "values of unregistered 'foo' kept");

    // Effective-value rules
    auto rows = EffectiveValues(p);
    std::string em_min, dE;
    for (const auto& row : rows) {
        if (row.name == "cuts%Em%min") em_min = row.effective;
        if (row.name == "dEbeam") dE = row.effective;
    }
    check(em_min == "-1000000", "cuts%Em opened when min=max (dbase.f:534-538), got " + em_min);
    check(dE == "0", "dEbeam effective in MeV");
    return simctest::finish("inp_reader_ctp_semantics");
}
