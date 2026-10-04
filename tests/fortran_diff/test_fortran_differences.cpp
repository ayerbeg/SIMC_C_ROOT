// tests/fortran_diff/test_fortran_differences.cpp — task P7.5 (D-031, D-032).
//
// One sub-test per item of docs/FORTRAN_DIFFERENCES.md that changes behaviour. For each corrected clear Fortran bug
// (class a) the test runs the C++ code twice: with fortran_compat = true it must give the Fortran result (taken from
// the Fortran source/reference where possible), with the default (corrected) the corrected result.
//
//   usage: test_fortran_differences <item> <reference dir> <simc_gfortran dir>
// Exit 0 pass, 1 fail, 77 skip (reference data missing).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "simc/core/Json.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/SpectrometerId.h"
#include "simc/input/InpReader.h"
#include "simc/setup/SetupBuilder.h"

using namespace simc;

namespace {

int g_checks = 0, g_fail = 0;
void expect(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) ++g_fail;
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what.c_str());
}

std::string g_ref, g_fortran;

// ------------------------------------------------------------------------------------------------ D-032
// One number per spectrometer; the retired code 6 is rejected at input (.inp and JSON) and by the set-up.
int testSpectrometerId() {
    const auto& t = SpectrometerTable();
    expect(t.size() == 7, "table has the 7 spectrometers");
    for (const auto& r : t) {
        expect(SpectrometerFromCode(r.code) == r.id && Code(r.id) == r.code, std::string("code round trip ") + r.key);
        expect(SpectrometerFromKey(r.key) == r.id, std::string("key round trip ") + r.key);
        expect(ArmCodeError(r.code, "electron").empty(), std::string("code accepted ") + r.key);
    }
    // dbase.f:247-263: phi = 3pi/2 for HMS, HRS-R, calo HMS side
    expect(IsBeamRight(SpectrometerId::HMS) && IsBeamRight(SpectrometerId::HRSR) &&
               IsBeamRight(SpectrometerId::CALO_HMS_SIDE) && !IsBeamRight(SpectrometerId::SHMS) &&
               !IsBeamRight(SpectrometerId::SOS) && !IsBeamRight(SpectrometerId::HRSL) &&
               !IsBeamRight(SpectrometerId::CALO_SHMS_SIDE),
           "beam side of every spectrometer as dbase.f:247-263");
    expect(!SpectrometerFromCode(6) && !SpectrometerFromCode(0) && !SpectrometerFromCode(9), "6, 0, 9 are not ids");
    expect(ArmCodeError(6, "electron") == "electron_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS",
           "message for electron_arm = 6");
    expect(ArmCodeError(6, "hadron") == "hadron_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS",
           "message for hadron_arm = 6");

    // JSON input
    for (const char* which : {"electron_arm", "hadron_arm"}) {
        SimcParameters p;
        Json j = {{"electron_arm", 5}, {"hadron_arm", 1}};
        j[which] = 6;
        std::string what;
        try {
            FromJson(j, p);
        } catch (const std::invalid_argument& e) {
            what = e.what();
        }
        expect(what.find(std::string(which) + " = 6 was the retired SHMS LSA tune; use 5 = SHMS") != std::string::npos,
               std::string("FromJson rejects ") + which + " = 6: '" + what + "'");
        SimcParameters ok;
        bool threw = false;
        try {
            FromJson(Json{{"electron_arm", 5}, {"hadron_arm", 1}}, ok);
        } catch (const std::exception&) {
            threw = true;
        }
        expect(!threw && ok.electron_arm == 5 && ok.hadron_arm == 1, "FromJson accepts 5/1");
    }

    // .inp input
    const std::string dir = "fortran_diff_d032";
    std::system(("mkdir -p " + dir).c_str());
    for (const char* which : {"electron_arm", "hadron_arm"}) {
        const std::string f = dir + "/" + which + "6.inp";
        {
            std::ofstream o(f);
            o << "begin parm x\n  electron_arm = " << (std::string(which) == "electron_arm" ? 6 : 5)
              << "\n  hadron_arm = " << (std::string(which) == "hadron_arm" ? 6 : 1) << "\nend parm x\n";
        }
        SimcParameters p;
        InpReader rd(p);
        const bool ok = rd.loadRun(f);
        expect(!ok && rd.errors().size() == 1 &&
                   rd.errors()[0] == std::string(which) + " = 6 was the retired SHMS LSA tune; use 5 = SHMS",
               std::string("InpReader rejects ") + which + " = 6");
    }
    // set-up from parameters built in code (both modes)
    for (bool compat : {false, true}) {
        SimcParameters p;
        p.electron_arm = 5;
        p.hadron_arm = 6;
        SetupOptions so;
        so.fortran_compat = compat;
        const RunSetup s = BuildRunSetup(p, so);
        expect(!s.ok() && s.stops.front().find("hadron_arm = 6 was the retired SHMS LSA tune") != std::string::npos,
               std::string("BuildRunSetup stops for hadron_arm = 6, fortran_compat=") + (compat ? "true" : "false"));
    }
    return 0;
}

const std::map<std::string, std::function<int()>>& items() {
    static const std::map<std::string, std::function<int()>> m = {
        {"D032_spectrometer_id", testSpectrometerId},
    };
    return m;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <item> [reference dir] [simc_gfortran dir]\n", argv[0]);
        return 2;
    }
    g_ref = argc > 2 ? argv[2] : "";
    g_fortran = argc > 3 ? argv[3] : "";
    const auto it = items().find(argv[1]);
    if (it == items().end()) {
        std::fprintf(stderr, "unknown item %s\n", argv[1]);
        return 2;
    }
    const int rc = it->second();
    if (rc == 77) return 77;
    std::printf("%s: %d checks, %d failures\n", argv[1], g_checks, g_fail);
    return g_fail ? 1 : 0;
}
