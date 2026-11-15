// src/core/SpectrometerId.cpp — see include/simc/core/SpectrometerId.h.
// The single table of spectrometer arm codes (the electron_arm/hadron_arm integers of the input), their names and
// which side of the beam they sit on. Used by the input validation (core/SimcParameters.cpp), the set-up (arm phi),
// the report writers and the spectrometer factory. The retired SHMS LSA code 6 is refused with an explanation
// (FORTRAN_DIFFERENCES.md D-032).
// Fortran: dbase.f:883-930 (codes and names), dbase.f:247-263 (arm phi).
#include "simc/core/SpectrometerId.h"

#include <stdexcept>

#include "simc/core/SimcParameters.h"

namespace simc {

namespace {

// The retired Fortran code (SHMS LSA tune, dbase.f:895-897). This file is the only place that knows it.
constexpr int kRetiredShmsLsaCode = 6;

const std::array<SpectrometerInfo, 7> kTable{{
    // id, code, key, label, stop_key, beam_right (beam right: phi = 3pi/2 in dbase.f:247-263, beam left: pi/2)
    {SpectrometerId::HMS, 1, "hms", "HMS", "hms", true},
    {SpectrometerId::SOS, 2, "sos", "SOS", "sos", false},
    {SpectrometerId::HRSR, 3, "hrsr", "HRS-R", "hrsr", true},
    {SpectrometerId::HRSL, 4, "hrsl", "HRS-L", "hrsl", false},
    {SpectrometerId::SHMS, 5, "shms", "SHMS", "shms", false},
    {SpectrometerId::CALO_HMS_SIDE, 7, "calo_hms", "calorimeter (HMS side of beam)", "calo", true},
    {SpectrometerId::CALO_SHMS_SIDE, 8, "calo_shms", "calorimeter (SHMS side of beam)", "calo", false},
}};

}  // namespace

const std::array<SpectrometerInfo, 7>& SpectrometerTable() { return kTable; }

std::optional<SpectrometerId> SpectrometerFromCode(int code) {
    for (const auto& r : kTable)
        if (r.code == code) return r.id;
    return std::nullopt;
}

std::optional<SpectrometerId> SpectrometerFromKey(const std::string& key) {
    for (const auto& r : kTable)
        if (key == r.key) return r.id;
    return std::nullopt;
}

const SpectrometerInfo& Info(SpectrometerId id) {
    for (const auto& r : kTable)
        if (r.id == id) return r;
    throw std::logic_error("SpectrometerId: value not in the table");
}

std::string ArmCodeError(int code, const char* which) {
    if (SpectrometerFromCode(code)) return "";
    const std::string name = std::string(which) + "_arm";
    if (code == kRetiredShmsLsaCode)
        return name + " = 6 was the retired SHMS LSA tune; use 5 = SHMS";
    return name + " = " + std::to_string(code) +
           " is not a spectrometer (1 = HMS, 2 = SOS, 3 = HRS-R, 4 = HRS-L, 5 = SHMS, 7/8 = calorimeter)";
}

// strict = false: only the retired code 6 is an error (unknown codes are left to the set-up, which reports them as the
// Fortran does); strict = true: every code outside the table is an error.
std::vector<std::string> ValidateArmCodes(const SimcParameters& p, bool strict) {
    std::vector<std::string> errors;
    const std::pair<int, const char*> arms[2] = {{p.electron_arm, "electron"}, {p.hadron_arm, "hadron"}};
    for (const auto& a : arms) {
        if (!strict && a.first != kRetiredShmsLsaCode) continue;
        const std::string e = ArmCodeError(a.first, a.second);
        if (!e.empty()) errors.push_back(e);
    }
    return errors;
}

}  // namespace simc
