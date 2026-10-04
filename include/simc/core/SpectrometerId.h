// include/simc/core/SpectrometerId.h
//
// ONE number per spectrometer (D-032, task P7.5). Every piece of SIMC++ that needs to know which spectrometer an arm
// code (.inp/JSON `electron_arm` / `hadron_arm`) means goes through this header: the input readers (validation), the
// set-up (arm phi, slop parameters), the target windows, the ntuple arm side, the default spectrometer plug-in, the
// .geni blocks and the run summary. No other file compares an arm code with a literal number.
//
// Fortran origin: dbase.f:883-930 (the announcement), infiles/inputfile_protoype.inp:22-23,
//   1 = HMS, 2 = SOS, 3 = HRS-R, 4 = HRS-L, 5 = SHMS (SSA tune), 7 = calorimeter on the HMS side of the beam,
//   8 = calorimeter on the SHMS side.
// The Fortran code 6 (SHMS with the LSA tune) is RETIRED: dbase.f:895-897/915-917 print "This spectrometer option no
// longer available" and STOP, and mc_shms.f:335-340 STOPs "LSA tune no longer in use!", although dbase.f:250/259,
// init.f:176/199, target.f:97/217 and results_write.f:90 still treat 6 as SHMS. SIMC++ rejects 6 when the input is
// read (InpReader::loadRun, FromJson, and BuildRunSetup for programmatic parameters) with the message of
// ArmCodeError(); 6 is not a SpectrometerId and is known only to the table in SpectrometerId.cpp.
// docs/FORTRAN_DIFFERENCES.md (item D-032).
#ifndef SIMC_CORE_SPECTROMETERID_H
#define SIMC_CORE_SPECTROMETERID_H

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace simc {

struct SimcParameters;

/// The spectrometer numbers of the Fortran input (electron_arm / hadron_arm). 6 is deliberately absent (retired).
enum class SpectrometerId : int {
    HMS = 1,
    SOS = 2,
    HRSR = 3,
    HRSL = 4,
    SHMS = 5,
    CALO_HMS_SIDE = 7,
    CALO_SHMS_SIDE = 8,
};

/// One row of the name <-> number table.
struct SpectrometerInfo {
    SpectrometerId id;
    int code;                 ///< the .inp/JSON number
    const char* key;          ///< lower-case key used in the summary JSON / .geni maps / plug-in names ("hms", ...)
    const char* label;        ///< human-readable name ("HMS", "HRS-R", "calorimeter (HMS side)")
    const char* stop_key;     ///< key of its stop counters in the summary/.geni ("calo" for both calorimeter sides:
                              ///< the Fortran has one set of caloSTOP_* counters, simc.f:527-531)
    bool beam_right;          ///< true: phi = 3pi/2 (dbase.f:247-263: HMS, HRS-R, calo HMS side); false: phi = pi/2
};

/// The table (all seven spectrometers, in code order).
const std::array<SpectrometerInfo, 7>& SpectrometerTable();

/// The spectrometer of an arm code, or nullopt for any code that is not in the table (0, 6, 9, ...).
std::optional<SpectrometerId> SpectrometerFromCode(int code);
/// The spectrometer of a key ("hms", "shms", ...), case-sensitive; nullopt if unknown.
std::optional<SpectrometerId> SpectrometerFromKey(const std::string& key);

const SpectrometerInfo& Info(SpectrometerId id);
inline int Code(SpectrometerId id) { return static_cast<int>(id); }
inline const char* Key(SpectrometerId id) { return Info(id).key; }
inline const char* Label(SpectrometerId id) { return Info(id).label; }
inline bool IsBeamRight(SpectrometerId id) { return Info(id).beam_right; }
inline const char* StopKey(SpectrometerId id) { return Info(id).stop_key; }
/// True if `code` is the arm code of `id`.
inline bool Is(int code, SpectrometerId id) { return code == Code(id); }

/// "" if `code` is a valid arm code; otherwise the message the input readers report. For the retired code 6:
/// "<which>_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS". `which` is "electron" or "hadron".
std::string ArmCodeError(int code, const char* which);
/// ArmCodeError for electron_arm and hadron_arm of a parameter set (empty = both valid).
/// `strict` = false accepts the unknown codes other than the retired one (they reach the Fortran set-up STOP
/// "I dont know what phi should be for the ... arm", dbase.f:252/261, which SIMC++ reproduces); only the retired code
/// is an input error then. The input readers use strict = false.
std::vector<std::string> ValidateArmCodes(const SimcParameters& p, bool strict = false);

}  // namespace simc

#endif  // SIMC_CORE_SPECTROMETERID_H
