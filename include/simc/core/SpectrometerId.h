/**
 * @file
 * @ingroup core
 * @brief simc::SpectrometerId: ONE number per spectrometer, and the table that maps the arm codes of the input to
 *        names, keys and properties.
 *
 * Every piece of SIMC++ that needs to know which spectrometer an arm code (`.inp`/JSON `electron_arm` /
 * `hadron_arm`) means goes through this header: the input readers (validation), the set-up (arm phi, slop parameters),
 * the target windows, the ntuple arm side, the default spectrometer plug-in, the `.geni` blocks and the run summary.
 * No other file compares an arm code with a literal number.
 *
 * Codes: 1 = HMS, 2 = SOS, 3 = HRS-R, 4 = HRS-L, 5 = SHMS (SSA tune), 7 = calorimeter on the HMS side of the beam,
 * 8 = calorimeter on the SHMS side. Code 6 (SHMS with the LSA tune) is RETIRED: the original program stops for it,
 * although parts of it still treat 6 as SHMS. SIMC++ rejects 6 when the input is read (InpReader::loadRun,
 * simc::FromJson, and simc::BuildRunSetup for programmatic parameters) with the message of simc::ArmCodeError();
 * 6 is not a SpectrometerId. Decision D-032 (docs/FORTRAN_DIFFERENCES.md), task P7.5.
 * @note Origin: Fortran SIMC dbase.f:883-930, infiles/inputfile_protoype.inp:22-23; code 6: dbase.f:895-897/915-917,
 *       mc_shms.f:335-340.
 */
#ifndef SIMC_CORE_SPECTROMETERID_H
#define SIMC_CORE_SPECTROMETERID_H

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace simc {

/** @addtogroup core
 *  @{ */

struct SimcParameters;

/// The spectrometer numbers of the Fortran input (electron_arm / hadron_arm). 6 is deliberately absent (retired).
enum class SpectrometerId : int {
    HMS = 1,              ///< Hall C High Momentum Spectrometer (beam right).
    SOS = 2,              ///< Hall C Short Orbit Spectrometer (beam left).
    HRSR = 3,             ///< Hall A right High Resolution Spectrometer (beam right).
    HRSL = 4,             ///< Hall A left High Resolution Spectrometer (beam left).
    SHMS = 5,             ///< Hall C Super High Momentum Spectrometer, SSA tune (beam left).
    CALO_HMS_SIDE = 7,    ///< Calorimeter (e.g. NPS) on the HMS side of the beam (beam right).
    CALO_SHMS_SIDE = 8,   ///< Calorimeter on the SHMS side of the beam (beam left).
};

/// One row of the name <-> number table.
struct SpectrometerInfo {
    SpectrometerId id;        ///< The spectrometer.
    int code;                 ///< the .inp/JSON number
    const char* key;          ///< lower-case key used in the summary JSON / .geni maps / plug-in names ("hms", ...)
    const char* label;        ///< human-readable name ("HMS", "HRS-R", "calorimeter (HMS side)")
    const char* stop_key;     ///< key of its stop counters in the summary/.geni ("calo" for both calorimeter sides:
                              ///< the Fortran has one set of caloSTOP_* counters, simc.f:527-531)
    bool beam_right;          ///< true: phi = 3pi/2 (dbase.f:247-263: HMS, HRS-R, calo HMS side); false: phi = pi/2
};

/**
 * @brief The spectrometer table.
 * @return All seven spectrometers, in code order (a static table; thread-safe).
 */
const std::array<SpectrometerInfo, 7>& SpectrometerTable();

/**
 * @brief The spectrometer of an arm code.
 * @param code The `.inp`/JSON arm number.
 * @return The spectrometer, or std::nullopt for any code that is not in the table (0, 6, 9, ...).
 */
std::optional<SpectrometerId> SpectrometerFromCode(int code);
/**
 * @brief The spectrometer of a key.
 * @param key Lower-case key, e.g. "hms", "shms" (case-sensitive).
 * @return The spectrometer, or std::nullopt if the key is unknown.
 */
std::optional<SpectrometerId> SpectrometerFromKey(const std::string& key);

/**
 * @brief The table row of a spectrometer.
 * @param id The spectrometer.
 * @return Its row of simc::SpectrometerTable() (static storage).
 */
const SpectrometerInfo& Info(SpectrometerId id);
/**
 * @brief The arm code of a spectrometer.
 * @param id The spectrometer.
 * @return The `.inp`/JSON number (1..5, 7, 8).
 */
inline int Code(SpectrometerId id) { return static_cast<int>(id); }
/**
 * @brief The lower-case key of a spectrometer.
 * @param id The spectrometer.
 * @return The key used in the summary JSON, the `.geni` maps and as plug-in name ("hms", ...); static string.
 */
inline const char* Key(SpectrometerId id) { return Info(id).key; }
/**
 * @brief The human-readable name of a spectrometer.
 * @param id The spectrometer.
 * @return E.g. "HMS", "HRS-R", "calorimeter (HMS side)"; static string.
 */
inline const char* Label(SpectrometerId id) { return Info(id).label; }
/**
 * @brief Tells on which side of the beam a spectrometer stands.
 * @param id The spectrometer.
 * @return True for beam right (central phi = 3 pi/2: HMS, HRS-R, calorimeter on the HMS side), false for beam left
 *         (phi = pi/2).
 */
inline bool IsBeamRight(SpectrometerId id) { return Info(id).beam_right; }
/**
 * @brief The key of the stop counters of a spectrometer.
 * @param id The spectrometer.
 * @return The key in the summary and `.geni` ("calo" for both calorimeter sides); static string.
 */
inline const char* StopKey(SpectrometerId id) { return Info(id).stop_key; }
/**
 * @brief Compares an arm code with a spectrometer.
 * @param code An arm code of the input.
 * @param id The spectrometer.
 * @return True if @p code is the arm code of @p id.
 */
inline bool Is(int code, SpectrometerId id) { return code == Code(id); }

/**
 * @brief The input error message for an arm code, if any.
 * @param code The arm code.
 * @param which "electron" or "hadron" (used in the message).
 * @return "" if @p code is a valid arm code; otherwise the message the input readers report. For the retired code 6:
 *         "<which>_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS".
 */
std::string ArmCodeError(int code, const char* which);
/**
 * @brief Checks electron_arm and hadron_arm of a parameter set.
 * @param p The parameters.
 * @param strict If false (what the input readers use), only the retired code 6 is an error: other unknown codes are
 *        accepted here and reach the set-up stop "I dont know what phi should be for the ... arm" (dbase.f:252/261),
 *        which SIMC++ reproduces. If true, every code outside the table is an error.
 * @return The simc::ArmCodeError messages; empty if both codes are valid.
 */
std::vector<std::string> ValidateArmCodes(const SimcParameters& p, bool strict = false);

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_SPECTROMETERID_H
