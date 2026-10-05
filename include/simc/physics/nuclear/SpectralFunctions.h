// include/simc/physics/nuclear/SpectralFunctions.h
//
// ISpectralFunction plug-ins (task P8.1), components of the D(e,e'p) / A(e,e'p) reaction plug-ins:
//
//   "theory"  TheorySpectralFunction: the init.f theory files (h2.theory for D, c12/fe56/au197.theory for A > 2), read
//             at start-up by BuildRunSetup (init.f:theory_init -> RunSetup::theory); weight = event.f:complete_main
//             l.1402-1430 (nuclear::TheorySfWeight).
//   "benhar"  BenharSpectralFunction: a Benhar table (benharsf_*.dat, use_benhar_sf = 1), read in init() by
//             sf_lookup_init (nuclear::SfTable::load); weight = targ%Z * transparency * sf_lookup_diff(Em, Pm)
//             (event.f:1399-1401).
//
// JSON (the "spectral_function" object of the reaction config; every key optional):
//   {"type": "theory" | "benhar", "file": "<file in the data directory>", "transparency": <Benhar only>}
// The file defaults to the Fortran choice from the target A (RunSetup::sf; a warning is printed when the Fortran
// falls back to the C-12 file for an A without its own file, D-033 F1). For "theory" the setup has already read the
// file (it also fixes the VERTEXedge limits), so a "file" given here must be the same one: SimcRun passes it to the
// setup (SetupOptions::sf_file) before the reaction is built.
#ifndef SIMC_PHYSICS_NUCLEAR_SPECTRALFUNCTIONS_H
#define SIMC_PHYSICS_NUCLEAR_SPECTRALFUNCTIONS_H

#include <iosfwd>
#include <string>

#include "simc/core/RunSetup.h"
#include "simc/physics/nuclear/SpectralFunctionTables.h"
#include "simc/stages/ISpectralFunction.h"

namespace simc {
namespace nuclear {

/// "theory": event.f:complete_main l.1402-1430 on RunSetup::theory (init.f:theory_init).
class TheorySpectralFunction : public ISpectralFunction {
public:
    std::string name() const override { return "theory"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    double weight(const Event& vertex) override;
    std::string file() const override { return t_ ? t_->file : file_; }

private:
    std::string file_;
    const TheoryData* t_{nullptr};
    bool heavy_{false};
};

/// "benhar": event.f:1399-1401, sf_lookup.f (the table and its SAVEd state live in this object).
class BenharSpectralFunction : public ISpectralFunction {
public:
    std::string name() const override { return "benhar"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    double weight(const Event& vertex) override;
    std::string file() const override { return file_; }

    /// The sf_lookup kernel (L1 tests).
    SfTable& table() { return table_; }
    double transparency() const { return transparency_; }
    void setLog(std::ostream* os) { log_ = os; }

private:
    std::string file_;
    double transparency_{0}, Z_{0};
    SfTable table_;
    std::ostream* log_;

public:
    BenharSpectralFunction();
};

}  // namespace nuclear
}  // namespace simc

#endif  // SIMC_PHYSICS_NUCLEAR_SPECTRALFUNCTIONS_H
