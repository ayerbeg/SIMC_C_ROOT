// include/simc/input/EffectiveValues.h
//
// Table of the values the Fortran actually uses after dbase.f has converted the raw
// .inp values (units, convert_to_logical, forced settings). It covers the per-variable
// rules of dbase.f:dbase_read and the two input-only rules of init.f:radc_init.
// It does NOT include quantities derived from several inputs at run time (edges and
// slops of limits_init, target_init averages, central event...) — those belong to RunSetup.
#ifndef SIMC_INPUT_EFFECTIVEVALUES_H
#define SIMC_INPUT_EFFECTIVEVALUES_H

#include <string>
#include <vector>

#include "simc/core/SimcParameters.h"

namespace simc {

struct EffectiveRow {
    std::string name;       ///< Fortran name
    std::string raw;        ///< value as read from the input (input units)
    std::string effective;  ///< value the Fortran uses after dbase.f/init.f conversions
    std::string units;      ///< units of `effective`
    std::string rule;       ///< Fortran rule applied (file:line), empty if none
};

std::vector<EffectiveRow> EffectiveValues(const SimcParameters& p);

/// Fixed-width text rendering of the table.
std::string FormatEffectiveTable(const std::vector<EffectiveRow>& rows);

}  // namespace simc

#endif  // SIMC_INPUT_EFFECTIVEVALUES_H
