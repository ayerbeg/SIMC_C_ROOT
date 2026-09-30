// include/simc/core/RunSetup.h
//
// RunSetup: quantities derived once from SimcParameters (dbase.f / init.f equivalents).
// This is the P1.4 subset needed by the output sinks; the SetupBuilder of later tasks
// completes it (edges, slops, target_init, radc_init, central event).
#ifndef SIMC_CORE_RUNSETUP_H
#define SIMC_CORE_RUNSETUP_H

#include <string>
#include <vector>

#include "simc/core/SimcParameters.h"

namespace simc {

/// Reaction-family flags. Fortran: dbase.f:dbase_read l.124-234 (doing_* logicals).
struct ReactionFlags {
    bool doing_eep{false}, doing_hyd_elast{false}, doing_deuterium{false}, doing_heavy{false};
    bool doing_pion{false}, doing_hydpi{false}, doing_deutpi{false}, doing_hepi{false};
    bool doing_kaon{false}, doing_hydkaon{false}, doing_deutkaon{false}, doing_hekaon{false};
    bool doing_delta{false}, doing_hyddelta{false}, doing_deutdelta{false}, doing_hedelta{false};
    bool doing_semi{false}, doing_semipi{false}, doing_semika{false};
    bool doing_hydsemi{false}, doing_deutsemi{false}, doing_hesemi{false};
    bool doing_rho{false}, doing_hydrho{false}, doing_deutrho{false}, doing_herho{false};
    bool doing_phsp{false}, doing_pizero{false}, doing_hplus{false}, do_fermi{false};
    bool using_tgt_field{false};
};

/// Fortran: dbase.f l.124-234 (without the STOP checks, which are reported in `problems`).
ReactionFlags DeriveReactionFlags(const SimcParameters& p, std::vector<std::string>* problems = nullptr);

/// Fortran: simc.f:report l.664-790 — the reaction title line printed in the .hist.
std::string ReactionTitle(const ReactionFlags& f, const SimcParameters& p);

/// One 1-D histogram axis of the Fortran diagnostic histograms (histograms.inc: min, bin; nHbins=50).
struct HistAxis {
    double min{0};
    double bin{1};
};

/// Axes of H%gen/H%geni/H%RECON (all three share them). Fortran: init.f:limits_init l.515-570.
struct DiagHistAxes {
    struct Arm { HistAxis delta, yptar, xptar; };
    Arm e, p;
    HistAxis Em, Pm;
};

struct RunSetup {
    SimcParameters params;          ///< raw input (.inp units)
    ReactionFlags flags;
    std::string run_name{"simc"};   ///< Fortran `base` (file stem)
    std::string output_dir{"."};    ///< Fortran 'outfiles/' and 'worksim/'

    // dbase.f:239-244
    double e_theta{0}, e_sin_th{0}, e_cos_th{0};   ///< rad
    double p_theta{0}, p_sin_th{0}, p_cos_th{0};

    DiagHistAxes hist;              ///< filled by limits_init equivalent (SetupBuilder)

    /// Fill params, flags and the spectrometer angles (direct dbase.f conversions only).
    static RunSetup FromParameters(const SimcParameters& p, const std::string& run_name = "simc");
};

}  // namespace simc

#endif  // SIMC_CORE_RUNSETUP_H
