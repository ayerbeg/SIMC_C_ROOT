// src/driver/SimcRun.cpp — see include/simc/driver/SimcRun.h.
// The run driver: the C++ counterpart of the simc.f main program. initialise() builds the RunSetup (setup/), resolves
// the pipeline (Pipeline.cpp: which reaction, physics model, radiator, target, arms and sinks), creates the stages by
// name through the Registry and computes the central event; step() is one try of the Fortran event loop (generate,
// montecarlo, complete_recon_ev, complete_main, pass_cuts, histograms/ntuple, counters); finish() normalises and
// hands the run summary to the sinks (.hist/.geni/.gen, summary JSON, ROOT tree). The parallel loop is in
// SimcRunParallel.cpp.
// Fortran: simc.f main program (l.1-623).
#include "simc/driver/SimcRun.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "simc/core/DataPath.h"
#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/RanluxRandom.h"
#include "simc/driver/FortranWriter.h"
#include "simc/driver/ReportSink.h"
#include "simc/material/FortranFormat.h"
#include "simc/setup/SetupBuilder.h"

namespace simc {
namespace driver {

namespace {

/// The driver's --quiet: every stage gets "verbose": false unless its config says otherwise (stage messages only;
/// no numerical effect).
Json withVerbose(Json cfg, bool quiet) {
    if (!cfg.is_object()) cfg = Json::object();
    if (quiet && !cfg.contains("verbose")) cfg["verbose"] = false;
    return cfg;
}

}  // namespace

// Fortran: simc.f l.226-238 (pass_cuts).
// FORTRAN-BUG(D-006): original Fortran simc.f:237: `recon%p%delta .ge. (SPedge%e%delta%max-slop%MC%p%delta%used)`;
// corrected because the hadron-arm delta must be tested against the hadron-arm limit SPedge%p%delta%max, like the
// other five hadron-arm tests (copy-paste of l.231). fortran_compat keeps the electron-arm limit.
bool PassCuts(const RunSetup& S, const Event& r) {
    const auto& E = S.SPedge;
    const auto& M = S.slop.MC;
    const double p_delta_max = S.fortran_compat ? E.e.delta.max : E.p.delta.max;
    return !(r.e.delta <= (E.e.delta.min + M.e.delta.used) || r.e.delta >= (E.e.delta.max - M.e.delta.used) ||
             r.e.yptar <= (E.e.yptar.min + M.e.yptar.used) || r.e.yptar >= (E.e.yptar.max - M.e.yptar.used) ||
             r.e.xptar <= (E.e.xptar.min + M.e.xptar.used) || r.e.xptar >= (E.e.xptar.max - M.e.xptar.used) ||
             r.p.delta <= (E.p.delta.min + M.p.delta.used) || r.p.delta >= (p_delta_max - M.p.delta.used) ||
             r.p.yptar <= (E.p.yptar.min + M.p.yptar.used) || r.p.yptar >= (E.p.yptar.max - M.p.yptar.used) ||
             r.p.xptar <= (E.p.xptar.min + M.p.xptar.used) || r.p.xptar >= (E.p.xptar.max - M.p.xptar.used));
}

SimcRun::SimcRun(SimcParameters params, const Json& pipeline, RunOptions options)
    : params_(std::move(params)), pipelineJson_(pipeline), opt_(std::move(options)) {
    if (opt_.override_ngen) params_.ngen = static_cast<int>(opt_.ngen);
    if (opt_.override_seed) params_.random_seed = static_cast<int>(opt_.seed);
    // dbase.f:119
    if (params_.random_seed == -1) params_.random_seed = static_cast<int>(std::time(nullptr));
}

SimcRun::~SimcRun() = default;

void SimcRun::setRandom(std::unique_ptr<IRandom> rng) {
    if (initialised_) throw std::logic_error("SimcRun::setRandom after initialise()");
    rng_ = std::move(rng);
    rngPtr_ = rng_.get();
}

IOutputSink& SimcRun::addSink(std::unique_ptr<IOutputSink> sink) {
    if (initialised_) throw std::logic_error("SimcRun::addSink after initialise()");
    if (!sink) throw std::invalid_argument("SimcRun::addSink: null sink");
    extraSinks_.push_back(std::move(sink));
    return *extraSinks_.back();
}

void SimcRun::log(const std::string& line) { (opt_.log ? *opt_.log : std::cout) << line << "\n"; }

// Fortran: simc.f l.53-160
void SimcRun::initialise() {
    if (initialised_) return;
    // ---- dbase_read, init.f (min_max_init, target_init, limits_init, radc_init), targetfac, luminosity.
    // target_init and extreme_trip_thru_target come from the pipeline's target model (ITargetLossProvider), through a
    // dedicated instance: the event-loop target below starts from its own fresh state (the target model may keep
    // per-event state that the set-up must not leave behind).
    {
        SetupOptions so;
        so.run_name = opt_.run_name;
        so.data_dir = opt_.data_dir.empty() ? DataDir("spectral") : opt_.data_dir;   // (data tree: core/DataPath.h)
        so.fortran_compat = opt_.fortran_compat;
        so.tree_layout = opt_.tree_layout;   // task U.2
        so.compat_target = opt_.compat_target;   // task U.R
        // D/A(e,e'p): a spectral function named in the JSON (D-033 F1) reaches the set-up, which reads the
        // theory file / chooses the Benhar file: "theory"/"benhar" set use_benhar_sf, "file" the file.
        const PipelineConfig::SpectralFunctionSpec sf = PipelineConfig::SpectralFunction(pipelineJson_);
        if (sf.type == "theory" || sf.type == "benhar") {
            const bool benhar = sf.type == "benhar";
            if (benhar != params_.use_benhar_sf)
                log(std::string(" NOTE: pipeline.reaction.spectral_function.type = ") + sf.type + " sets use_benhar_sf = " +
                    (benhar ? "T" : "F"));
            params_.use_benhar_sf = benhar;
        }
        so.sf_file = sf.file;
        std::unique_ptr<ITarget> startupTarget = Registry<ITarget>::create(PipelineConfig::TargetType(pipelineJson_));
        setup_ = BuildRunSetup(params_, so, startupTarget.get());
    }
    for (const auto& m : setup_.messages) log(" " + m);
    if (!setup_.ok()) {
        std::string all;
        for (const auto& s : setup_.stops) all += (all.empty() ? "" : "; ") + s;
        throw FortranStop("STOP during setup (legacy SIMC behaviour): " + all);
    }
    setup_.run_name = opt_.run_name;
    setup_.output_dir = opt_.output_dir;
    log(" Use random seed = " + std::to_string(params_.random_seed));
    log(setup_.fortran_compat ? std::string(" SIMC++: fortran_compat = true (exact Fortran behaviour of upstream ") +
                                    (setup_.compat294f662() ? "294f662" : "59ccdf6") + ", known Fortran bugs kept)"
                              : " SIMC++: fortran_compat = false (clear legacy-SIMC bugs corrected, docs/translation/markdown/06_catalogue.md)");
    if (!params_.random_state_file.empty() && params_.random_state_file.find_first_not_of(' ') != std::string::npos)
        log(" NOTE: random_state_file '" + params_.random_state_file +
            "' is not used by SIMC++ (the RANLUX state is not portable; D-002)");

    // ---- the pipeline
    pipe_ = PipelineConfig::Defaults(setup_);
    pipe_.apply(pipelineJson_);
    if (pipe_.reaction.type.empty())
        throw NotPorted("reaction '" + ReactionTitle(setup_.flags, params_) +
                        "': no event generator for these reaction flags and this target (no doing_* target class "
                        "is set, e.g. targ%A < 1, or (e,e'rho) with A >= 4)");
    const bool quiet = opt_.quiet_stages;

    // ---- the stages, created by name and wired through the StageContext (interfaces only)
    StageContext ctx(setup_);
    if (!rngPtr_) {
        // Mt64Random lives in the STATIC library simcfw_random: constructing it here also makes the linker keep its
        // object file, hence its SIMC_REGISTER(IRandom, "mt64") registration.
        if (pipe_.rng.type == "mt64" && !Registry<IRandom>::contains("mt64"))
            rng_ = std::make_unique<Mt64Random>();
        else if (pipe_.rng.type == "ranlux" && !Registry<IRandom>::contains("ranlux"))
            rng_ = std::make_unique<RanluxRandom>();   // (same reason: simcfw_random is STATIC)
        else
            rng_ = Registry<IRandom>::create(pipe_.rng.type);
        rng_->configure(pipe_.rng.config, params_);
        rng_->init(ctx);
        rngPtr_ = rng_.get();
    }
    target_ = Registry<ITarget>::create(pipe_.target.type);
    target_->configure(withVerbose(pipe_.target.config, quiet), params_);
    target_->init(ctx);

    radiator_ = Registry<IRadiator>::create(pipe_.radiator.type);
    radiator_->configure(withVerbose(pipe_.radiator.config, quiet), params_);
    radiator_->init(ctx);

    ctx.target = target_.get();
    ctx.radiator = radiator_.get();

    beam_ = Registry<IBeam>::create(pipe_.beam.type);
    beam_->configure(withVerbose(pipe_.beam.config, quiet), params_);
    beam_->init(ctx);

    // the physics model before the reaction: the reaction assembles the event weight with it (it calls the model's
    // cross section in completeMain), so the model must exist when the reaction's init() wires it in
    model_ = Registry<IPhysicsModel>::create(pipe_.physics.type);
    model_->configure(withVerbose(pipe_.physics.config, quiet), params_);
    model_->init(ctx);
    ctx.physics = model_.get();

    reaction_ = Registry<IReaction>::create(pipe_.reaction.type);
    reaction_->configure(withVerbose(pipe_.reaction.config, quiet), params_);
    reaction_->init(ctx);

    MonteCarlo::ArmChoice hadron{pipe_.hadron_arm.type, quiet ? withVerbose(pipe_.hadron_arm.config, true)
                                                              : pipe_.hadron_arm.config};
    MonteCarlo::ArmChoice electron{pipe_.electron_arm.type, quiet ? withVerbose(pipe_.electron_arm.config, true)
                                                                  : pipe_.electron_arm.config};
    mc_ = std::make_unique<MonteCarlo>(setup_, *target_, hadron, electron);
    if (quiet) mc_->setLog(nullptr);

    for (const auto& s : pipe_.sinks) {
        std::unique_ptr<IOutputSink> sink = Registry<IOutputSink>::create(s.type);
        sink->configure(s.config, params_);
        sink->init(ctx);
        sinks_.push_back(std::move(sink));
    }
    for (auto& s : extraSinks_) sinks_.push_back(std::move(s));
    extraSinks_.clear();

    // ---- simc.f:85 calculate_central (after radc_init; central%rad from the radiator, simc.f:1268-1281)
    central_ = calculateCentral(setup_, *reaction_, *radiator_, *model_);
    if (setup_.flags.doing_semi) {
        // simc.f:1289-1303 (calculate_central): the SIDIS central kinematics, printed list-directed as the Fortran
        // does (the numbers go to stdout only)
        const Event& v0 = central_.vertex0;
        auto item = [&](const char* label, double x) { log(" " + std::string(label) + " " + fortran::listReal(x)); };
        log(" SIDIS: central event kinematics");
        item("Pt", std::sqrt(v0.pt2) / 1.e3);
        item("z", v0.zhad);
        item("thpq", v0.theta_pq * 180.0 / 3.1415926536);
        item("Q2", v0.Q2);
        item("xbj", v0.xbj);
        item("Ein", v0.Ein);
        item("A", setup_.targ.A);
        item("Eprime", v0.e.E);
        item("Ppi", v0.p.P);
        item("The", v0.e.theta);
        item("multiplicity", central_.main0.ntup.sigcm);
        item("lab cross section (nb/Gev2/sr2)", central_.central.sigcc * 1000.0 * 1000.0 * 1000.0);
    }

    // simc.f:120-151 (using_tgt_field): the target-field / arm angles, printed before trgInit
    for (const auto& l : mc_->targetFieldLog()) log(l);

    // ---- report bookkeeping
    report_ = RunReport{};
    report_.run_name = opt_.run_name;
    report_.ngen = params_.ngen;
    report_.random_seed = params_.random_seed;
    report_.rng_type = rngPtr_ ? rngPtr_->name() : std::string();
    report_.random_state_file = params_.random_state_file;
    report_.central = central_.central;
    report_.totals.slop = setup_.slop;   // slop%*%used from limits_init, lo/hi accumulated by limits_update

    for (auto& s : sinks_) s->beginRun(setup_);

    // simc.f:160-166
    log(" ");
    report_.begin_time = FortranCtime(static_cast<long>(std::time(nullptr)));
    ev_ = EventState{};
    ninform_ = 1000;
    initialised_ = true;
}

bool SimcRun::done() const {
    if (!initialised_) return false;
    const long n = std::labs(static_cast<long>(params_.ngen));
    return report_.totals.nevent >= n;
}

// simc.f:174-183
void SimcRun::resetForTry() {
    NtupVars& nt = ev_.main.ntup;
    nt.Mh2_final = setup_.Mh2;   // Mh2_final = Mh2 (COMMON, carried in ntup)
    nt.radphot = 0.0;
    nt.radarm = 0.0;
    nt.decdist = 0.0;
    nt.resfac = 0.0;
    nt.sigcm = 0.0;
    nt.krel = 0.0;
    nt.mm = 0.0;
    nt.mmA = 0.0;
    nt.t = 0.0;
}

// Fortran: simc.f l.169-352 (one pass of the do while loop)
void SimcRun::step() {
    if (!initialised_) initialise();
    RunTotals& T = report_.totals;
    EventState& ev = ev_;
    MainEvent& main = ev.main;
    const RunSetup& S = setup_;
    IRandom& rng = *rngPtr_;

    resetForTry();
    ++T.ntried;
    ev.ntried = T.ntried;

    // --- 1. bookkeeping of the try
    // l.195-200 ... keep the human interested
    if (T.ntried % ninform_ == 1) {
        fortran::Record r;
        r.x(1).lit("Generating Event").i(T.ntried, 9).lit(" ... ").i(T.nevent, 8).lit("  successes so far - Monitor:");
        r.e(T.wtcontribute * S.luminosity / static_cast<double>(T.ntried), 11, 4);
        log(r.str());
        if (T.ntried >= 5000) ninform_ = 20000;
    }

    // --- 2. the event chain: generate the vertex, transport it through both arms, reconstruct, weight
    model_->beginEvent(T.nevent);
    reaction_->beginTry(T.ntried);
    // l.203 generate
    bool success = reactions::generateEvent(*beam_, *reaction_, main, ev.vertex, ev.orig, rng);
    // l.211 montecarlo
    if (success) success = mc_->run(ev.orig, main, ev.recon, rng);
    // l.216-218 complete_recon_ev
    if (success) success = reaction_->completeRecon(ev.recon);
    // l.224 complete_main
    if (success) success = reaction_->completeMain(false, main, ev.vertex, central_.vertex0, ev.recon);

    // --- 3. cuts: pass_cuts decides the resolution sums and wtcontr; hard_cuts also rejects the event
    // l.226-238 pass_cuts (PassCuts, top of this file; FORTRAN-BUG(D-006) there)
    ev.pass_cuts = PassCuts(S, ev.recon);
    // l.240-243
    if (params_.hard_cuts) {
        if (!ev.pass_cuts) success = false;
        if (S.flags.doing_eep && (ev.recon.Em > S.cuts.Em.max)) success = false;
    }
    ev.success = success;
    // l.245
    if (success) T.sum_sigcc = T.sum_sigcc + main.sigcc;

    // --- 4. output and counters
    // l.250-287 histograms (geni every try; RECON/gen for successes) and l.333 results_ntu_write: the sinks
    for (auto& s : sinks_) s->fill(main, ev.vertex, ev.orig, ev.recon, success);

    if (success) {
        // l.295-298
        T.ncontribute = T.ncontribute + 1;
        if (!radiator_->hadronRadiates()) T.ncontribute_no_rad_proton = T.ncontribute_no_rad_proton + 1;
        // l.311-329
        if (ev.pass_cuts) {
            const Event& r = ev.recon;
            const Event& v = ev.vertex;
            T.npasscuts = T.npasscuts + 1;
            T.wtcontribute = T.wtcontribute + main.weight;
            const double de = r.e.delta - main.SP.e.delta, xe = r.e.xptar - v.e.xptar, ye = r.e.yptar - v.e.yptar,
                         ze = r.e.z - main.SP.e.z;
            // FORTRAN-BUG(D-057 F18): original Fortran simc.f:316-317, 324-325 use vertex%p%xptar/yptar, the RHO's lab
            // theta/phi for doing_rho (the .hist 'Hadron arm' resolution lines); corrected: the decay pion (orig%p).
            const ArmFull& hp = GeneratedHadronArm(S.flags.doing_rho && !S.fortran_compat, v, ev.orig);
            const double dp = r.p.delta - main.SP.p.delta, xp = r.p.xptar - hp.xptar, yp = r.p.yptar - hp.yptar,
                         zp = r.p.z - main.SP.p.z;
            T.sumerr.e.delta = T.sumerr.e.delta + de;
            T.sumerr.e.xptar = T.sumerr.e.xptar + xe;
            T.sumerr.e.yptar = T.sumerr.e.yptar + ye;
            T.sumerr.e.ytar = T.sumerr.e.ytar + ze;
            T.sumerr.p.delta = T.sumerr.p.delta + dp;
            T.sumerr.p.xptar = T.sumerr.p.xptar + xp;
            T.sumerr.p.yptar = T.sumerr.p.yptar + yp;
            T.sumerr.p.ytar = T.sumerr.p.ytar + zp;
            // (x)**2 -> x*x
            T.sumerr2.e.delta = T.sumerr2.e.delta + de * de;
            T.sumerr2.e.xptar = T.sumerr2.e.xptar + xe * xe;
            T.sumerr2.e.yptar = T.sumerr2.e.yptar + ye * ye;
            T.sumerr2.e.ytar = T.sumerr2.e.ytar + ze * ze;
            T.sumerr2.p.delta = T.sumerr2.p.delta + dp * dp;
            T.sumerr2.p.xptar = T.sumerr2.p.xptar + xp * xp;
            T.sumerr2.p.yptar = T.sumerr2.p.yptar + yp * yp;
            T.sumerr2.p.ytar = T.sumerr2.p.ytar + zp * zp;
        }
        // l.332-333 limits_update (Egamma_used from the radiator)
        reactions::limitsUpdate(main, ev.vertex, ev.orig, ev.recon, S.flags, radiator_->photonEnergies(), T.contrib,
                                T.slop, S.flags.doing_rho && !S.fortran_compat);   // D-057 F18
    }

    // l.341-347: ngen > 0 counts successes, ngen < 0 counts tries (|ngen| tries)
    if (params_.ngen < 0) {
        T.nevent = T.nevent + 1;
    } else {
        if (success) T.nevent = T.nevent + 1;
    }
}

std::vector<StopCounters> SimcRun::collectStops() const {
    std::vector<StopCounters> v;
    ISpectrometer* h = &mc_->hadronArm();
    ISpectrometer* e = &mc_->electronArm();
    StopCounters sh;
    sh.spectrometer = h->name();
    sh.role = (h == e) ? "both" : "hadron";
    sh.counts = h->stopCounters();
    v.push_back(sh);
    if (e != h) {
        StopCounters se;
        se.spectrometer = e->name();
        if (se.spectrometer == sh.spectrometer) {   // two objects of one type (different configs): keep keys unique
            v.front().spectrometer += "_hadron";
            se.spectrometer += "_electron";
        }
        se.role = "electron";
        se.counts = e->stopCounters();
        v.push_back(se);
    }
    return v;
}

// Fortran: simc.f l.354-620
void SimcRun::finish() {
    if (finished_) return;
    if (!initialised_) initialise();
    RunTotals& T = report_.totals;
    {
        fortran::Record r;
        r.x(1).lit("---> Last Event ").i(T.ntried, 9).lit(" ...").i(T.nevent, 9).lit("  successes");
        log(r.str());
    }
    report_.end_time = FortranCtime(static_cast<long>(std::time(nullptr)));
    report_.norm = Normalise(setup_, T, reaction_->genVolume());
    report_.stops = mergedStopsSet_ ? mergedStops_ : collectStops();
    finished_ = true;

    // l.446: if (ngen.eq.0) goto 1000 — no output files
    if (params_.ngen == 0) {
        log(" ngen = 0: no output (simc.f:446)");
        return;
    }
    // The report prints the COMMON state as it is at the end of the run (simc.f:844 targ%Mrec_struck, overwritten per
    // event for the pion Delta final states): the reaction may update its part of it (IReaction::finalTargetState).
    RunSetup reportSetup = setup_;
    (finalStateReaction_ ? finalStateReaction_ : reaction_.get())->finalTargetState(reportSetup.targ);
    const RunSummary summary = MakeRunSummary(reportSetup, report_);
    for (auto& s : sinks_) {
        s->endRun(summary);
        if (auto* rs = dynamic_cast<IReportSink*>(s.get())) rs->writeReport(reportSetup, report_);
    }
    log(" Come back soon!");
}

// The whole run: initialise, the event loop (sequential, or the parallel loop of SimcRunParallel.cpp), finish.
void SimcRun::run() {
    initialise();
    const int n = threadCount();
    if (n > 1) {
        runThreaded(n);
        return;
    }
    while (!done()) step();
    finish();
}

}  // namespace driver
}  // namespace simc
