// src/driver/SimcRun.cpp — see include/simc/driver/SimcRun.h.
// Fortran: simc.f main program (l.1-623).
#include "simc/driver/SimcRun.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "simc/core/Registry.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/driver/FortranWriter.h"
#include "simc/driver/ReportSink.h"
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
    // dedicated instance: the event-loop target below starts from its own fresh state, as before P7.4.
    {
        SetupOptions so;
        so.run_name = opt_.run_name;
        so.data_dir = opt_.data_dir;
        std::unique_ptr<ITarget> startupTarget = Registry<ITarget>::create(PipelineConfig::TargetType(pipelineJson_));
        setup_ = BuildRunSetup(params_, so, startupTarget.get());
    }
    for (const auto& m : setup_.messages) log(" " + m);
    if (!setup_.ok()) {
        std::string all;
        for (const auto& s : setup_.stops) all += (all.empty() ? "" : "; ") + s;
        throw std::runtime_error("Fortran STOP during the setup: " + all);
    }
    setup_.run_name = opt_.run_name;
    setup_.output_dir = opt_.output_dir;
    log(" Use random seed = " + std::to_string(params_.random_seed));
    if (!params_.random_state_file.empty() && params_.random_state_file.find_first_not_of(' ') != std::string::npos)
        log(" NOTE: random_state_file '" + params_.random_state_file +
            "' is not used by SIMC++ (the RANLUX state is not portable; D-002)");

    // ---- the pipeline
    pipe_ = PipelineConfig::Defaults(setup_);
    pipe_.apply(pipelineJson_);
    if (pipe_.reaction.type.empty())
        throw NotPorted("reaction '" + ReactionTitle(setup_.flags, params_) +
                        "' is not ported yet (only H(e,e'p) and doing_phsp; P8+)");
    const bool quiet = opt_.quiet_stages;

    // ---- the stages, created by name and wired through the StageContext (interfaces only)
    StageContext ctx(setup_);
    if (!rngPtr_) {
        // Mt64Random lives in the STATIC library simcfw_random: constructing it here also makes the linker keep its
        // object file, hence its SIMC_REGISTER(IRandom, "mt64") registration.
        if (pipe_.rng.type == "mt64" && !Registry<IRandom>::contains("mt64"))
            rng_ = std::make_unique<Mt64Random>();
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

    reaction_ = Registry<IReaction>::create(pipe_.reaction.type);
    reaction_->configure(withVerbose(pipe_.reaction.config, quiet), params_);
    reaction_->init(ctx);

    model_ = Registry<IPhysicsModel>::create(pipe_.physics.type);
    model_->configure(withVerbose(pipe_.physics.config, quiet), params_);
    model_->init(ctx);

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

    ctx_ = reactions::MainContext::FromSetup(setup_);

    // ---- simc.f:85 calculate_central (after radc_init; central%rad from the radiator, simc.f:1268-1281)
    central_ = calculateCentral(setup_, *reaction_, *radiator_, *model_, ctx_);

    // ---- report bookkeeping
    report_ = RunReport{};
    report_.run_name = opt_.run_name;
    report_.ngen = params_.ngen;
    report_.random_seed = params_.random_seed;
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

    // l.195-200 ... keep the human interested
    if (T.ntried % ninform_ == 1) {
        fortran::Record r;
        r.x(1).lit("Generating Event").i(T.ntried, 9).lit(" ... ").i(T.nevent, 8).lit("  successes so far - Monitor:");
        r.e(T.wtcontribute * S.luminosity / static_cast<double>(T.ntried), 11, 4);
        log(r.str());
        if (T.ntried >= 5000) ninform_ = 20000;
    }

    model_->beginEvent(T.nevent);
    // l.203 generate
    bool success = reactions::generateEvent(*beam_, *reaction_, main, ev.vertex, ev.orig, rng);
    // l.211 montecarlo
    if (success) success = mc_->run(ev.orig, main, ev.recon, rng);
    // l.216-218 complete_recon_ev
    if (success) success = reaction_->completeRecon(ev.recon);
    // l.224 complete_main
    if (success) success = reactions::completeMain(false, main, ev.vertex, central_.vertex0, ev.recon, *model_, ctx_);

    // l.226-238 pass_cuts. FORTRAN-QUIRK(D-006, simc.f:234): the hadron delta upper bound uses SPedge%e%delta%max.
    {
        const auto& E = S.SPedge;
        const auto& M = S.slop.MC;
        const Event& r = ev.recon;
        ev.pass_cuts = !(r.e.delta <= (E.e.delta.min + M.e.delta.used) || r.e.delta >= (E.e.delta.max - M.e.delta.used) ||
                         r.e.yptar <= (E.e.yptar.min + M.e.yptar.used) || r.e.yptar >= (E.e.yptar.max - M.e.yptar.used) ||
                         r.e.xptar <= (E.e.xptar.min + M.e.xptar.used) || r.e.xptar >= (E.e.xptar.max - M.e.xptar.used) ||
                         r.p.delta <= (E.p.delta.min + M.p.delta.used) ||
                         r.p.delta >= (E.e.delta.max - M.p.delta.used) ||   // FORTRAN-QUIRK(D-006)
                         r.p.yptar <= (E.p.yptar.min + M.p.yptar.used) || r.p.yptar >= (E.p.yptar.max - M.p.yptar.used) ||
                         r.p.xptar <= (E.p.xptar.min + M.p.xptar.used) || r.p.xptar >= (E.p.xptar.max - M.p.xptar.used));
    }
    // l.240-243
    if (params_.hard_cuts) {
        if (!ev.pass_cuts) success = false;
        if (S.flags.doing_eep && (ev.recon.Em > S.cuts.Em.max)) success = false;
    }
    ev.success = success;
    // l.245
    if (success) T.sum_sigcc = T.sum_sigcc + main.sigcc;

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
            const double dp = r.p.delta - main.SP.p.delta, xp = r.p.xptar - v.p.xptar, yp = r.p.yptar - v.p.yptar,
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
                                T.slop);
    }

    // l.341-347
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
    report_.stops = collectStops();
    finished_ = true;

    // l.446: if (ngen.eq.0) goto 1000 — no output files
    if (params_.ngen == 0) {
        log(" ngen = 0: no output (simc.f:446)");
        return;
    }
    const RunSummary summary = MakeRunSummary(setup_, report_);
    for (auto& s : sinks_) {
        s->endRun(summary);
        if (auto* rs = dynamic_cast<IReportSink*>(s.get())) rs->writeReport(setup_, report_);
    }
    log(" Come back soon!");
}

void SimcRun::run() {
    initialise();
    while (!done()) step();
    finish();
}

}  // namespace driver
}  // namespace simc
