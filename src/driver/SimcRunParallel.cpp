// src/driver/SimcRunParallel.cpp — the parallel event loop of SimcRun (docs/manual/markdown/reference/PERFORMANCE.md).
//
// Not a Fortran feature: the Fortran SIMC (and SimcRun::step) runs one try after the other on one random stream. Here
// N worker SimcRun objects, each with its own stages built from the same input (set-up, target, radiator, beam, physics,
// reaction, spectrometers) and its own mt64 stream Mt64Random::seedStream(seed, k), run their share of the events at
// the same time; the run's (master) object keeps the set-up, the log and the output sinks:
//   1. master.initialise(): the set-up, the messages, the sinks (beginRun)
//   2. worker k = 0..N-1 (built one after the other in this thread): pipeline without sinks, ngen_k = ngen/N (+1 for
//      k < ngen%N; the sign of ngen kept: successes or tries), log silenced, the worker sinks of IParallelSink
//   3. N threads run `while (!worker.done()) worker.step()`
//   4. merge in the fixed worker order 0..N-1: RunTotals (MergeTotals), stop counters (summed per key), sinks
//      (IParallelSink::mergeWorkerSink); the end-of-run target state (IReaction::finalTargetState) is the last worker's
//   5. master.finish(): normalisation and the outputs, exactly as for one thread
// Nothing depends on the thread timing, so a given (seed, N) always gives the same outputs. The workers share no
// mutable state: every stage is per worker; the only process-wide caches (NtupleLayout column lists, the Registry) are
// filled before the threads start (and guarded). ROOT: every worker owns its TFile/TTree/TH1 (the ROOT sinks call
// ROOT::EnableThreadSafety() in makeWorkerSink).
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <exception>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <streambuf>
#include <thread>

#include "simc/core/random/Mt64Random.h"
#include "simc/driver/FortranWriter.h"
#include "simc/driver/SimcRun.h"
#include "simc/stages/IParallelSink.h"

namespace simc {
namespace driver {

namespace {

/// A std::ostream that discards everything (the workers' driver log: the master prints the run's messages once).
class NullBuffer : public std::streambuf {
protected:
    int overflow(int c) override { return c; }
};
struct NullStream : std::ostream {
    NullStream() : std::ostream(&buf) {}
    NullBuffer buf;
};

}  // namespace

int SimcRun::threadCount() const {
    if (opt_.threads > 0) return opt_.threads;
    if (pipe_.threads > 0) return pipe_.threads;
    return 1;
}

std::string SimcRun::threadsRefusal(int n) const {
    if (n <= 1) return "";
    const std::string what = "--threads " + std::to_string(n) + " (more than one event-loop thread) ";
    if (setup_.fortran_compat)
        return what + "cannot be used with --fortran-compat: the exact Fortran reproduction needs the single random "
                      "stream in the Fortran order (run with --threads 1)";
    if (rng_ == nullptr || dynamic_cast<Mt64Random*>(rng_.get()) == nullptr || pipe_.rng.type != "mt64")
        return what + "needs the rng \"mt64\" (the pipeline rng is \"" + pipe_.rng.type +
               "\"): the RANLUX and scripted sources replay one sequential stream (run with --threads 1)";
    for (const auto& s : sinks_)
        if (dynamic_cast<IParallelSink*>(s.get()) == nullptr)
            return what + "needs output sinks that can merge worker outputs (IParallelSink); sink '" + s->name() +
                   "' cannot (run with --threads 1)";
    return "";
}

void SimcRun::runThreaded(int n) {
    initialise();
    if (const std::string why = threadsRefusal(n); !why.empty()) throw std::invalid_argument(why);
    if (params_.ngen == 0) {   // simc.f:446: no events, no outputs
        finish();
        return;
    }
    const long total = std::labs(static_cast<long>(params_.ngen));
    const int nw = static_cast<int>(std::min<long>(n, total));   // (no idle worker when ngen < n)
    const long sign = params_.ngen < 0 ? -1 : 1;
    const std::uint64_t base = Mt64Random::BaseSeed(pipe_.rng.config, params_);
    log(" SIMC++: event loop on " + std::to_string(nw) + " threads (mt64 streams seed_seq(" + std::to_string(base) +
        ", worker); docs/manual/markdown/reference/PERFORMANCE.md)");

    // ---- 2. the workers, built one after the other in this thread
    NullStream nullLog;
    Json workerPipeline = pipelineJson_.is_object() ? pipelineJson_ : Json::object();
    workerPipeline["sinks"] = Json::array();
    workerPipeline.erase("threads");
    std::vector<std::unique_ptr<SimcRun>> workers;
    std::vector<std::vector<IOutputSink*>> workerSinks(static_cast<std::size_t>(nw));   // [k][master sink index]
    for (int k = 0; k < nw; ++k) {
        RunOptions wo = opt_;
        wo.threads = 1;
        wo.override_ngen = true;
        wo.ngen = sign * (total / nw + (k < total % nw ? 1 : 0));
        wo.override_seed = false;   // (params_ already holds the run's seed)
        wo.log = &nullLog;
        auto w = std::make_unique<SimcRun>(params_, workerPipeline, wo);
        auto rng = std::make_unique<Mt64Random>();
        rng->seedStream(base, static_cast<std::uint64_t>(k));
        w->setRandom(std::move(rng));
        for (const auto& s : sinks_) {
            std::unique_ptr<IOutputSink> ws = dynamic_cast<IParallelSink&>(*s).makeWorkerSink(k, nw);
            workerSinks[static_cast<std::size_t>(k)].push_back(ws ? &w->addSink(std::move(ws)) : nullptr);
        }
        w->initialise();   // (stage warnings: as the master's, std::cout unless --quiet)
        workers.push_back(std::move(w));
    }

    // ---- 3. the threads
    std::vector<std::atomic<long>> tried(static_cast<std::size_t>(nw)), made(static_cast<std::size_t>(nw));
    std::vector<std::exception_ptr> errors(static_cast<std::size_t>(nw));
    std::atomic<bool> abort{false};
    std::atomic<int> running{nw};
    std::mutex m;
    std::condition_variable cv;
    std::vector<std::thread> threads;
    for (int k = 0; k < nw; ++k) {
        threads.emplace_back([&, k] {
            SimcRun& w = *workers[static_cast<std::size_t>(k)];
            try {
                while (!w.done() && !abort.load(std::memory_order_relaxed)) {
                    w.step();
                    tried[static_cast<std::size_t>(k)].store(w.totals().ntried, std::memory_order_relaxed);
                    made[static_cast<std::size_t>(k)].store(w.totals().nevent, std::memory_order_relaxed);
                }
            } catch (...) {
                errors[static_cast<std::size_t>(k)] = std::current_exception();
                abort = true;
            }
            {
                std::lock_guard<std::mutex> lk(m);
                --running;
            }
            cv.notify_all();
        });
    }
    // progress lines (the Fortran "Generating Event" monitor, summed over the workers) every 10 s
    {
        std::unique_lock<std::mutex> lk(m);
        while (!cv.wait_for(lk, std::chrono::seconds(10), [&] { return running.load() == 0; })) {
            long nt = 0, ne = 0;
            for (int k = 0; k < nw; ++k) {
                nt += tried[static_cast<std::size_t>(k)].load(std::memory_order_relaxed);
                ne += made[static_cast<std::size_t>(k)].load(std::memory_order_relaxed);
            }
            fortran::Record r;
            r.x(1).lit("Generating Event").i(nt, 9).lit(" ... ").i(ne, 8).lit("  successes so far (" +
                                                                                std::to_string(nw) + " threads)");
            log(r.str());
        }
    }
    for (auto& t : threads) t.join();
    for (const auto& e : errors)
        if (e) std::rethrow_exception(e);   // (the first failing worker in the worker order)

    // ---- 4. merge, worker order 0..nw-1
    RunTotals& T = report_.totals;
    for (const auto& w : workers) MergeTotals(T, w->totals());
    mergedStops_ = workers.front()->collectStops();
    for (std::size_t k = 1; k < workers.size(); ++k) {
        const std::vector<StopCounters> s = workers[k]->collectStops();
        for (std::size_t a = 0; a < s.size() && a < mergedStops_.size(); ++a)
            for (const auto& [key, v] : s[a].counts) mergedStops_[a].counts[key] += v;
    }
    mergedStopsSet_ = true;
    for (std::size_t i = 0; i < sinks_.size(); ++i) {
        auto& ps = dynamic_cast<IParallelSink&>(*sinks_[i]);
        for (std::size_t k = 0; k < workers.size(); ++k)
            if (IOutputSink* ws = workerSinks[k][i]) ps.mergeWorkerSink(*ws);
    }
    finalStateReaction_ = workers.back()->reaction_.get();
    report_.threads = nw;

    // ---- 5. the outputs (the workers stay alive until here: finalStateReaction_ points into the last one)
    finish();
    finalStateReaction_ = nullptr;
}

}  // namespace driver
}  // namespace simc
