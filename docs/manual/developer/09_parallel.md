# 9. Parallelism and performance internals

User-level description, numbers and the L2 validation of the threaded runs: [PERFORMANCE.md](../../PERFORMANCE.md).
This chapter is the implementation view.

## 9.1 The parallel event loop (`SimcRun::runThreaded`, src/driver/SimcRunParallel.cpp)

Not a Fortran feature. `SimcRun::run()` calls `runThreaded(n)` when `threadCount() > 1` (`RunOptions::threads`, i.e.
`--threads N`, else the pipeline `"threads"`, else 1). Sequence:

1. **Master** `initialise()`: the set-up, the messages, the master stages and sinks (`beginRun`). The master's stages
   are used only for the central event and the report; the master runs no tries.
2. **Refusals** (`threadsRefusal(n)`, `std::invalid_argument`, `simc++` exit 2): `fortran_compat` (the exact
   reproduction needs one stream in the Fortran order), an rng other than `mt64` or a source set with `setRandom`
   (RANLUX and scripted streams are sequential), a sink that does not implement `IParallelSink`.
3. **Workers**, built one after the other in the master thread: `nw = min(n, |ngen|)` `SimcRun` objects with the same
   parameters and pipeline but no registry sinks and a silent driver log; worker k gets `ngen_k = ngen/nw` (+1 for
   `k < ngen % nw`, sign kept: successes or tries), its own `Mt64Random` seeded with `seedStream(base, k)`, and one
   worker sink per master sink (`makeWorkerSink(k, nw)`, added with `addSink`; `nullptr` = none). Each worker's
   `initialise()` builds its **own** stages (target, radiator, beam, model, reaction, spectrometers) and set-up.
4. **Threads**: `std::thread` per worker running `while (!done()) step()`; nothing is exchanged during the loop except
   two relaxed atomics per worker (tries, events) for the progress line the master prints every 10 s. An exception in a
   worker sets an abort flag; after the join the first failing worker's exception (in worker order) is rethrown.
5. **Merge, in the fixed worker order 0..nw−1**: `MergeTotals` (driver/RunReport.h: ntried, nevent, ncontribute,
   npasscuts, wtcontribute, sum_sigcc and the resolution sums added; the contrib/slop ranges combined as min/max with the
   `update_range` comparisons, `reactions::mergeContrib/mergeSlops`), the stop counters (summed per key), every sink
   (`mergeWorkerSink`). The end-of-run target state (`finalTargetState`) is the last worker's.
6. **Master** `finish()`: `normfac = luminosity * genvol / ntried` from the merged totals, then every output exactly as
   in a one-thread run; the summary gets `run.threads = nw`.

**Streams**: `Mt64Random::seedStream(base, k)` seeds `std::mt19937_64` with `std::seed_seq{base & 0xffffffff, base >>
32, k, 0x53494d43}`; `base = Mt64Random::BaseSeed(config, params)` is the seed the one-thread run would use. A stream
depends only on (base, k), not on N, and differs from the one-thread stream (`seed(base)`).

**Determinism**: no result depends on thread timing (fixed split, fixed streams, fixed merge order), so (seed, N) always
gives the same files (`threads_reproducibility_eep/_poltar`); `--threads 1` is the sequential loop bit for bit
(`threads_one_is_sequential`); refusals: `threads_refusal`. The threaded output is statistically equivalent to, not
bitwise equal to, the one-thread output (validated by L2, PERFORMANCE.md §3).

**Sinks** (`IParallelSink` implementations): `RootTreeSink` (worker temporary files `<run>.root.w<k>.tmp`, entries
appended in worker order, files removed), `GenHistSink` and `FortranHistSink` (histograms summed; TH1D errors in
quadrature), `RunSummarySink` (`nullptr`: works at `endRun` from the merged report).

**Thread safety**: one `SimcRun` = one Fortran process; every stage and every COMMON/SAVEd mirror is per worker. The only
process-wide mutable state is the Registry (read-only after static initialisation; workers are built in the master
thread) and the `NtupleLayout` column lists (guarded by a mutex). ROOT: per-worker `TFile`/`TTree`/`TH1D` without
directory, `ROOT::EnableThreadSafety()` called in `makeWorkerSink`. Stage warnings (Fortran `write(6,*)` mirrors) go to
`std::cout` from each worker and may interleave (`--quiet` silences them). ThreadSanitizer (clang, `-fsanitize=thread`)
runs of six representative inputs report no race (P15.4). A new stage with a `static` mutable cache breaks this
(chapter 8 §8.7).

**Cost**: each worker builds its stages (0.1-0.3 s, mostly the pion MAID tables) and holds its own tables (pion: about
+15 MB per worker).

## 9.2 Single-thread hot spots

Profile (callgrind, PERFORMANCE.md §1): before P15.4 libm `tgamma` was 20-56 % of radiative runs (29 calls per try in
init.f:radc_init_ev/basicrad_init_ev and radc.f:peaked_rad_weight, on 12 distinct arguments), then the COSY `transp`
maps, `trip_thru_target` (energy loss: log/pow), `generate_rad`, the target-field Runge-Kutta (`TargetField::rk4`), the
MAID table parsing at start-up, and `std::map` lookups in `GenHistSink::fill`.

Bitwise-neutral optimisations (P15.4), each proved by the golden set (`p153_golden.sh`, 326 identical files) and the
replay/RANLUX tests:

| optimisation | where | why it cannot change a bit |
|---|---|---|
| `tgamma` memo (16 slots, key = argument bit pattern) | `radiation::TgammaMemo` in `RadcState` | pure function: the cached value is the value the call returns |
| COSY power tables | `CosyEngine::transp`/recon (`PowerTable`) | the same `powi(x, k)` computed once per call and variable; products formed in the same order |
| MAID tables with `std::from_chars` | `pion::MaidTables::load` | correctly rounded like `strtod`; anything else falls back to `strtod` |
| `GenHistSink` histograms resolved once | `GenHistSink::inc(slot, ...)` | same TH1D calls without map lookups |

Rules for future optimisations: an optimisation in a kernel must leave every bit unchanged (verify with
`p153_golden.sh check` and the `replay|ranlux` tests); no fast-math, no `-march`, no reordering of floating-point
expressions, no replacement of libm calls by approximations. Remaining ideas are listed in PERFORMANCE.md §4 (a faster
correctly rounded gamma, parallel worker construction).
