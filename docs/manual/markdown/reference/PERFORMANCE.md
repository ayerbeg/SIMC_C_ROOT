# SIMC++ performance (task P15.4)

Two independent things make SIMC++ faster without touching the exact-Fortran reproduction:

1. **Single-thread optimisations that change no bit** of any output (section 2). Every `--fortran-compat` run, every
   RANLUX replay and every default (corrected, mt64) run with a given seed gives the same files as before.
2. **An opt-in parallel event loop**, `simc++ --threads N` (section 3). It uses N independent random streams, so its
   output is statistically equivalent to, not bitwise equal to, the one-thread run; a given (seed, N) is reproducible.

Machine of the numbers below: Intel i7-3770 (4 cores, 8 hardware threads, 3.4 GHz), 15 GB, Fedora 44 (glibc 2.43,
gcc 16.2.1, `-O3`, no `-march`/`-ffast-math`: both would change the floating-point results), ROOT 6.40.04.

## 1. Profile

callgrind (Release + `-g`) on short runs of six representative inputs (`reference/runs/<case>`, seed 4242, corrected
mode), before the optimisations (commit 0f7e4ef). Inclusive share of the instructions:

| case | hot spots (inclusive) |
|---|---|
| eep H rad (`eep_h_elastic_rad_A`) | libm `tgamma` 40% (all from `radc_init_ev`/`basicrad_init_ev`), COSY `transp` 11%, ROOT start-up (`TFile::Open`) 9%, `trip_thru_target` 8% |
| pion dec1 (`pion_h_k1_rad_dec1_A`) | `tgamma` 22%, MAID table parsing (`strtod`, start-up) 17%, COSY `transp` 25% (SHMS with decay), `generate_rad` 12% |
| SIDIS (`sidis_h_k1_pip_rad_A`) | `tgamma` 30%, COSY `transp` 24%, `trip_thru_target` 8%, `GenHistSink::fill` 5% (std::map lookups) |
| rho (`rho_h_k1_rad_A`) | `tgamma` 56%, `generate_rad` 33% (of which `peaked_rad_weight` 11%), energy loss 10% |
| calo pi0 (`p14f2_excl_h_ng2_A`) | `tgamma` 27%, MAID parsing (start-up) 25%, `trip_thru_target` 7% |
| poltar field (`shipped_poltar_exclusivepion_example`) | target-field Runge-Kutta (`TargetField::rk4`/`field`) 25%, `tgamma` 22%, MAID parsing 20% |

The radiative-correction set-up of every try (init.f:radc_init_ev → basicrad_init_ev, radc.f:peaked_rad_weight) calls
`gamma()` 29 times on only 12 distinct arguments (1+lambda(i), 1+bt(i), 1+g(i), ...), and glibc 2.43's correctly rounded
`tgamma` costs several thousand instructions per call: it was the first hot spot of every radiative run.

## 2. Single-thread optimisations (bitwise identical)

| change | where | why it cannot change a bit |
|---|---|---|
| `tgamma` memo (16 slots, key = the argument's bit pattern) | `radiation::TgammaMemo` in `RadcState` (kernel/RadcState.h); RadcInit.cpp, Radc.cpp | `tgamma` is a pure function: the cached value IS the value the call would return (12 instead of 29 calls per try) |
| COSY power tables | `optics::CosyEngine::transp`/`recon` (`PowerTable`) | `powi(x, k)` computed once per call and variable with the same `powi`, then looked up; the products are formed in the same order |
| MAID tables parsed with `std::from_chars` | `pion::MaidTables::load` | correctly rounded like `strtod` (same double); anything else falls back to `strtod` |
| `GenHistSink` histograms resolved once | `GenHistSink::inc(slot, ...)` | same TH1D calls, without two `std::map<std::string>` lookups per fill |

Proof: `validation/scripts/p153_golden.sh` (20 inputs covering every reaction family and arm type and the target field,
both modes, plus one `simc_compare`) gives 326 identical files against the golden set of commit c29875c, and every
RANLUX / replay CTest (bitwise against the Fortran) passes unchanged (73/73 of the `replay|ranlux|radc|transp...`
subset, then the full `check_all`).

After them (callgrind, same runs): eep H rad 4.18e9 → 3.04e9 instructions (tgamma 40 % → 23 %), rho 8.62e9 → 5.68e9
(tgamma 56 % → 34 %); the remaining hot spots are the 12 distinct tgamma calls per try, generate_rad, trip_thru_target
(energy loss, log/pow) and the COSY maps.

Single-thread wall time per try (and the Fortran's, from the reference `.hist` BEGIN/END times; the Fortran runs were
made on the same machine, possibly under other load):

| case (ngen, tries/event) | before (0f7e4ef) µs/try | after, 1 thread µs/try | gain | Fortran µs/try | SIMC++ (1 thread) vs Fortran |
|---|---|---|---|---|---|
| eep H rad (50000, 21) | 23.5 | 16.3 | 1.44x | 27.7 | 1.7x faster |
| pion dec1 (50000, 21) | 40.8 | 27.1 | 1.50x | 91.1 | 3.4x faster |
| SIDIS (20000, 148) | 10.5 | 7.2 | 1.45x | 21.8 | 3.0x faster |
| rho (1000, 3106) | 19.1 | 14.2 | 1.35x | 39.9 | 2.8x faster |
| calo pi0 (NPS) (10000, 287) | 8.1 | 5.8 | 1.38x | 18.6 | 3.2x faster |
| poltar field (5000, 225) | 49.1 | 40.7 | 1.21x | 59.4 | 1.5x faster |

(`validation/scripts/p154_timing.sh`, seed 4242, corrected mode, whole program including the start-up of about 0.5 s;
before = commit 0f7e4ef built the same way. The Fortran times come from runs made for the references, with larger ngen
and sometimes several runs at the same time on this machine, and include the Fortran's `.bin` ntuple writing: an
indication, not a benchmark.)

## 3. Parallel event loop: `--threads N`

```
simc++ run.inp --threads 8              # or "pipeline": {"threads": 8}; the command-line option wins
```

* **Design** (`SimcRun::runThreaded`, src/driver/SimcRunParallel.cpp). The run's (master) `SimcRun` builds the set-up,
  prints the messages and opens the output sinks as always. Then N worker `SimcRun`s are built one after the other (same
  parameters and pipeline, no registry sinks, silent driver log), each with its OWN stages (target, radiator, beam,
  physics model, reaction, spectrometers: no stage object and no Fortran-COMMON mirror is shared) and its own worker
  sinks. The N threads run `while (!done()) step()`; nothing is exchanged during the loop. After the join the master
  merges, in the fixed worker order 0..N-1: the run accumulators (`MergeTotals`: ntried, nevent, ncontribute,
  npasscuts, wtcontribute, sum_sigcc, the resolution sums; the contrib/slop ranges as min/max), the spectrometer stop
  counters (summed per key) and the sinks (`IParallelSink::mergeWorkerSink`). `finish()` then computes the
  normalisation (normfac = luminosity * genvol / ntried from the merged totals, genvol from the set-up) and writes every
  output exactly as for one thread. The `.hist` "Monitor"/progress line is printed every 10 s from the workers' counters.
* **Event split.** Worker k makes `ngen/N` (+1 for k < `ngen%N`) of the successes (`ngen > 0`) or of the tries
  (`ngen < 0`); the merged counts are exactly `ngen`. No more workers than |ngen| are started.
* **Random streams.** Worker k: `std::mt19937_64` seeded with `std::seed_seq{seed & 0xffffffff, seed >> 32, k,
  0x53494d43}` (`Mt64Random::seedStream`), seed = the run's mt64 seed (pipeline `rng.seed`, else `random_seed`). A
  stream depends only on (seed, k), not on N, and differs from the one-thread stream (`mt19937_64::seed(seed)`).
* **Reproducibility.** The merge order and every worker's stream are fixed, so a given (seed, N) always gives the same
  output files (CTest `threads_reproducibility_eep`, `threads_reproducibility_poltar`). `--threads 1` is the
  sequential loop, bit for bit (`threads_one_is_sequential`).
* **Outputs.** `<run>.root`: the workers fill their own trees in temporary files `<run>.root.w<k>.tmp`, appended to the
  run's tree in worker order and removed. `.gen.root`, `.gen`: the workers' histograms summed (TH1D: contents summed,
  errors in quadrature, entries summed). `.hist`, `.geni`, summary: from the merged totals; the summary JSON gets
  `"run": {"threads": N}` (absent for one thread). The end-of-run COMMON state that the report prints
  (`IReaction::finalTargetState`, targ%Mrec_struck of the pion Delta final states) is the last worker's.
* **Refused** (exit 2, with a message): `--fortran-compat` (the exact Fortran reproduction needs the one random stream in
  the Fortran order) and any rng other than `mt64` (RANLUX, scripted: one sequential stream); sinks that do not
  implement `IParallelSink` (CTest `threads_refusal`).
* **Thread safety.** Every stage, Fortran-COMMON mirror and SAVEd-local mirror is per worker (one `SimcRun` = one
  Fortran process). The only process-wide state is the Registry (read-only after static initialisation; the workers are
  built in the master thread) and the column lists of `NtupleLayout` (filled once, now under a mutex). ROOT: every
  worker owns its `TFile`/`TTree` and its `TH1D`s (no directory); the ROOT sinks call `ROOT::EnableThreadSafety()`
  before the threads start. The stage warnings (`write(6,*)` mirrors) go to `std::cout` from every worker, as in a
  one-thread run, so with `--threads` they may be interleaved (`--quiet` silences them).
* **Cost.** Each worker builds its stages (about 0.2 s, mostly the pion MAID tables) and holds its own tables (pion:
  +15 MB per worker); a pion run on 8 threads uses about 450 MB.

Speed-up (same inputs as section 2; wall time of the whole program including start-up; the machine has 4 physical
cores, so beyond 4 threads only hyper-threading helps):

| case | 1 thread | 2 threads | 4 threads | 8 threads | ms/event at 8 threads |
|---|---|---|---|---|---|
| eep H rad | 17.0 s (1.00x) | 9.6 s (1.77x) | 5.8 s (2.90x) | 5.5 s (3.07x) | 0.111 |
| pion dec1 | 28.7 s (1.00x) | 16.1 s (1.78x) | 9.5 s (3.02x) | 8.9 s (3.21x) | 0.179 |
| SIDIS | 21.4 s (1.00x) | 11.8 s (1.81x) | 6.9 s (3.09x) | 6.3 s (3.39x) | 0.315 |
| rho | 44.0 s (1.00x) | 22.1 s (1.99x) | 13.0 s (3.38x) | 11.1 s (3.97x) | 11.060 |
| calo pi0 (NPS) | 16.7 s (1.00x) | 9.5 s (1.76x) | 5.9 s (2.81x) | 6.0 s (2.76x) | 0.605 |
| poltar field | 45.7 s (1.00x) | 24.6 s (1.86x) | 14.9 s (3.07x) | 10.5 s (4.34x) | 2.108 |

The start-up (set-up, ROOT, and about 0.2 s per worker) is part of these short runs; for runs of minutes the speed-up
approaches the 4-thread / 8-thread asymptote (about 3.4x / 4x on this 4-core machine).

L2 validation (`validation/scripts/p154_l2.sh`): for each case a `--threads 8` and a `--threads 1` run with the
reference's ngen, `simc_compare` of t8 vs the Fortran reference A, t8 vs t1 and t1 vs A:

| case (ngen) | t8 vs Fortran A | t1 vs Fortran A | t8 vs t1 | wall t8 / t1 (s, loaded machine) |
|---|---|---|---|---|
| eep H rad (200000) | PASS 1/214 | PASS 0/214 | PASS 0/214 | 41 / 167 |
| pion dec1 (200000) | FAIL 4/271, 1 gross: `sigcc_central` (D-038 F5) | FAIL 3/273, the same gross | PASS 4/274 | 63 / 236 |
| SIDIS (100000) | FAIL 6/271, gross: `wtcontr` (D-006), `decdist`, `pfermi` (D-051 F3) | FAIL 10/271, the same gross | PASS 3/273 | 53 / 218 |
| rho (100000) | FAIL 5/280, 2 gross: `wtcontr` (D-006), `pfermi` (D-051 F3) | FAIL 2/280, the same gross (*) | PASS 3/282 (*) | 1105 / 6851 (*) |

`run-archive:<label>` names a run folder of the SIMC++ project's run archive (the `simc_runs/` directory next to the source tree, where `tools/check_all.sh` and the run script write `<timestamp>_<label>/` folders); `reference/...` is relative to the Fortran reference data `$SIMC_REFERENCE`.

(*) rho t1 = the accepted P13.2 corrected run (run-archive:P13.2_L2_corrected_h_k1_rad: same seed and ngen,
SIMC++ 75f8dbd before the optimisations; the one-thread outputs are unchanged since, as the pion/SIDIS t1 runs show).
F/N = statistical failures (p < 0.01) / tests. Every gross failure against the Fortran is a documented correction of
the default (corrected) mode, present identically in the one-thread run (FORTRAN_DIFFERENCES.md: §9 pion, §12 SIDIS,
§13 rho); `--threads` cannot run in the compat mode, so the threaded runs are compared in the corrected mode. The
one-thread runs (t1) of pion and SIDIS are bitwise identical (`.hist`, tree) to the accepted P9.2 / P12.2 corrected L2
runs of the same seed. Runs: run-archive:P15.4_L2_\<case\>_t{1,8}/ (compare_vs_A/, compare_vs_t1/).

ThreadSanitizer (clang 22, `-fsanitize=thread`): `--threads 4` runs of eep, pion dec1, SIDIS, rho, calo pi0 and the
polarized-target field input report no data race.

## 4. Not done (possible next steps)

* `tgamma` is still the largest single cost of radiative runs (12 calls per try). A faster correctly-rounded gamma
  would keep the bits but is a numerical library of its own; a non-correctly-rounded one would change results.
* The workers are built sequentially (0.1-0.3 s each); building them in parallel would need a thread-safety audit of
  every set-up reader.
* `-march=native` (FMA contraction) and `-ffast-math` are deliberately not used: they change the floating-point results
  and would break the bitwise Fortran comparisons.
