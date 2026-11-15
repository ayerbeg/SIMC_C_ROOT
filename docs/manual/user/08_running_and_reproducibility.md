# 8. Running modes and reproducibility

## 8.1 Seeds

The seed comes from, in this order: the pipeline `"rng": {"seed": N}`, then `--seed N`, then the input's `random_seed`.
`--seed` overrides `random_seed`, and the `.hist` prints the overridden value. With the default generator:

* `random_seed > 0`: `std::mt19937_64` seeded with it;
* `random_seed = 0` or absent (and other values ≤ 0 except −1): the fixed seed 5489, so the run is still reproducible;
* `random_seed = −1`: the current time, as in the Fortran (dbase.f:119). The log line ` Use random seed = ...` shows the value.

A given input, seed, mode and number of threads always produce the same files.

## 8.2 Random-number generators: `mt64` and `ranlux`

| generator | pipeline | properties |
|---|---|---|
| `mt64` (default) | `"rng": "mt64"` or `{"type": "mt64", "seed": N}` | std::mt19937_64. Fast; the project's statistical validation (L2) uses it (decision D-002) |
| `ranlux` | `"rng": {"type": "ranlux"}`, optionally `"seed": N` or `"state": [25 integers]` | a bit-exact port of the Fortran RANLUX (`call_ranlux.f`, luxury level 3, seeded with RLUXGO(3, seed, 0, 0) from `random_seed`, as the Fortran). About 25 % slower (D-058). `"state"` takes the 25 integers of a Fortran `*_start_random_state.dat` |

`ranlux` gives the same random numbers as the Fortran. It reproduces a Fortran run only in the Fortran-compatible mode,
because several corrections of the default mode change the number of random numbers drawn (for example D-057 F2 for
rho runs).

## 8.3 Corrected mode and `--fortran-compat`

| | default (corrected) | `--fortran-compat` |
|---|---|---|
| clear Fortran bugs ([chapter 10](10_differences_from_fortran.md)) | corrected | kept |
| ambiguous cases (expert list) | Fortran | Fortran |
| reproduces a Fortran run exactly (with `ranlux`) | no | yes |
| `--threads N > 1` | yes | refused |
| use it for | physics | validation, and comparing with an existing Fortran result |

In a JSON input: `"run": {"fortran_compat": true}`. The log says which mode is active.

## 8.4 Reproducing a Fortran run exactly

Run the same `.inp` with RANLUX and `--fortran-compat`, in one thread. Example with the reference run of the converted
`test_eep_h` input (1000 events, with an `extra_dbase_file`):

```bash
echo '{"pipeline": {"rng": {"type": "ranlux"}}}' > ranlux.json
R=$SIMC_REF/runs/eep_h_test_pctnames            # the Fortran run: .inp, extra file, .hist, .geni, .gen, .root
simc++ $R/eep_h_test_pctnames.inp -o exact --fortran-compat --pipeline ranlux.json > exact.log
diff <(sed 2,3d exact/eep_h_test_pctnames.hist) <(sed 2,3d $R/eep_h_test_pctnames.hist) && echo ".hist: identical apart from the time lines"
cmp exact/eep_h_test_pctnames.geni $R/eep_h_test_pctnames.geni && cmp exact/eep_h_test_pctnames.gen $R/eep_h_test_pctnames.gen && echo ".geni, .gen: identical"
simc_root_identity exact/eep_h_test_pctnames.root $R/eep_h_test_pctnames.root | tail -1
```

Output:

```
.hist: identical apart from the time lines
.geni, .gen: identical
exact/eep_h_test_pctnames.root vs .../reference/runs/eep_h_test_pctnames/eep_h_test_pctnames.root: IDENTICAL
```

Notes:

* The Fortran run must use the same input, including its `random_seed` (none = the RANLUX default seed), and must be
  made with the same toolchain class: the bitwise references depend on glibc's libm (decision D-044).
* Run the same command **without** `--fortran-compat` and the `.hist` differs in the corrected lines (here two lines
  where the Fortran writes integers as `T`/`F`, D-018).
* The project checked this on every runnable shipped input (39/39), the spect_mode and π0 probes and the
  polarized-target examples (CTest label `regression_shipped`). The script
  `validation/scripts/p15_shipped_regression.sh` runs the same procedure on a reference folder.

## 8.5 Threads: `--threads N`

`--threads N` (or `"pipeline": {"threads": N}`; the option wins) runs the event loop on N threads
(`docs/PERFORMANCE.md` §3):

* each thread has its own copy of every stage and its own mt64 random stream, derived from (seed, thread number);
* the events are split: thread k makes ngen/N events (+1 for k < ngen mod N), so the total is exactly ngen (successes
  for ngen > 0, tries for ngen < 0);
* the outputs are merged in a fixed thread order: the tree entries thread after thread, the histograms and counters
  summed, the normalisation computed from the merged totals;
* **a given seed and N always give the same files**, but they are statistically equivalent to, not equal to, the
  one-thread run, or a run with another N. `--threads 1` is exactly the sequential run;
* the summary JSON records `"run": {"threads": N}`;
* the stage warnings of the threads may be interleaved in the log (`--quiet` silences them).

```
$ simc++ run.json -o out_t4 --ngen 5000 --threads 4 > out_t4.log; grep -E "Ntried|Ncontribute" out_t4/run.hist
            Ntried         =     102516
            Ncontribute    =       5000
$ simc++ run.json -o out_t4b --ngen 5000 --threads 4 > /dev/null && simc_root_identity out_t4/run.root out_t4b/run.root | tail -1
out_t4/run.root vs out_t4b/run.root: IDENTICAL
```

**Why not with `--fortran-compat` or RANLUX:** the exact Fortran reproduction needs the single random stream, consumed in
the Fortran order, and N threads use N streams. Both combinations are refused before the run starts (exit status 2):

```
$ simc++ run.json -o out_compat --ngen 5000 --threads 4 --fortran-compat
simc++: --threads 4 (more than one event-loop thread) cannot be used with --fortran-compat: the exact Fortran reproduction needs the single random stream in the Fortran order (run with --threads 1)
```

The L2 comparison of 8-thread runs with the Fortran gives the same verdicts as the one-thread runs
(`docs/PERFORMANCE.md` §3).

## 8.6 Other run options

* `--ngen N` overrides `ngen` (use a small N for a test; a negative N gives |N| tries, which always ends).
* `--quiet` silences the stage messages (the Fortran-style progress lines stay).
* `--print-pipeline` prints the plug-ins chosen and continues with the run.
* `--run-name`, `-o`: output naming ([4.2.2](04_inputs.md#422-run)).
