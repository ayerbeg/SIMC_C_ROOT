# 11. Performance tips

Measured on the development machine (Intel i7-3770, 4 cores / 8 hardware threads, GCC 16.2 `-O3`; `docs/PERFORMANCE.md`).

## 11.1 Speed

One thread, whole program including the start-up of about 0.5 s:

| case | µs per try | SIMC++ (1 thread) vs Fortran |
|---|---|---|
| H(e,e'p) radiative | 16.3 | 1.7× faster |
| π+ H, decay on | 27.1 | 3.4× faster |
| SIDIS | 7.2 | 3.0× faster |
| rho | 14.2 | 2.8× faster |
| π0 calorimeter (NPS) | 5.8 | 3.2× faster |
| polarized-target field | 40.7 | 1.5× faster |

(The Fortran times come from the reference runs and are indicative only.) The cost of a run is (number of tries) × (cost per
try). The number of tries per accepted event is Ntried/Ncontribute in the `.hist`: about 21 for the H(e,e'p) example,
about 3000 for rho.

## 11.2 Tips

* **Use threads.** `--threads 4` gave 2.8-3.4× and `--threads 8` 2.8-4.3× on the 4-core machine (`docs/PERFORMANCE.md`
  §3). Beyond the number of physical cores only hyper-threading helps. Threads are not allowed with `--fortran-compat` or
  RANLUX ([8.5](08_running_and_reproducibility.md#85-threads---threads-n)).
* **Memory with threads.** Each thread holds its own tables; a pion run on 8 threads uses about 450 MB. Each thread also
  builds its stages at the start (0.1-0.3 s per thread), which dominates very short runs.
* **Set tight `SPedge` limits.** The generation volume comes from them, plus the slop. Wide-open defaults (±1e10) waste
  almost every try, and the run may never end ([12](12_faq.md#the-run-never-ends)).
* **Test with negative `ngen`.** `--ngen -100000` makes exactly 100000 tries, whatever the acceptance.
* **mt64 is faster than RANLUX** (by about 25 %). Use RANLUX only to reproduce a Fortran run.
* **Leave out outputs you do not need.** For example, a `pipeline.sinks` list without `gen_hist` (or `Nntu = 0` for no tree).
  The `.hist` is written by `fortran_hist`.
* **Do not build with `-march=native` or `-ffast-math`.** They change the floating-point results; the bitwise tests and
  the exact reproduction of Fortran runs then fail.
* Radiative runs spend most of their time in the radiative-correction set-up of each try (the gamma function), the
  target energy loss and the COSY maps. Runs without radiation are faster.
