# Comparing SIMC++ with the Fortran SIMC (user guide)

Two kinds of comparison are possible:

* **exact**: RANLUX + `--fortran-compat` against a Fortran run of the same input. The files must be identical
  ([8.4](../../manual/markdown/user/08_running_and_reproducibility.md#84-reproducing-legacy-simc-results-exactly));
* **statistical**: any SIMC++ run (default mode, mt64, any seed, any number of threads) against a Fortran run, with
  `simc_compare`. This guide covers it. The full description of the statistics is in `validation/README.md`.

## 1. Building `simc_compare`

`simc_compare` is a separate CMake project in `validation/` (user manual [2.3](../../manual/markdown/user/02_installation.md#23-building)):

```bash
cmake -S $SIMCPP_SRC/validation -B $SIMCPP_BUILD-validation -DCMAKE_PREFIX_PATH=$ROOTSYS -DCMAKE_BUILD_TYPE=Release
cmake --build $SIMCPP_BUILD-validation -j8
```

## 2. Running it

```
simc_compare --a SPEC --b SPEC --out DIR [--config CFG.json] [--label-a NAME] [--label-b NAME]
             [--skip v1,v2] [--nbins N] [--prefer-hist] [--no-pdf] [--expect PASS|FAIL|IDENTICAL]
             [--mode-a compat|corrected] [--mode-b compat|corrected]
```

A `SPEC` is a SIMC++ `<run>.summary.json`, a Fortran `<run>.hist` (the tool reads `.hist`, `.geni` and `.root` next to
it), a `<run>.root`, or a directory with one run. Example: a full SIMC++ run of `eep_h_elastic_rad_A` (default mode,
seed 4242, 4 threads) against the Fortran reference run (seed 12345):

```bash
simc++ $SIMC_FORTRAN/infiles/eep_h_elastic_rad_A.inp -o l2 --seed 4242 --threads 4 > l2.log
$SIMCPP_BUILD-validation/simc_compare --a $SIMC_REFERENCE/runs/eep_h_elastic_rad_A --b l2/eep_h_elastic_rad_A.summary.json --out cmp \
                       --config $SIMCPP_SRC/validation/config/compare_default.json
```

It writes `cmp/summary.txt`, `cmp/summary.json`, `cmp/compare.pdf` (one page per variable: weighted and unweighted
overlays with B/A ratios) and the two parsed run summaries. The exit status is 0 for PASS (or when `--expect` is met),
1 for FAIL, 2 for a usage or I/O error.

## 3. Reading the verdict

The top of `summary.txt` from the example:

```
VERDICT: PASS   tests counted N=214  failed (p<0.01) F=0  expected 2.14  fraction 0.0000  P(X>=F|N,alpha)=1.000e+00 (need >= 1.000e-03)  gross failures=0 (det mismatch or p<1.000e-05)
worst p = 1.485e-02  (ssdeltai [chi2_ww])   largest |pull| = 1.859  (hms.scin/Ntried)

RUN-LEVEL TESTS
quantity                       method          A       errA             B       errB     pull          p  status
acceptance Ncontribute/Ntried  binomial  0.0481844   0.000105     0.0481577   0.000105    0.179      0.858  ok
yield sum(Weight) [Weight/Ntried] pull  2.48294e-05   5.64e-08   2.48231e-05   5.64e-08   0.0798      0.936  ok
wtcontr                        pull        207.354      0.472       207.306      0.472   0.0726      0.942  ok
...
TREE VARIABLES (p-values; chi2_abs: yield/tried incl. normalisation, chi2_ww/uu: TH1::Chi2Test shape, ks_w: unbinned weighted KS)
variable   status              lo          hi    chi2_abs     chi2_ww     chi2_uu        ks_w  result
hsdelta    compared        -7.005        7.05        0.17       0.147        0.15       0.266  ok
...
```

* **Run-level tests**: the acceptance Ncontribute/Ntried (binomial), the weighted yield per try, `wtcontr`,
  `AVERAGE.sigcc` (pulls), the deterministic quantities (`sigcc_central`, charge, `targetfac`, luminosity, `genvol`,
  `Ngen`, which must agree to the print precision; a mismatch is a **gross failure**), and every stop counter
  (binomial per counter, plus a χ² of each spectrometer's partition of its trials).
* **Tree variables**: for each branch, four tests on weighted histograms and the unbinned weighted KS test. Events with
  NaN/Inf values are counted and excluded.
* **Verdict**: with N tests at α = 0.01 about N·α failures are expected by chance. **PASS** if the number of failures F is
  compatible with that (binomial tail probability ≥ 0.001) **and** there is no gross failure (no deterministic
  mismatch, no p < 1e-5). Otherwise FAIL. A handful of failures with p around 0.001-0.01 is normal; a single gross
  failure is not.

Configuration (`validation/config/compare_default.json`, all keys optional): `alpha`, `global_alpha`, `gross_p`,
`nbins`, per-variable binning, `variables`, `skip`, the normalisation of each side (`norm_a`/`norm_b`: per try,
normfac or a number), `det_rel_tol`.

Comparing two runs of the **same** code with the same seed must give `IDENTICAL` (`--expect IDENTICAL`).

## 4. Mode of each side: `wtcontr`

The σ of `wtcontr` uses the events of the pass_cuts window of each run. That window differs between the Fortran and the
corrected mode (D-006). The tool takes the mode of a SIMC++ JSON from its `run.fortran_compat`; a Fortran `.hist` counts as
compat. For a summary without that key, pass `--mode-a/--mode-b compat|corrected`. A warning appears when the window does
not reproduce the run's own wtcontr within 0.2 % (usually a wrong mode).

## 5. Expected differences in the default mode

The default mode corrects clear Fortran bugs, so some comparisons against Fortran runs **fail on purpose**. The project's
L2 runs (`docs/translation/markdown/FORTRAN_DIFFERENCES.md` §7-§17) found these typical cases:

| run type | what fails | why |
|---|---|---|
| any run where recon hadron delta reaches the electron delta limit (e.g. HMS hadron with a different delta window; spect_mode runs) | `wtcontr` (pull up to +27 σ in spect_mode runs), resolutions | D-006 pass_cuts: the Fortran uses the electron delta max for the hadron |
| pion, kaon, Delta | `sigcc_central` (gross: deterministic mismatch) | D-038 F5 / D-043 F13 / D-048 F9: the central cross section is computed at the real central kinematics |
| kaon with decay OFF | weights, `decdist` | D-043 F5: survival probability with the kaon momentum |
| kaon | `saghai` column | D-043 F1-F3 (0 in the Fortran) |
| SIDIS, rho | `pfermi` (Fortran NaN), `wtcontr` (D-006), `decdist` (survival momentum) | D-051 F3, F6 |
| Delta | `sigcm`, Ntried and wtcontr with radiation | D-048 F2, F3 |
| HRS | acceptance, FP path, decay | D-061 F1, F1e |
| rho | random stream (no event-by-event comparison), energy loss of the ρ | D-057 F1, F2 |

For each failure, check the variable against [the list of differences](user_differences.md). If it is not listed there,
look further: run the corrected and the compat SIMC++ runs with the same seed, which differ only by the corrections, or
repeat the Fortran-mode comparison with RANLUX. Since D-058, the project's first check of a doubtful Fortran-mode flag
is a RANLUX run: if it is identical to the Fortran, the flag was a fluctuation.
