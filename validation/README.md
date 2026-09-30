# SIMC++ validation tools

## simc_compare (P1.5): are two SIMC runs statistically the same?

Standalone CMake project (it does not touch the top-level build):

```bash
export LD_LIBRARY_PATH=/home/ayerbe/Builds/ROOT6_install/lib:/usr/local/lib
cmake -S validation -B validation/build -DCMAKE_PREFIX_PATH=/home/ayerbe/Builds/ROOT6_install
cmake --build validation/build -j
(cd validation/build && ctest --output-on-failure)      # self-tests, about 40 s
validation/scripts/run_compare_selftests.sh             # same, outputs into reference/comparisons/
```

### Usage

```
simc_compare --a SPEC --b SPEC --out DIR [--config CFG.json]
             [--label-a NAME] [--label-b NAME] [--perturb VAR:SHIFT]...
             [--skip v1,v2] [--nbins N] [--prefer-hist] [--no-pdf]
             [--expect PASS|FAIL|IDENTICAL] [--expect-flagged v1,v2]
simc_compare --dump-summary SPEC OUT.json
```

`SPEC` is one of the following:

* `<run>.summary.json`: the run summary JSON (below). The tree is found through its `tree.file`,
  or else `<run>.root` next to it.
* `<run>.hist`: the Fortran text output. The tool reads `<run>.hist`, `<run>.geni` and `<run>.root`.
* `<run>.root`: the tool uses `<run>.summary.json` if present, otherwise `<run>.hist`/`.geni`.
* A directory that holds one run.

When both a summary JSON and a `.hist` exist, the JSON is used. `--prefer-hist` reverses this.
The `.hist`/`.geni` parser is the fallback until every producer writes the JSON.

Outputs in `--out`:

* `summary.txt`: verdict, run-level table and per-variable table.
* `summary.json`: the same content in machine-readable form, with every test.
* `compare.pdf`: summary pages, then one page per variable. The left panel is weighted and
  normalised to yield per tried event. The right panel is unweighted with unit area. Each panel
  has a B/A ratio below it.
* `runA.summary.json`, `runB.summary.json`: the run-level numbers as parsed.

Exit codes:

* 0: PASS, or the `--expect` expectation was met.
* 1: FAIL, or the expectation was not met.
* 2: usage or I/O error.

`--perturb VAR:SHIFT` is a test option. It adds `SHIFT`, in branch units, to `VAR` in run B.
For example, `hsdelta:0.1` is +0.1 % in HMS delta.

### Configuration (`config/compare_default.json`, all keys optional)

| key | default | meaning |
|---|---|---|
| `alpha` | 0.01 | a test fails if p < alpha |
| `global_alpha` | 0.001 | minimum binomial tail probability of the failure count |
| `gross_p` | 1e-5 | any counted test with p below this makes the verdict FAIL |
| `nbins` | 50 | bins per variable |
| `quantile_lo`/`quantile_hi` | 0.005/0.995 | auto range, taken from the combined (A+B) unweighted quantiles |
| `binning` | {} | per-variable `{"nbins":N,"min":x,"max":y}` |
| `variables` | [] | restrict to these branches (empty means every common branch) |
| `skip` | [] | branches not compared |
| `norm_a`/`norm_b` | `per_tried` | `per_tried` (W/Ntried), `normfac` (W*normfac/nevent), or a number (W*number) |
| `det_rel_tol` | 1e-6 | relative tolerance for deterministic scalars, added to the print precision |
| `tree_name`, `weight_branch`, `sigcc_branch` | h10, Weight, sigcc | |
| `tests` | all on | `chi2_abs`, `chi2_ww`, `chi2_uu`, `ks_weighted` |
| `stop_per_counter` | true | binomial pull for each stop counter |

### Statistics

In this section, N is Ntried, W_i is the `Weight` branch value, and the sums run over the tree
(successful events). Each pull is `(A-B)/sqrt(sA^2+sB^2)` and has a two-sided normal p-value.

* **Acceptance** `Ncontribute/Ntried`. This is a two-proportion z-test with the pooled
  proportion p = (kA+kB)/(NA+NB) and sigma^2 = p(1-p)(1/NA+1/NB). When ngen>0 the number of
  successes is fixed and Ntried is random (inverse sampling). The estimator k/N then has the
  same variance, a(1-a)/N, to leading order.
* **Nco_no_rad_prot/Ncontribute**: the same binomial test.
* **Stop counters** (`.geni`): the tool gives a binomial pull for every counter as a fraction of
  Ntried. It also runs a 2xK homogeneity chi^2 on each spectrometer's partition of its trials.
  The categories are all stop reasons plus `successes`, plus an `unaccounted` category equal to
  trials minus their sum. `trials`, `hut` and `q2_mid_dup` are not categories. ndf = populated
  categories - 1.
* **Weighted yield per tried event**. The event weight is w_i = W_i if the event is in the tree
  and 0 otherwise, for all N tried events. Then:
  `Y = S1/N`, `sigma_Y^2 = (S2 - S1^2/N)/N^2`, with `S1 = sum W`, `S2 = sum W^2`.
  This is the variance of the mean of the N per-try weights. It includes the fluctuation in the
  number of accepted events. With `norm = normfac`, the same result is multiplied by
  lumi*genvol (normfac*N/nevent). With a numeric scale s: Y = s*S1 and sigma = s*sqrt(S2).
* **wtcontr** (.hist) = sum_{pass_cuts} W * normfac/nevent (simc.f:313,398,920). Its sigma is
  wtcontr times the relative error of the tree yield. The tree holds all successes, while
  wtcontr uses only the `pass_cuts` subset. The tool warns if the two differ by more than 0.1 %.
* **AVERAGE.sigcc** = sum_{successes} sigcc / nevent (simc.f:248,959). For ngen>0, sigma is
  sd(sigcc)/sqrt(n). For ngen<0 (nevent = Ntried) it uses the per-tried formula above.
* **normfac**. For ngen>0 it is random through Ntried: sigma_rel = sqrt((1-acc)/Ncontribute).
  It is reported as info only and not counted, because it carries the same information as the
  acceptance test. For ngen<0 it is deterministic.
* **Deterministic scalars** (`sigcc_central`, `charge_mC`, `targetfac`, `luminosity_ub`,
  `luminosity_GeV2`, `genvol`, `Ngen`). These pass if
  `|A-B| <= det_rel_tol*max(|A|,|B|) + halfulp_A + halfulp_B`. halfulp is half a unit in the
  last printed digit: e16.6 or e16.8 in the .hist, and 0 for full-precision JSON. A mismatch is
  a gross failure.
* **Per tree variable**, for events with finite values and weights. Under- and overflow are
  included in all χ² tests.
  * `chi2_abs`: absolute comparison of the weighted histograms normalised to yield per tried
    event, `sum_bins (a-b)^2/(sa^2+sb^2)`, with `sa^2 = scale^2 * sum W^2` in the bin. ndf is the
    number of bins with non-zero error. This tests shape and normalisation together.
  * `chi2_ww`: `TH1::Chi2TestX(...,"WW UF OF")`, the Gagunashvili weighted/weighted
    homogeneity test (shape only, normalisation free).
  * `chi2_uu`: `TH1::Chi2TestX(...,"UU UF OF")` on unweighted counts (shape only).
  * `ks_w`: an **unbinned weighted two-sample KS test** (not TH1::KolmogorovTest).
    D = sup|F_A - F_B| over the weighted ECDFs. The effective sizes are
    n_eff = (sum w)^2/sum w^2, and p = TMath::KolmogorovProb(D*sqrt(nA nB/(nA+nB))).
  * A branch that is constant in both runs is not histogrammed. Its constant value is compared
    instead (a mismatch is a gross failure).
* **Verdict**. N is the number of counted tests, and F is the number with p < alpha.
  Expected = N*alpha. P_bin = P(X >= F | Binomial(N, alpha)) (TMath::BinomialI).
  **PASS** if P_bin >= global_alpha and there are no gross failures (no deterministic mismatch
  and no counted p < gross_p). Tests on the same quantity (for example `hsdelta`/`hsdeltai`, or
  the four tests of one variable) are correlated. The binomial count is therefore approximate
  and slightly conservative in the tail. In the Fortran A-vs-B runs the p-values are close to
  uniform (mean p = 0.45-0.57 for each method).

### .geni stop-counter mapping (by physical meaning)

Stop counters are stored under the Fortran COMMON variable name without its prefix
(`hSTOP_slit_hor` -> `slit_hor`). The mapping follows the variable actually written, not the
printed label. The following Fortran label mistakes are known:

* simc.f:477 (SOS) and simc.f:489/501 (HRS): the label says `Slit hor/vert` but the line writes
  `slit_vert, slit_hor`.
* simc.f:519-520 (SHMS): the `Q3 ... mid` slot prints `shmsSTOP_q2_mid`. This is stored as
  `q2_mid_dup`. `q3_mid` is never printed, so it falls into `unaccounted`.
* simc.f:526 (SHMS): the label `DC1, DC2, Scin, Cal` has only `dc1, dc2` written.

HMS keys: `trials slit_hor slit_vert slit_oct Q1_in Q1_mid Q1_out Q2_in Q2_mid Q2_out Q3_in Q3_mid
Q3_out D1_in D1_out hut dc1 dc2 scin cal successes`.
SHMS keys: `trials hb_in hb_men hb_mex hb_out slit_hor slit_vert slit_oct q1_in q1_men q1_mid q1_mex
q1_out q2_* q3_in q3_men q2_mid_dup q3_mex q3_out D1_in D1_flr D1_mid1..D1_mid7 D1_mex D1_out hut
dc1 dc2 s1 s2 s3 cal successes`.
SOS keys: `trials slit_* quad_in quad_mid quad_out bm01_in bm01_out bm02_in bm02_out exit hut dc1
dc2 scin cal successes`.
HRSr/HRSl keys: `trials slit_* Q1_* Q2_* D1_in D1_out Q3_* hut dc1 dc2 s1 s2 cal`.
Calo keys: `trials slit_hor slit_vert successes`.

## Run summary JSON (what the C++ writer should emit)

The authoritative schema is `docs/summary.schema.json` (owned by P1.4). This reader is tolerant:
it looks up each key by name, with a few aliases, up to three levels deep. The minimal content it
needs is shown below. Examples converted from the Fortran references are in
`validation/examples/*.summary.json`, generated with `simc_compare --dump-summary <run>.hist out.json`.

```json
{
  "schema": "simc_run_summary", "schema_version": 1,
  "run":      {"tag": "eep_h_elastic_rad_A", "producer": "simc++", "random_seed": 12345},
  "counters": {"ngen": 200000, "ntried": 4150725, "ncontribute": 200000, "nco_no_rad_prot": 0},
  "normalization": {"wtcontr": 207.35421, "sigcc_central": 5.79704e-4, "sigcc_average": 5.94383e-4,
                    "charge_mC": 1.0, "targetfac": 2.67954e-9, "luminosity_ub": 3.73199e8,
                    "luminosity_GeV2": 1453.16, "genvol": 0.0224744, "normfac": 404142.0},
  "flags":    {"doing_phsp": false, "using_rad": true},
  "stops":    {"hms":  {"trials": 3628935, "slit_hor": 2789053, "...": 0, "successes": 200594},
               "shms": {"trials": 200594, "...": 0, "successes": 200000}},
  "tree":     {"file": "eep_h_elastic_rad_A.root", "name": "h10", "weight_branch": "Weight"}
}
```

* Counters are integers. `ngen` has the Fortran sign convention (ngen<0 means |ngen| tries).
* Normalisation values are plain numbers, or `{"value": x}`, in the .hist units:
  `wtcontr = wtcontribute/nevent` (simc.f:920), luminosity in ub^-1 and GeV^2.
* Accepted aliases: `Ntried`, `ncontribute_no_rad_proton`, `central_sigcc`, `average_sigcc`,
  `charge`, `luminosity`, `seed`, `stop_counters`, `root_file`, `tree_file`.
* `stops.<spec>` is lower-case (`hms sos hrsr hrsl shms calo`). Its counter keys follow the
  mapping above. The keys `trials` and `successes` are required for the partition test.
* A relative `tree.file` is resolved relative to the JSON file.
* If `producer` contains "fortran" or "hist2json", values are assumed to carry the .hist print
  precision (e16.6, or e16.8 for wtcontr) in the deterministic comparisons.
