# 7. Following upstream SIMC

SIMC++ follows the official Fortran SIMC, [JeffersonLab/simc_gfortran](https://github.com/JeffersonLab/simc_gfortran).
The file `UPSTREAM_SIMC` at the top of the source tree records the upstream commit that `--fortran-compat` reproduces,
together with the deliberate exceptions. The translation itself (chapters 1-6) was made and proved against commit
294f662 (2026-01-16). This chapter describes the first update, to commit 59ccdf6 (2026-10-01; decisions D-080 to
D-082), and the procedure for the next updates.

## 7.1 The update 294f662 → 59ccdf6

Upstream added 27 commits between the two versions (D-080). The Fortran impact analysis (`reference/U.F1_report.md`) built the
new version in a separate checkout (`simc_gfortran_upstream`), ran it on the existing inputs and on new ones, and
classified each commit. The frozen tree `simc_gfortran` stays at 294f662: all the references of chapters 1-6 depend on
it.

| upstream commit(s) | change | effect on existing inputs | SIMC++ (task) |
|---|---|---|---|
| a26f085, f48a112, 75b5586, 5a71480, b93f558, 805083c, e9895ec | first nuclear-elastic attempt, reverted | none | nothing to port |
| 4df7351, 60c2047 | A(e,e')A nuclear elastic (`doing_nuc_elast`, nuclear_elastic.f); new `.hist` line `doing_nuc_elast= F` | `.hist` line only | ported as at 20d46fc (U.1; findings U-A1, U-B2, U-DEV1) |
| 7f252d6, 4fa5ae1, ddbbb3e, e3f70c5 | CTP portability, .gitignore | none | nothing to port |
| 20d46fc | the STOP for H which_pion = 3 also tests `doing_pizero` (dbase.f:832) | H(e,e'π−)Δ++ runs | already corrected (D-038 F1); compat follows it (U.R) |
| e4dd175, d693325 | one ntuple layout for exclusive, SIDIS and rho; vertex z/pt² always computed; pfermi guard | layout of every tree | upstream layout = default, `--tree-layout legacy` (U.2); z/pt², pfermi (U.1) |
| 4446170, 3df944d | STOP if `worksim/` is missing | none (SIMC++ creates the output directory) | nothing to port |
| a85dc6b, 40e2867, e189bda | Benhar spectral function for Cu-64 (`benharsf_64.dat`) | A = 64 runs that use the SF table | ported (U.1) |
| 39d7de0, 73f5134 | which_pion 2/3: Z/N factors; `tgtweight` only for which_pion 0/10 and 1/11 | A ≥ 3 which_pion 2/3: weight × 1/Z | not followed (U-A3, issue #75) |
| 1fae2d7, 26521e4 | transp.f:158: K → ππ daughter mass Mk → Mpi | kaon runs with decay | already corrected (D-006); compat follows it (U.R) |
| 65d140b, 59ccdf6 | util/recon_hcana: `W_recon` sentinel 1e38 | not SIMC output | nothing to port |

### What changed in SIMC++

* **U.1**: the Cu-64 spectral function, the nuclear-elastic reaction (`doing_nuc_elast`, user manual chapter 5), the
  `.hist` line, z/pt² for every reaction, the SIDIS z > 1 warning, pfermi = 0 instead of NaN in both modes.
* **U.2**: the 59ccdf6 tree layout is the default; `--tree-layout legacy` (JSON `"output": {"layout": "legacy"}`) writes
  the 294f662 layout (user manual chapter 7).
* **U.R**: `--fortran-compat` now reproduces 59ccdf6 (`UPSTREAM_SIMC`). The option `--compat-target 294f662` (JSON
  `"run": {"compat_target": "294f662"}`) keeps the old behaviour of the two findings that upstream fixed, so the
  references of the frozen tree can still be reproduced bit for bit. The test harnesses choose the target from the
  reference: a `.hist` with the line `doing_nuc_elast=` is a 59ccdf6 reference, otherwise 294f662
  (`tests/l1/ranlux_fullrun.sh`, `validation/scripts/p15_shipped_regression.sh`); the compat regression baselines are
  in the 294f662 behaviour (`tests/regression/output_identity.sh` adds `--compat-target 294f662`).

### Findings fixed upstream

Three findings of the translation were fixed by upstream itself (D-081). The catalogue marks them "fixed upstream in
&lt;commit&gt;" (chapter 6, §6.5).

| finding | 294f662 | fixed in | SIMC++ default | `--fortran-compat` (59ccdf6) | `--compat-target 294f662` |
|---|---|---|---|---|---|
| D-006 kaon decay mass | transp.f:158 `m_final = Mk` | 1fae2d7 | Mpi | Mpi | Mk |
| D-038 F1 Δ++ STOP | dbase.f:823-826 STOP without `doing_pizero` | 20d46fc | runs | runs | STOP |
| D-051 F3 SIDIS pfermi NaN | results_write.f:208-209 0/0 | e4dd175 | 0 | 0 | 0 (since U.1) |

### New findings in the upstream code

| id | Category | what | SIMC++ |
|---|---|---|---|
| U-A1 | A | e4dd175 rebuilt event.f/NtupleInit.f/results_write.f from a copy older than 4df7351 and lost the nuclear-elastic lines: the shipped `ee_calcium48.inp` never succeeds | the complete 20d46fc code is kept, in both modes |
| U-A2 | A | kaon with Nntu = 1 STOPs at the first event (66 values written, 65 tags) | no STOP: the kaon tree has 66 columns |
| U-A4 | A (latent) | π0 with the target field writes 84 values into arrays of 80 | the tree is sized from the layout |
| U-A3 | B | 39d7de0: for A ≥ 3, which_pion 2/3 lose the factor Z (weight × 1/Z) | old weight kept in both modes, pending the author's answer |
| U-B1 | B | the new eep columns Wi and epsiloni are always 0 | kept |
| U-B2 | B | the nuclear-elastic model (W1 = 0, Gaussian form factor, W with Mp) | kept |
| U-D1 | D | `.hist` line, z/pt² always computed, SIDIS z > 1 warning, dead init.f:344 clause | kept |
| U-DEV1 | DEV | `doing_nuc_elast` combined with another reaction | refused at set-up |

U-A1, U-A2, U-A3 and U-A4 were reported upstream in
[JeffersonLab/simc_gfortran#75](https://github.com/JeffersonLab/simc_gfortran/issues/75) (2026-11-19). The pion 1/Z
decision (U-A3) waits for the answer there; the user decides when it comes.

### Exceptions: where `--fortran-compat` does not reproduce 59ccdf6

* **U-A1 nuclear elastic**: SIMC++ keeps the complete 20d46fc code (the 59ccdf6 code never produces an event). Reference:
  `reference/runs/up_ee_calcium48_20d46fc`.
* **U-A2, U-A4**: no STOP and no overflow; SIMC++ writes the tree.
* **U-A3 pion 1/Z**: the 294f662 weight for A ≥ 3, which_pion 2/3 (the runs `up_upf1_pion_c12_wp{2,3}_A` differ by the
  factor 1/6 for C-12 and are not identity tests).

### Validation of the compat target

RANLUX + `--fortran-compat` full runs against the 59ccdf6 runs (`l1_simcrun_ranlux_up_*`, label `upstream_layout`):
`.hist` (except the time lines), `.geni`, `.gen` text-identical and the whole tree bitwise identical. U.R added
`up_upf1_pion_h_wp3_A` (H(e,e'π−)Δ++, which 294f662 STOPs) and `up_upf1_kaon_h_k1_rad_dec1_nontu` (kaon decay in
flight, Nntu = 0 because of U-A2: `.hist`/`.geni`/`.gen` only). The 294f662 references keep passing with
`--compat-target 294f662`.

## 7.2 How to follow the next upstream changes

1. **Check.** `tools/check_upstream_simc.sh` lists the upstream commits after the one in `UPSTREAM_SIMC` and the files
   they touch ("Up to date." when there are none).
2. **Fortran impact analysis** (fortran-reference): update the upstream checkout, build it, run the existing inputs with
   the old and the new binary, and classify each commit: no effect, new feature, changed behaviour, fix of a known
   finding, new finding. Produce the references (`reference/runs/up_*`, L1 tables) and a report like
   `reference/U.F1_report.md`.
3. **Port** (one coherent group per task): follow the new Fortran; keep the corrections of the catalogue. Upstream fixes
   of known findings: mark the finding "fixed upstream in &lt;commit&gt;" (D-081) and make `--fortran-compat` follow the fix,
   keeping the old behaviour behind `--compat-target` when old references need it. New upstream bugs: new findings in
   FORTRAN_DIFFERENCES.md (`U-*`), reported upstream.
4. **Validate**: RANLUX full-run identity against the new references, L1 tables, the regression baselines (regenerate
   them only when the change is proved), `ctest -L docs`, `tools/check_all.sh`.
5. **Record**: update `UPSTREAM_SIMC` (full sha, date, the exceptions); `tools/check_upstream_simc.sh` must then report
   0 newer commits. Update this chapter, the catalogue and the manuals (both editions) and rebuild the PDFs.
