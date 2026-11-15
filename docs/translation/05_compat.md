# 5. Fortran compatibility layer

Policy (D-031): **clear Fortran bugs are corrected by default**; the run option `fortran_compat` restores the exact
Fortran behaviour so that the bitwise tests (L1 tables, captures, RANLUX full runs, the compat regression baselines) can
still prove the translation. The complete list of differences is [FORTRAN_DIFFERENCES.md](../FORTRAN_DIFFERENCES.md).

## 5.1 Classes and markers

| class | meaning | default mode | `--fortran-compat` | marker in the code |
|---|---|---|---|---|
| a | clear Fortran bug | corrected | Fortran | `// FORTRAN-BUG(D-xxx): original Fortran <file:line>: <code>; corrected because <reason>` |
| a-np | clear bug in a path that is not ported | note only | — | (none left at P15.D1) |
| b | ambiguous, maybe intended (the "expert list") | Fortran | Fortran | `// FORTRAN-DIFF-CANDIDATE(D-xxx): ... Question: ...` |
| c | open question with the user (acceptance mode, D-009/D-028) | Fortran | Fortran | `// FORTRAN-OPEN(D-031 flag 1)` |
| d | faithful odd behaviour, not a bug (precision, formats, dead code, guards) | Fortran | Fortran | `// FORTRAN-QUIRK(...)` |
| dev | deliberate design deviation of SIMC++ (JSON input, exit codes, no `.bin`, arm code 6 refused, ...) | SIMC++ | SIMC++ | FORTRAN_DIFFERENCES.md §6 |

Who decides: findings come from the Fortran reports (`reference/P*.F*_report.md`) or from the C++ reading; the director
classifies them in a decision (`D-0xx`: e.g. D-057 for rho, D-061 SOS/HRS, D-067 target field). An agent that suspects
a bug keeps the Fortran behaviour and reports it; it does not correct on its own.

## 5.2 How the switch reaches the code

There is no global. The flag is set once and passed down explicitly:

```
simc++ --fortran-compat  |  JSON "run": {"fortran_compat": true}
  → driver::RunOptions::fortran_compat                         (CommandLine.cpp)
  → SetupOptions::fortran_compat → RunSetup::fortran_compat    (BuildRunSetup; set-up corrections use it directly)
  → every stage reads StageContext::setup.fortran_compat in init()
  → kernels get it explicitly from their stage:
        function arguments            SpectrometerAngles(..., fortran_compat), rotateHaxis/rotateVaxis(..., fortran_compat)
        argument structs              ExtremeTripArgs::fortran_compat, proton::DeForestState::fortran_compat,
                                      pion::PeepiConfig::fortran_compat, the Saghai/peeK/peepiX/peerho/peedelta configs
        kernel setters                optics::CosyEngine::setFortranCompat (via McHms/McShms/McSos/McHrs::setFortranCompat),
                                      calo::McCalo::setFortranCompat
        driver                        PassCuts(setup, recon), MonteCarlo, PizeroCalorimeter, RunReport, the sinks
                                      (NtupleLayout, FortranHistSink, GenHistSink: report-level corrections)
```

The run log prints the mode (` SIMC++: fortran_compat = true/false ...`), and the summary JSON records it
(`run.fortran_compat`, D-054), which `simc_compare` uses to choose the pass_cuts window (chapter 2).
`--fortran-compat` also refuses `--threads N > 1` (the exact reproduction needs one stream in the Fortran order).

Where the switches live (files that read the flag): src/setup (5 files), src/driver (10), src/output (3),
src/material, src/optics (2), src/radiation (2), src/reactions (+ pion, rho, sidis), src/physics/{proton, nuclear,
pion, kaon, delta, sidis, rho}, src/spectrometers/{hms, shms, sos, hrs, calo}. `grep -rn fortran_compat src include`
lists them.

## 5.3 Adding a correction (class a)

1. **Decision.** The finding must have a D-number and class a (DECISIONS.md). Note the Fortran file:line, what the
   Fortran does, what the correction does, and whether the correction changes the number of random draws (D-058).
2. **Code.** Keep the Fortran path and add the corrected one behind the flag the kernel already receives (or pass the
   flag through its config struct / setter; never read a global):

   <!-- from: src/driver/SimcRun.cpp -->
   ```cpp
   // FORTRAN-BUG(D-006): original Fortran simc.f:237: `recon%p%delta .ge. (SPedge%e%delta%max-slop%MC%p%delta%used)`;
   // corrected because the hadron-arm delta must be tested against the hadron-arm limit SPedge%p%delta%max, like the
   // other five hadron-arm tests (copy-paste of l.231). fortran_compat keeps the electron-arm limit.
   bool PassCuts(const RunSetup& S, const Event& r) {
       const auto& E = S.SPedge;
       const auto& M = S.slop.MC;
       const double p_delta_max = S.fortran_compat ? E.e.delta.max : E.p.delta.max;
   ```

   (`driver::PassCuts`, the first correction of the project.) The compat path must stay bitwise the
   Fortran, including its random draws.
3. **Both-mode test.** Add an item to `tests/fortran_diff/test_fortran_differences.cpp` (a function registered in
   `items()` under the id, e.g. `"D061_F1e_hrs_vdc_rotation"`) that shows **both** behaviours: `fortran_compat = true`
   gives the Fortran result (ideally against a reference table or run), the default gives the corrected one. Add the id
   to `SIMC_FORTRAN_DIFF_ITEMS` in tests/fortran_diff/CMakeLists.txt (CTest `fortran_diff_<id>`, labels
   `fortran_diff;unit`). Corrections visible only in a full run can instead be tested by a RANLUX full run in both modes
   (tests/l1, `ranlux_fullrun.sh`).
4. **The L1 / replay tests keep running in compat** (they pass `fortran_compat = true` or `--fortran-compat`); check
   that they still pass bitwise.
5. **FORTRAN_DIFFERENCES.md row** in §2 (class a): `| id | Fortran file:line | original behaviour | C++ default |
   effect (outputs, size, configurations) | test |`, plus the measured effect in the per-reaction L2 section when the
   correction changes distributions (an L2 run in both modes, `simc_compare` corrected vs Fortran).
6. **Baselines.** A correction changes the corrected-mode regression baselines (`regression_corrected_*` and the
   `*_corrected_*` families): regenerate them (chapter 2 §2.6) and state in the commit which families changed and why.
   The compat baselines must not change.
7. **User manual**: one line in the user-facing summary of corrections.

## 5.4 Ambiguous cases, quirks and the acceptance mode

* Class b: keep the Fortran in both modes, write `FORTRAN-DIFF-CANDIDATE(D-xxx)` with the question, add a §4 row in
  FORTRAN_DIFFERENCES.md. These items form the expert list (P15.X) for the user and a physicist.
* Class d: keep, write `FORTRAN-QUIRK(...)` with the reason (e.g. `FORTRAN-QUIRK(D-012 F2)`), add a §5 row when it
  can surprise a user.
* Class c: the acceptance mode (`doing_phsp`, D-009/D-028) stays Fortran with `FORTRAN-OPEN(D-031 flag 1)` until the
  user answers.
