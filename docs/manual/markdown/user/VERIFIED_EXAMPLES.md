# Verified examples (transcript)

Every command shown in this manual was run as recorded below. Code and manual text: commit `b045ce6` (branch
`translation`); this transcript file was added in the next commit. Toolchain: GCC 16.2.1, ROOT 6.40.04, glibc 2.43,
Release build installed with `cmake --install`.

**Later wording change.** Three diagnostic messages were reworded afterwards ("Fortran STOP during the setup" → "STOP during setup (legacy SIMC behaviour)", "IGNORED, as by the Fortran" → "IGNORED (legacy SIMC behaviour; ...)", and the corrected-mode banner). Those lines below show the current wording; they were re-run on the final commit. No result changed.

**Paths.** The outputs are real, but the machine-specific paths of the run were rewritten into the layout of
[2.1](02_installation.md#21-directory-layout-used-in-this-manual): `~/src/simcpp` ($SIMCPP_SRC), `~/build/simcpp`
($SIMCPP_BUILD), `~/build/simcpp-validation` (the comparison tool's build), `~/opt/simcpp` ($SIMCPP_PREFIX),
`~/src/simc_gfortran` ($SIMC_FORTRAN), `~/simc_reference` ($SIMC_REFERENCE), `~/opt/root` ($ROOTSYS) and `~/work`
(the working directory). The environment was:

```
export ROOTSYS=~/opt/root LD_LIBRARY_PATH=$ROOTSYS/lib
export SIMCPP_SRC=~/src/simcpp SIMCPP_BUILD=~/build/simcpp SIMCPP_PREFIX=~/opt/simcpp
export SIMC_FORTRAN=~/src/simc_gfortran SIMC_REFERENCE=~/simc_reference
export EX=$SIMCPP_PREFIX/share/simc/docs/manual/markdown/user/examples
export PATH=$SIMCPP_PREFIX/bin:$ROOTSYS/bin:/usr/bin:/bin
```

Each block ran in a fresh `bash -c` in one working directory; its output and non-zero exit status follow it. The build
commands of chapter 2 were run with this layout without compiler warnings. The full test suite is run by
`tools/check_all.sh` on the final commit. The `ctest -L docs` run at the end was made before this file existed, so the
link test then reported only the not-yet-existing `VERIFIED_EXAMPLES.md` links, and `docs_no_machine_paths` reported
files of other documents that were still being edited. The `$SIMC_FORTRAN` and `$SIMC_REFERENCE` commands (chapter 8.4
and `docs/translation/markdown/user_comparing_with_fortran.md`) need the optional legacy data.

```
$ simc++ --help
usage: simc++ <file.inp|file.json> [-o OUTDIR] [--seed N] [--ngen N] [--pipeline FILE]
               [--run-name NAME] [--data-dir DIR] [--infiles-dir DIR] [--quiet] [--print-pipeline]
               [--fortran-compat] [--threads N] [--print-data-dir]
  --fortran-compat  exact Fortran behaviour (known Fortran bugs kept; bitwise validation)
  --threads N       event loop on N threads (independent mt64 streams per thread; a given seed and N
                    always give the same output; not with --fortran-compat or rng ranlux). Default 1.
  --print-data-dir  print the default data directory (env SIMC_DATA_DIR overrides it) and exit

$ simc++ --print-data-dir
~/opt/simcpp/share/simc/data

$ SIMC_DATA_DIR=/data/simc simc++ --print-data-dir
/data/simc

$ simc++ $EX/h_elastic_rad.inp -o out --ngen 5000 > out.log; echo "exit status $?"
exit status 0

$ ls out
h_elastic_rad.gen
h_elastic_rad.geni
h_elastic_rad.gen.root
h_elastic_rad.hist
h_elastic_rad.root
h_elastic_rad.summary.json
h_elastic_rad.summary.txt

$ grep -E "fortran_compat|random seed" out.log
 Use random seed = 12345
 SIMC++: fortran_compat = false (clear legacy-SIMC bugs corrected, docs/translation/markdown/06_catalogue.md)

$ grep -E "^ *(Ngen|Ntried|Ncontribute)|wtcontr|normfac" out/h_elastic_rad.hist
            Ngen (request) =       5000
            Ntried         =     102995
            Ncontribute    =       5000
               MeV: wtcontr=   0.20906661E+03
                   normfac =     0.407176E+06       

$ cat > peek_h10.C <<'MACRO'
  // peek_h10.C: open a SIMC++ tree, print its size and the weighted Em distribution
  void peek_h10(const char* file = "out/h_elastic_rad.root") {
    TFile f(file);
    TTree* t = (TTree*)f.Get("h10");
    std::cout << t->GetEntries() << " events, " << t->GetListOfBranches()->GetEntries() << " branches" << std::endl;
    t->Scan("hsdelta:ssdelta:W:Em:Weight", "", "", 5);
    TH1D h("h", "E_{m} (weighted);E_{m} [GeV]", 100, -0.02, 0.2);
    t->Project("h", "Em", "Weight");
    std::cout << "sum of weights in [-0.02, 0.2] GeV: " << h.Integral() << ", mean Em " << h.GetMean() << " GeV" << std::endl;
  }
  MACRO
  root -l -b -q peek_h10.C
Processing peek_h10.C...
5000 events, 46 branches
************************************************************************
*    Row   *   hsdelta *   ssdelta *         W *        Em *    Weight *
************************************************************************
*        0 * 6.8000464 * -2.967325 * 1.0069015 * 0.0488890 * 0.0002927 *
*        1 * 2.4381718 * -0.939918 * 0.9532427 * 0.0063532 * 0.0004068 *
*        2 * -2.744465 * -9.402965 * 1.4675965 * 0.8642277 * 0.0004998 *
*        3 * -1.449645 * 0.4658077 * 0.9453506 * 0.0016173 * 0.0005580 *
*        4 * 0.9987937 * -1.216251 * 0.9978361 * 0.0712616 * 0.0004664 *
************************************************************************
sum of weights in [-0.02, 0.2] GeV: 2.24626, mean Em 0.0176472 GeV

$ inp2json -o run.json $EX/h_elastic_rad.inp && head -c 300 run.json; echo
{
  "_source": {
    "files": [
      "~/opt/simcpp/share/simc/docs/manual/markdown/user/examples/h_elastic_rad.inp"
    ],
    "generator": "inp2json (SIMC++ P1.3)",
    "not_registered": {},
    "note": "parameters = every dbase.f:regallvars variable, in .inp unit

$ simc++ run.json -o out_json --ngen 5000 > out_json.log && diff <(sed 2,3d out/h_elastic_rad.hist) <(sed 2,3d out_json/run.hist) && echo ".hist identical (apart from the time lines)"
.hist identical (apart from the time lines)

$ simc_root_identity out/h_elastic_rad.root out_json/run.root
  h10                      TTree  5000 entries x 46 leaves: identical
out/h_elastic_rad.root vs out_json/run.root: IDENTICAL

$ simc++ $SIMCPP_PREFIX/share/simc/data/config/h_elastic_rad.json -o out_cfg --ngen 5000 > out_cfg.log; echo "exit status $?"; ls out_cfg | head -3
exit status 0
h_elastic_rad.gen
h_elastic_rad.geni
h_elastic_rad.gen.root

$ $SIMCPP_BUILD/bin/simc++ $SIMCPP_SRC/docs/manual/markdown/user/examples/h_elastic_rad.inp -o out_readme --ngen 5000 > /dev/null; echo "exit status $?"
exit status 0

$ simc++ run.json -o out_t4 --ngen 5000 --threads 4 > out_t4.log; echo "exit status $?"; grep -E "Ntried|Ncontribute" out_t4/run.hist
exit status 0
            Ntried         =     102516
            Ncontribute    =       5000

$ simc++ run.json -o out_t4b --ngen 5000 --threads 4 > /dev/null && simc_root_identity out_t4/run.root out_t4b/run.root | tail -1
out_t4/run.root vs out_t4b/run.root: IDENTICAL

$ simc++ run.json -o out_compat --ngen 5000 --threads 4 --fortran-compat
simc++: --threads 4 (more than one event-loop thread) cannot be used with --fortran-compat: the exact Fortran reproduction needs the single random stream in the Fortran order (run with --threads 1)
[exit 2]

$ simc++ run.json -o out_ideal --ngen 200 --pipeline $SIMCPP_PREFIX/share/simc/examples/ideal_box/hadron_arm.json --print-pipeline --quiet | sed -n '/"spectrometers"/,/^  }/p'
  "spectrometers": {
    "electron": {
      "type": ""
    },
    "hadron": {
      "delta_max": 8.0,
      "delta_min": -8.0,
      "type": "ideal_box",
      "xptar_max": 0.05,
      "xptar_min": -0.05,
      "yptar_max": 0.03,
      "yptar_min": -0.03
    }
  },

$ simc++ $EX/legacy_names.inp -o legacy > legacy.log 2> legacy.err; echo "exit status $?"; cat legacy.err
exit status 3
simc++: WARNING: 46 input variable(s) are not SIMC variables and are IGNORED (legacy SIMC behaviour; their defaults are used):
  EXPER.charge, spec.e.P, spec.e.theta, spec.p.P, spec.p.theta, targ.A
  targ.Z, targ.mass_amu, targ.mrec_amu, targ.rho, targ.thick, targ.angle
  targ.abundancy, targ.can, SPedge.e.delta.min, SPedge.e.delta.max, SPedge.e.yptar.min, SPedge.e.yptar.max
  SPedge.e.xptar.min, SPedge.e.xptar.max, SPedge.p.delta.min, SPedge.p.delta.max, SPedge.p.yptar.min, SPedge.p.yptar.max
  SPedge.p.xptar.min, SPedge.p.xptar.max, gen.xwid, gen.ywid, targ.fr_pattern, targ.fr1
  targ.fr2, targ.xoffset, targ.yoffset, targ.zoffset, spec.e.offset.x, spec.e.offset.y
  spec.e.offset.z, spec.e.offset.xptar, spec.e.offset.yptar, spec.p.offset.x, spec.p.offset.y, spec.p.offset.z
  spec.p.offset.xptar, spec.p.offset.yptar, cuts.Em.min, cuts.Em.max
  Legacy '.' names (targ.A, spec.e.P, ...) must be converted to '%' names with simc_gfortran/infiles/convert_inputfile.sh.
simc++: STOP during setup (legacy SIMC behaviour): STOP bad targ.can value
simc++: note: the input has 46 ignored variable(s) (listed above), so their defaults were used — the likely cause of this STOP. Convert legacy '.' names with simc_gfortran/infiles/convert_inputfile.sh.

$ sed -E ':a; s/^(\s*(targ|spec|SPedge|gen|cuts|EXPER)[A-Za-z_%]*)\./\1%/; ta' $EX/legacy_names.inp > converted.inp
  simc++ converted.inp -o converted_out > converted.log 2>&1; echo "exit status $?"; grep -c "IGNORED" converted.log; grep -E "Ntried|Ncontribute" converted_out/converted.hist
exit status 0
0
            Ntried         =    4158458
            Ncontribute    =     200000

$ echo '{"pipeline": {"rng": {"type": "ranlux"}}}' > ranlux.json
  R=$SIMC_REFERENCE/runs/eep_h_test_pctnames
  simc++ $R/eep_h_test_pctnames.inp -o exact --fortran-compat --pipeline ranlux.json > exact.log; echo "exit status $?"
  diff <(sed 2,3d exact/eep_h_test_pctnames.hist) <(sed 2,3d $R/eep_h_test_pctnames.hist) && echo ".hist: identical apart from the time lines"
  cmp exact/eep_h_test_pctnames.geni $R/eep_h_test_pctnames.geni && cmp exact/eep_h_test_pctnames.gen $R/eep_h_test_pctnames.gen && echo ".geni, .gen: identical"
  simc_root_identity exact/eep_h_test_pctnames.root $R/eep_h_test_pctnames.root | tail -1
exit status 0
.hist: identical apart from the time lines
.geni, .gen: identical
exact/eep_h_test_pctnames.root vs ~/simc_reference/runs/eep_h_test_pctnames/eep_h_test_pctnames.root: IDENTICAL

$ R=$SIMC_REFERENCE/runs/eep_h_test_pctnames
  simc++ $R/eep_h_test_pctnames.inp -o exact_corr --pipeline ranlux.json > /dev/null; diff <(sed 2,3d exact_corr/eep_h_test_pctnames.hist) <(sed 2,3d $R/eep_h_test_pctnames.hist) | head -12
65c65
<                using_Eloss= T        using_Coulomb= F        deForest_flag= 0
---
>                using_Eloss= T        using_Coulomb= F        deForest_flag= F
79c79
<       use_expon=  0   include_hard=  T    calc_spence=  T
---
>       use_expon=  F   include_hard=  T    calc_spence=  T

$ cat > arm6.json <<'JSON'
  {"parameters": {"electron_arm": 6, "hadron_arm": 1}}
  JSON
  simc++ arm6.json -o arm6
simc++: input error: arm6.json: parameters: electron_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS
[exit 2]

$ python3 -c '
  import json
  j = json.load(open("run.json"))
  j["parameters"]["ngen"] = 2000
  j["parameters"]["Ngne"] = 5           # a typo: reported and ignored
  j["run"] = {"name": "edited", "output_dir": "edited_out", "fortran_compat": False}
  json.dump(j, open("edited.json", "w"), indent=1)'
  simc++ edited.json > edited.log; echo "exit status $?"; ls edited_out
simc++: warning: unknown parameter 'Ngne' (not registered in dbase.f:regallvars), ignored
exit status 0
edited.gen
edited.geni
edited.gen.root
edited.hist
edited.root
edited.summary.json
edited.summary.txt

$ simc++ $SIMC_FORTRAN/infiles/eep_h_elastic_rad_A.inp -o l2 --seed 4242 --threads 4 > l2.log; echo "exit status $?"
exit status 0

$ $SIMCPP_BUILD-validation/simc_compare --a $SIMC_REFERENCE/runs/eep_h_elastic_rad_A --b l2/eep_h_elastic_rad_A.summary.json --out cmp \
                         --config $SIMCPP_SRC/validation/config/compare_default.json > cmp.log; echo "exit status $?"
  ls cmp; sed -n '1,25p' cmp/summary.txt
exit status 0
compare.pdf
runA.summary.json
runB.summary.json
summary.json
summary.txt
simc_compare  A=A  B=B
A: ~/simc_reference/runs/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist  [fortran-hist]  tree ~/simc_reference/runs/eep_h_elastic_rad_A/eep_h_elastic_rad_A.root (200000 entries)
B: l2/eep_h_elastic_rad_A.summary.json  [summary-json]  tree l2/eep_h_elastic_rad_A.root (200000 entries)
normalisation: A Weight/Ntried, B Weight/Ntried
VERDICT: PASS   tests counted N=214  failed (p<0.01) F=0  expected 2.14  fraction 0.0000  P(X>=F|N,alpha)=1.000e+00 (need >= 1.000e-03)  gross failures=0 (det mismatch or p<1.000e-05)
worst p = 1.485e-02  (ssdeltai [chi2_ww])   largest |pull| = 1.859  (hms.scin/Ntried)

RUN-LEVEL TESTS
quantity                                      method                   A       errA             B       errB     pull          p  status
acceptance Ncontribute/Ntried                 binomial         0.0481844   0.000105     0.0481577   0.000105    0.179      0.858  ok
yield sum(Weight) [Weight/Ntried]             pull           2.48294e-05   5.64e-08   2.48231e-05   5.64e-08   0.0798      0.936  ok
wtcontr                                       pull               207.354      0.472       207.306      0.472   0.0726      0.942  ok  sigma = wtcontr * relerr: A pass_cuts subset, compat window; B pass_cuts subset, corrected window
AVERAGE.sigcc                                 pull           0.000594383   3.67e-07   0.000594476   3.67e-07   -0.179      0.858  ok
normfac                                       pull                404142        882        403918        881     0.18      0.858  ok (info)  Var(Ntried)/Ntried^2=(1-acc)/Ncontribute; not counted (duplicates acceptance)
sigcc_central                                 det            0.000579704   1.08e-09   0.000579704          -        -          1  ok  tol 1.08e-09
charge_mC                                     det                      1      6e-06             1          -        -          1  ok  tol 6e-06
targetfac                                     det            2.67954e-09   7.68e-15   2.67954e-09          -        -          1  ok  tol 7.68e-15
luminosity_ub                                 det            3.73199e+08        873   3.73199e+08          -        -          1  ok  tol 873
luminosity_GeV2                               det                1453.16    0.00645       1453.16          -        -          1  ok  tol 0.00645
genvol                                        det              0.0224744   7.25e-08     0.0224744          -        -          1  ok  tol 7.25e-08
Ngen (request)                                det                 200000          -        200000          -        -          1  ok
hms.trials/Ntried                             binomial          0.874289   0.000163      0.874455   0.000163   -0.721      0.471  ok
hms.slit_hor/Ntried                           binomial          0.671944    0.00023      0.671843    0.00023    0.308      0.758  ok
hms.slit_vert/Ntried                          binomial          0.146327   0.000173      0.146589   0.000174    -1.07      0.286  ok
hms.slit_oct/Ntried                           binomial        0.00695011   4.08e-05    0.00697492   4.08e-05    -0.43      0.667  ok

$ (cd $SIMCPP_BUILD && ctest -L docs)
Test project ~/build/simcpp
    Start 434: docs_user_manual_parameters
1/7 Test #434: docs_user_manual_parameters ......   Passed    0.11 sec
    Start 435: docs_user_manual_links
2/7 Test #435: docs_user_manual_links ...........   Passed    0.04 sec
    Start 436: docs_no_machine_paths
3/7 Test #436: docs_no_machine_paths ............***Failed    0.03 sec
    Start 437: docs_dev_manual
4/7 Test #437: docs_dev_manual ..................   Passed    0.04 sec
    Start 438: docs_dev_manual_links
5/7 Test #438: docs_dev_manual_links ............   Passed    0.04 sec
    Start 439: docs_translation
6/7 Test #439: docs_translation .................   Passed    0.19 sec
    Start 440: docs_translation_links
7/7 Test #440: docs_translation_links ...........   Passed    0.04 sec

86% tests passed, 1 tests failed out of 7

Label Time Summary:
docs    =   0.48 sec*proc (7 tests)

Total Test time (real) =   0.57 sec

The following tests FAILED:
	436 - docs_no_machine_paths (Failed)                    docs
Errors while running CTest
Output from these tests are in: ~/build/simcpp/Testing/Temporary/LastTest.log
Use "--rerun-failed --output-on-failure" to re-run the failed cases verbosely.
[exit 8]
```
