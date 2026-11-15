# Verified examples (transcript)

Every command shown in this manual was run as recorded below (task P15.D2). Commit of the code and the manual text:
`8db5f23` (branch `translation`); this transcript file was added in the next commit. Toolchain: GCC 16.2.1, ROOT
6.40.04, Fedora 44 (glibc 2.43), Release build installed with `cmake --install build --prefix $PREFIX`.

Environment of the run (`.tmp_P15.D2/verify/run_cmds.sh`; each block runs in a fresh `bash -c` in one work directory,
output and non-zero exit status shown under it):

```
ROOTSYS=/home/ayerbe/Builds/ROOT6_install   LD_LIBRARY_PATH=$ROOTSYS/lib   (not /usr/local/lib)
PREFIX=/home/ayerbe/NewDisk/SIMC++_Project/.tmp_P15.D2/install      PATH=$PREFIX/bin:$ROOTSYS/bin:/usr/bin:/bin
SIMC_SRC=/home/ayerbe/NewDisk/SIMC++_Project/SIMC++_Claude          (the source tree)
SIMC_FORTRAN=/home/ayerbe/NewDisk/simc_gfortran                     SIMC_REF=/home/ayerbe/NewDisk/reference
```

The build commands of chapter 2 were run with these paths (`cmake -S SIMC++_Claude -B build ...`, `cmake --build build
-j8`, the same for `validation` into `build-val`, `cmake --install`), with no compiler warning. The complete test suite
(`ctest`, chapter 2.5) is run by `tools/check_all.sh` on the final commit (its result is in the task report). The
`ctest -L docs` run at the end of the transcript was made before this file existed, so its link test reports the
then-missing `VERIFIED_EXAMPLES.md` links only.

Where the manual abbreviates an output (`...`), the full output is here. The tiny differences between the chapters and this
transcript are only paths printed by the tools.

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
/home/ayerbe/NewDisk/SIMC++_Project/.tmp_P15.D2/install/share/simc/data

$ SIMC_DATA_DIR=/data/simc simc++ --print-data-dir
/data/simc

$ simc++ $SIMC_FORTRAN/infiles/eep_h_elastic_rad_A.inp -o out --ngen 5000 > out.log; echo "exit status $?"
exit status 0

$ ls out
eep_h_elastic_rad_A.gen
eep_h_elastic_rad_A.geni
eep_h_elastic_rad_A.gen.root
eep_h_elastic_rad_A.hist
eep_h_elastic_rad_A.root
eep_h_elastic_rad_A.summary.json
eep_h_elastic_rad_A.summary.txt

$ grep -E "fortran_compat|random seed|Ntried|Ncontribute" out.log
 Use random seed = 12345
 SIMC++: fortran_compat = false (clear Fortran bugs corrected, docs/FORTRAN_DIFFERENCES.md)

$ grep -E "^ *(Ngen|Ntried|Ncontribute)|wtcontr|normfac" out/eep_h_elastic_rad_A.hist
            Ngen (request) =       5000
            Ntried         =     102995
            Ncontribute    =       5000
               MeV: wtcontr=   0.20906661E+03
                   normfac =     0.407176E+06       

$ inp2json -o run.json $SIMC_FORTRAN/infiles/eep_h_elastic_rad_A.inp && head -c 300 run.json; echo
{
  "_source": {
    "files": [
      "/home/ayerbe/NewDisk/simc_gfortran/infiles/eep_h_elastic_rad_A.inp"
    ],
    "generator": "inp2json (SIMC++ P1.3)",
    "not_registered": {},
    "note": "parameters = every dbase.f:regallvars variable, in .inp units; extra_dbase_file is already merged",
    

$ simc++ run.json -o out_json --ngen 5000 > out_json.log && diff <(sed 2,3d out/eep_h_elastic_rad_A.hist) <(sed 2,3d out_json/run.hist) && echo ".hist identical (apart from the time lines)"
.hist identical (apart from the time lines)

$ simc_root_identity out/eep_h_elastic_rad_A.root out_json/run.root
  h10                      TTree  5000 entries x 46 leaves: identical
out/eep_h_elastic_rad_A.root vs out_json/run.root: IDENTICAL

$ cat > peek_h10.C <<'MACRO'
  // peek_h10.C: open a SIMC++ tree, print its size and the weighted Em distribution
  void peek_h10(const char* file = "out/eep_h_elastic_rad_A.root") {
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

$ simc++ run.json -o out_t4 --ngen 5000 --threads 4 > out_t4.log; echo "exit status $?"; grep -E "Ntried|Ncontribute" out_t4/run.hist
exit status 0
            Ntried         =     102516
            Ncontribute    =       5000

$ simc++ run.json -o out_t4b --ngen 5000 --threads 4 > /dev/null && simc_root_identity out_t4/run.root out_t4b/run.root | tail -1
out_t4/run.root vs out_t4b/run.root: IDENTICAL

$ simc++ run.json -o out_compat --ngen 5000 --threads 4 --fortran-compat
simc++: --threads 4 (more than one event-loop thread) cannot be used with --fortran-compat: the exact Fortran reproduction needs the single random stream in the Fortran order (run with --threads 1)
[exit 2]

$ simc++ run.json -o out_ideal --ngen 200 --pipeline $PREFIX/share/simc/examples/ideal_box/hadron_arm.json --print-pipeline --quiet | sed -n '/"spectrometers"/,/^  }/p'
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

$ simc++ $SIMC_FORTRAN/infiles/test_eep_h.inp -o legacy > legacy.log 2> legacy.err; echo "exit status $?"; cat legacy.err; tail -2 legacy.log
exit status 3
simc++: WARNING: 46 input variable(s) are not SIMC variables and are IGNORED, as by the Fortran (their defaults are used):
  EXPER.charge, spec.e.P, spec.e.theta, spec.p.P, spec.p.theta, targ.A
  targ.Z, targ.mass_amu, targ.mrec_amu, targ.rho, targ.thick, targ.angle
  targ.abundancy, targ.can, SPedge.e.delta.min, SPedge.e.delta.max, SPedge.e.yptar.min, SPedge.e.yptar.max
  SPedge.e.xptar.min, SPedge.e.xptar.max, SPedge.p.delta.min, SPedge.p.delta.max, SPedge.p.yptar.min, SPedge.p.yptar.max
  SPedge.p.xptar.min, SPedge.p.xptar.max, gen.xwid, gen.ywid, targ.fr_pattern, targ.fr1
  targ.fr2, targ.xoffset, targ.yoffset, targ.zoffset, spec.e.offset.x, spec.e.offset.y
  spec.e.offset.z, spec.e.offset.xptar, spec.e.offset.yptar, spec.p.offset.x, spec.p.offset.y, spec.p.offset.z
  spec.p.offset.xptar, spec.p.offset.yptar, cuts.Em.min, cuts.Em.max
  Legacy '.' names (targ.A, spec.e.P, ...) must be converted to '%' names with simc_gfortran/infiles/convert_inputfile.sh.
simc++: Fortran STOP during the setup: STOP bad targ.can value
simc++: note: the input has 46 ignored variable(s) (listed above), so their defaults were used — the likely cause of this STOP. Convert legacy '.' names with simc_gfortran/infiles/convert_inputfile.sh.
cuts.Em.max not registered, registering as int
 targ.Mrec replaced by the expectation M - Mp (dbase.f:283-284)

$ mkdir -p converted && cp $SIMC_FORTRAN/infiles/test_eep_h.inp $SIMC_FORTRAN/infiles/extra_hydrogen.inp converted/ && cd converted && for f in test_eep_h.inp extra_hydrogen.inp; do sh $SIMC_FORTRAN/infiles/convert_inputfile.sh $f; done && rm -f temp.inp && grep -c '%' test_eep_h.inp extra_hydrogen.inp
test_eep_h.inp:27
extra_hydrogen.inp:20

$ simc++ converted/test_eep_h.inp -o converted_out > converted.log 2>&1; echo "exit status $?"; grep -c "IGNORED" converted.log; grep -E "Ntried|Ncontribute" converted_out/test_eep_h.hist
exit status 0
0
            Ntried         =       4069
            Ncontribute    =       1000

$ echo '{"pipeline": {"rng": {"type": "ranlux"}}}' > ranlux.json
  R=$SIMC_REF/runs/eep_h_test_pctnames
  simc++ $R/eep_h_test_pctnames.inp -o exact --fortran-compat --pipeline ranlux.json > exact.log; echo "exit status $?"
  diff <(sed 2,3d exact/eep_h_test_pctnames.hist) <(sed 2,3d $R/eep_h_test_pctnames.hist) && echo ".hist: identical apart from the time lines"
  cmp exact/eep_h_test_pctnames.geni $R/eep_h_test_pctnames.geni && cmp exact/eep_h_test_pctnames.gen $R/eep_h_test_pctnames.gen && echo ".geni, .gen: identical"
  simc_root_identity exact/eep_h_test_pctnames.root $R/eep_h_test_pctnames.root | tail -1
exit status 0
.hist: identical apart from the time lines
.geni, .gen: identical
exact/eep_h_test_pctnames.root vs /home/ayerbe/NewDisk/reference/runs/eep_h_test_pctnames/eep_h_test_pctnames.root: IDENTICAL

$ R=$SIMC_REF/runs/eep_h_test_pctnames
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

$ VAL=/home/ayerbe/NewDisk/SIMC++_Project/.tmp_P15.D2/build-val
  $VAL/simc_compare --a $SIMC_REF/runs/eep_h_elastic_rad_A --b l2/eep_h_elastic_rad_A.summary.json --out cmp \
                    --config $SIMC_SRC/validation/config/compare_default.json > cmp.log; echo "exit status $?"
  ls cmp; sed -n '1,25p' cmp/summary.txt
exit status 0
compare.pdf
runA.summary.json
runB.summary.json
summary.json
summary.txt
simc_compare  A=A  B=B
A: /home/ayerbe/NewDisk/reference/runs/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist  [fortran-hist]  tree /home/ayerbe/NewDisk/reference/runs/eep_h_elastic_rad_A/eep_h_elastic_rad_A.root (200000 entries)
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

$ (cd /home/ayerbe/NewDisk/SIMC++_Project/.tmp_P15.D2/build && ctest -L docs)
Test project /home/ayerbe/NewDisk/SIMC++_Project/.tmp_P15.D2/build
    Start 434: docs_user_manual_parameters
1/4 Test #434: docs_user_manual_parameters ......   Passed    0.12 sec
    Start 435: docs_user_manual_links
2/4 Test #435: docs_user_manual_links ...........***Failed    0.04 sec
    Start 436: docs_dev_manual
3/4 Test #436: docs_dev_manual ..................   Passed    0.24 sec
    Start 437: docs_dev_manual_links
4/4 Test #437: docs_dev_manual_links ............   Passed    0.08 sec

75% tests passed, 1 tests failed out of 4

Label Time Summary:
docs    =   0.49 sec*proc (4 tests)

Total Test time (real) =   0.56 sec

The following tests FAILED:
	435 - docs_user_manual_links (Failed)                   docs
Errors while running CTest
Output from these tests are in: /home/ayerbe/NewDisk/SIMC++_Project/.tmp_P15.D2/build/Testing/Temporary/LastTest.log
Use "--rerun-failed --output-on-failure" to re-run the failed cases verbosely.
[exit 8]
```
