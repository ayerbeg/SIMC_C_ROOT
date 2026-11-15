# 3. Quick start

The commands below were run as shown on the commit of this manual. The full transcript is in
[VERIFIED_EXAMPLES.md](VERIFIED_EXAMPLES.md). They assume:

```bash
export ROOTSYS=/home/ayerbe/Builds/ROOT6_install  LD_LIBRARY_PATH=$ROOTSYS/lib
export PATH=$PREFIX/bin:$ROOTSYS/bin:$PATH          # $PREFIX = the installation prefix (chapter 2)
export SIMC_FORTRAN=/home/ayerbe/NewDisk/simc_gfortran   # a simc_gfortran tree, for its infiles/
```

In a build tree, replace `simc++` with `build/bin/simc++`.

## 3.1 Run a Fortran input file

`eep_h_elastic_rad_A.inp` is H(e,e'p) at Ebeam = 10.6 GeV with the SHMS on the electron and the HMS on the proton,
radiative corrections on, `random_seed = 12345`, 200000 events requested. For a quick look, ask for 5000:

```
$ simc++ $SIMC_FORTRAN/infiles/eep_h_elastic_rad_A.inp -o out --ngen 5000 > out.log; echo "exit status $?"
exit status 0
$ grep -E "fortran_compat|random seed" out.log
 Use random seed = 12345
 SIMC++: fortran_compat = false (clear Fortran bugs corrected, docs/FORTRAN_DIFFERENCES.md)
```

The log has the same set-up messages and progress lines as the Fortran ("Generating Event ... successes so far -
Monitor: ..."). It ends with `Come back soon!`. A run takes about 16 µs per try on one core for this input; here
102995 tries gave 5000 accepted events ([chapter 11](11_performance.md)).

## 3.2 What it writes

```
$ ls out
eep_h_elastic_rad_A.gen
eep_h_elastic_rad_A.geni
eep_h_elastic_rad_A.gen.root
eep_h_elastic_rad_A.hist
eep_h_elastic_rad_A.root
eep_h_elastic_rad_A.summary.json
eep_h_elastic_rad_A.summary.txt
```

| file | content |
|---|---|
| `.hist` | the Fortran run report: inputs, central event, normalisation (`wtcontr`, `normfac`, luminosity, `genvol`), limits |
| `.geni` | the spectrometer stop counters (where the tries were lost) |
| `.gen` | the Fortran diagnostic histograms (text); `.gen.root` the same as ROOT TH1D |
| `.root` | the event tree `h10`, one entry per accepted event (written when `Nntu > 0`) |
| `.summary.json`, `.summary.txt` | the run summary in machine-readable form (the input of `simc_compare`) |

The normalisation is in the `.hist`, as in the Fortran:

```
$ grep -E "^ *(Ngen|Ntried|Ncontribute)|wtcontr|normfac" out/eep_h_elastic_rad_A.hist
            Ngen (request) =       5000
            Ntried         =     102995
            Ncontribute    =       5000
               MeV: wtcontr=   0.20906661E+03
                   normfac =     0.407176E+06
```

[Chapter 7](07_outputs.md) explains every output and every branch.

## 3.3 Look at the tree in ROOT

The `Weight` branch is the event weight. As in the Fortran, `Weight` × normfac / nevent is the yield of one event for the
simulated charge (nevent = Ncontribute when ngen > 0) ([7.1.1](07_outputs.md#711-weights-and-normalisation)). A small macro:

```cpp
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
```

```
$ root -l -b -q peek_h10.C
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
```

Here the electron is in the SHMS (`electron_arm = 5`, a beam-left spectrometer), so `hs*` holds the proton (HMS) and `ss*`
the electron ([7.1.2](07_outputs.md#712-arm-prefixes-hs-and-ss)). Event 2 has radiated: W = 1.47 GeV, Em = 0.86 GeV.

## 3.4 The same run from JSON

```
$ inp2json -o run.json $SIMC_FORTRAN/infiles/eep_h_elastic_rad_A.inp
$ simc++ run.json -o out_json --ngen 5000 > out_json.log && diff <(sed 2,3d out/eep_h_elastic_rad_A.hist) <(sed 2,3d out_json/run.hist) && echo ".hist identical (apart from the time lines)"
.hist identical (apart from the time lines)
$ simc_root_identity out/eep_h_elastic_rad_A.root out_json/run.root
  h10                      TTree  5000 entries x 46 leaves: identical
out/eep_h_elastic_rad_A.root vs out_json/run.root: IDENTICAL
```

(The output names follow the input file name, here `run`; use `--run-name` or `"run": {"name": ...}` to change them.)
Edit `run.json` to change parameters, add a `pipeline` section to change plug-ins, or add `"run": {"fortran_compat":
true}` for the exact Fortran mode ([chapter 4](04_inputs.md)).

## 3.5 Next steps

* More events, faster: `--threads N` ([chapter 8](08_running_and_reproducibility.md), [chapter 11](11_performance.md)).
* Reproduce a Fortran run exactly: [8.4](08_running_and_reproducibility.md#84-reproducing-a-fortran-run-exactly).
* Compare with a Fortran run statistically: [chapter 9](09_comparing_with_fortran.md).
* An old input stops with `bad targ.can value`: [4.1.3](04_inputs.md#413-legacy--names-and-convert_inputfilesh).
