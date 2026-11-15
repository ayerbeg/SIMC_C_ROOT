# 12. FAQ and troubleshooting

## The run stops with `bad targ.can value`

The input uses legacy `.` names (`targ.A`, `spec.e.P`, ...). They are ignored, as by the Fortran, so the target is left at
its defaults. Convert the file and its `extra_dbase_file` with `convert_inputfile.sh`
([4.1.3](04_inputs.md#413-legacy--names-and-convert_inputfilesh)). The warning above the STOP lists the ignored names.
For a liquid target (Z < 2.4), `targ%can` must be 1, 2 or 3.

## `simc++: Fortran STOP during the setup: STOP ...` (exit 3)

SIMC++ reproduces the Fortran input checks. Common ones:

| message | cause |
|---|---|
| `I dont know what phi should be for the electron arm` (or hadron) | `electron_arm`/`hadron_arm` missing or not one of 1-5, 7, 8 |
| `Pi- production from Hydrogen not allowed!`, `Coherent production from Hydrogen/Deuterium not allowed!`, `Sigma- production from Hydrogen not allowed!` | reaction and target do not fit (`which_pion`, `which_kaon`, `targ%A`) |
| `pizero_ngamma not set correctly(must be 1 or 2), stopping` | π0 run without `pizero_ngamma` |
| `I don't have ANY idea what (e,e'm)X we're doing!!!` | `doing_semi` without `doing_pion` or `doing_kaon` |
| `Moron! one_tail>3 turns radiation off, but using_rad wants it on.` | `one_tail` and `using_rad` |
| `Idiot! check setting of deForest_flag` | \|deForest_flag\| > 1 |
| `BAD targ.angle ...` | solid target turned by more than about 58° |

## `simc++: input error` or `configuration error` (exit 2)

The input was refused before the run: arm code 6 (use 5 for the SHMS), π0 from the neutron on hydrogen, `doing_decay` for a
hadron that is not a π or K, rho on A ≥ 4 or combined with another reaction, `using_tgt_field` without `tgt_field_file`,
`--threads` with `--fortran-compat` or RANLUX, an unreadable JSON file. The message says which.

## The run never ends

With `ngen > 0`, SIMC++ (like the Fortran) runs until it has `ngen` **accepted** events. If almost nothing is accepted
(wrong angles or momenta, wide-open `SPedge` defaults, a kinematics outside the acceptance), it runs for a very long time.
Test with a negative `ngen` (`--ngen -10000` = 10000 tries) and look at Ncontribute and the `.geni` stop counters. The
progress lines ("... successes so far") show the rate.

## All weights or the luminosity are 0

Check `targ%abundancy` (in %, 100 for a pure target: the default 0 gives zero luminosity), `EXPER%charge`, and for
A(e,e'p) with Benhar spectral functions `transparency`. Delta runs have no cross section at all
([5.5](05_reactions.md#55-eepπ-in-the-delta-region-doing_delta)).

## Data files not found

`simc++ --print-data-dir` shows where the tables are looked for. Set `SIMC_DATA_DIR` to a directory with the layout of
`data/` (`matrices/`, `pion/`, `kaon/`, `sidis/`, `spectral/`, `field/`), or use `--data-dir` for the spectral and field
files ([2.4](02_installation.md#24-where-the-data-files-are-found)). An installed `simc++` uses `$PREFIX/share/simc/data`.

## `error while loading shared libraries` or crashes in ROOT at start-up

`LD_LIBRARY_PATH` must contain the `lib` of the ROOT 6.40 that SIMC++ was built with, and **not** a directory with an older
ROOT or libtbb (on the development machine: `/usr/local/lib`) ([2.1](02_installation.md#21-requirements)).

## NaN or constant columns in the tree

| column | when | why |
|---|---|---|
| `*xfp *xpfp *yfp *ypfp *ytar` = 0 | `spect_mode` ≠ 0, for the switched-off arm | no spectrometer simulated ([6.4](06_spectrometers.md#64-spect_mode-switching-the-spectrometer-monte-carlo-off)) |
| `pfermi` NaN | SIDIS/rho on H or A ≥ 3 with `--fortran-compat` | Fortran 0/0 (0 in the default mode) |
| `th_tarq` NaN | `using_tgt_field` with `targ_pol = 0` | Fortran acos(0/0), both modes |
| `saghai` = 0 | kaon runs with `--fortran-compat` | Fortran precision bug (D-043 F1) |
| `sigcm` garbage / NaN | Delta runs with `--fortran-compat` | D-048 F3 |
| `sigcent`, `centjac` = 0 | SIDIS/rho | never assigned in the Fortran |
| `xcal_gamma*` = −1e10 | π0 runs | the photon missed the calorimeter |
| `Weight`, `siglab` without cross section | Delta runs | D-048 F4 |

`simc_compare` counts and excludes non-finite entries.

## My numbers differ from a Fortran run

* Different random generator and seed: the runs are statistically equivalent, not identical. Compare them with
  `simc_compare` ([chapter 9](09_comparing_with_fortran.md)).
* The default mode corrects Fortran bugs; see [10.1](10_differences_from_fortran.md#101-corrections-of-the-default-mode-class-a)
  and [9.5](09_comparing_with_fortran.md#95-expected-differences-in-the-default-mode).
* For an exact reproduction use RANLUX, `--fortran-compat` and one thread ([8.4](08_running_and_reproducibility.md#84-reproducing-a-fortran-run-exactly)).

## Two runs with `--threads` give different numbers than without

That is expected: N threads use N random streams. A given seed and N are reproducible; `--threads 1` equals the
sequential run ([8.5](08_running_and_reproducibility.md#85-threads---threads-n)).

## Warnings about `pt_b_param`, `sigc_*`, `random_state_file`, unknown names

These inputs are read and ignored, as by the Fortran (or, for `random_state_file`, not supported). Remove them from the
input to silence the warnings. Unknown names are usually typos or legacy names.

## The `.hist` "Input Spectrometer Limits" differ from my input

The Fortran prints the limits **after** widening them by the spectrometer slop; SIMC++ does the same.

## Can I use my own spectrometer, reaction or cross section?

Yes: every stage is a plug-in. See the developer manual (`docs/manual/developer/`) and `docs/EXTENDING.md`; the
`examples/ideal_box` and `examples/out_of_tree_plugin` directories are complete examples.
