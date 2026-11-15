# 4. Numerical conventions

SIMC++ reproduces the Fortran **bit for bit** wherever a Fortran reference exists (L1 tables, captures, RANLUX full
runs). That is only possible if every rule of this chapter is followed. When an L1 test differs in the last bit,
the cause is almost always one of them.

## 4.1 Units

One convention everywhere inside the framework (core/Units.h): MeV, MeV/c, MeV², cm, rad, slopes dx/dz and dy/dz,
g/cm², g/cm³, delta in %. Cross sections are in the units of the Fortran routine they replace (documented per model).
Exceptions:

* **input units** (`SimcParameters`: degrees, mrad, mg/cm², ...) are converted once by the set-up with the dbase.f
  formulas; `RunSetup` holds internal units, `RunSetup::params` the raw input;
* **output units** (ntuple: GeV, mrad, ... as results_write.f; `.hist` formats) are converted by the output code;
* kernels keep the Fortran's own literals where the Fortran writes them (brem.f works in GeV with `/1000.`), so a port
  stays line-by-line comparable.

The named factors of `simc::units` are exact and must be used with the same operation as the Fortran (`x / 1000.`
stays a division, never a multiplication by 0.001), so naming a literal never changes a bit.

## 4.2 Precision: double, real*4, REAL(16)

The Fortran is compiled with `-fdefault-real-8`: every `REAL` and every untyped literal (`0.62`, `1./137.0359895`) is
a **double**. Use double literals of the same decimal value. But `-fdefault-real-8` also promotes `DOUBLE PRECISION`
to **REAL(16)** (binary128), which has visible consequences:

* **REAL(16) literals with a `d` exponent** (`1.d-27`, `1.d6`) make the expression they appear in a binary128
  computation, rounded to double at the assignment. SIMC++ reproduces them with `__float128` (`__extension__ typedef
  __float128 quad;`), building the constant exactly (e.g. 1e-27 = 1/(5²⁷·2²⁷) with an exact integer denominator):
  pion_coll_absorb (src/spectrometers/hms/McHmsColl.cpp), the target factor (src/setup/SetupBuilder.cpp l.170-177),
  the COSY recon guard (src/optics/CosyEngine.cpp l.729), rho_physics.f divisions (src/physics/rho/RhoPhysics.cpp),
  physics_delta.f (src/physics/delta/DeltaModel.cpp), fdss.f (src/physics/sidis/Fdss.cpp), COMPLEX(16) in the kaon
  model (src/physics/kaon/KaonCrossSections.cpp, glibc `hypotf128`). The objdump of the Fortran object
  (`__multf3`, `__divtf3`, `__trunctfdf2`) tells which operations are quad.
* **REAL(16) actual arguments passed to real*8 dummies** (an ABI bug of the Fortran source: the callee reads the
  first 8 bytes of the 16-byte value) give garbage values: the HRS transport lengths and the 45° VDC rotation
  (D-061 F1/F1e: zd = −3.13e232, rotation 0°), the kaon `fint` arguments (D-043), the Delta `sig_blok` arguments
  (D-048 F3). These are class-a bugs: corrected by default, reproduced in `--fortran-compat` (the garbage bit patterns
  are computed exactly, e.g. `McHrs::transpLength`).
* **real*4 (float)** where the Fortran declares it: CERNLIB `lfit` (sums in double from float inputs, results rounded
  to float, optics/Lfit.h) and the HMS/SHMS/SOS/HRS drift-chamber arrays fed to it, so the focal-plane values are
  float values; `COMPLEX*8` of sig_multipole (single-precision arithmetic written out by hand,
  src/physics/pion/SigMultipole.cpp); the kaon `fint` tables; the ntuple (`RootTreeSink` stores `Float_t`, as
  make_root_tree's real*4 `ntup_out`).
* **CTP input literals are floats** (D-018): CTP stores every real literal of a `.inp` through `float`, then evaluates
  expressions in double. `InpReader` rounds every non-integer literal to float before evaluating
  (src/input/InpReader.cpp l.81-88). A JSON input is double (it is a different input format).

## 4.3 Powers

* `x**2` (and `x**2.`) → `x*x` (or `sq(x)`). gcc ≤ 15 folded `std::pow(x, 2.)` into `x*x`, gcc 16 calls libm `pow`
  under `-fmath-errno`, which can differ by 1 ulp. The unit test `lint_no_pow_square` (tests/unit/lint_no_pow_square.sh)
  fails on any `pow(x, 2)` in src/, include/, tools/.
* `x**n` with an integer variable or `integer*2` exponent → `optics::powi(x, n)`, the libgcc `__powidf2` algorithm
  (repeated squaring, `1/y` for negative n; optics/CosyEngine.h). Not `std::pow`.
* real exponents (`x**0.61`, `pow(Q,3.)` in fofa) → `std::pow`, as gfortran calls libm `pow`.

## 4.4 Operation order, NaN signs and Fortran intrinsics

* Keep the **operation order** of the source: `a*b*c` is `(a*b)*c`, sums of dot products are `x+y+z` in that order,
  `-(a)*b` is the product then the sign flip. Do not "simplify" an expression of a kernel.
* When a NaN can reach an output (garbage arguments, D-048 F3; `acos` of > 1; 0/0), the **sign of the NaN** matters for
  bitwise comparisons and follows the x86 rule: an operation with two NaN operands returns the first operand's NaN, and
  `-x` flips the sign bit. Match the operand order of the **object code**, not only of the source: gfortran may commute
  operands, e.g. `exp(-11.5*abs(t))` is compiled as `exp(-(abs(t)*11.5))` (src/physics/pion/PionCrossSections.cpp
  l.33-58 explains the cases found with `objdump -d physics_pion.o`).
* gfortran `MAX(a, b)` / `MIN(a, b)` for reals are `b if b > a or a is NaN` — not `std::max` (which differs for NaN). Use
  the local helpers `fmaxF` / `fminF` (src/reactions/EepReactionBase.cpp, src/driver/RunReportInternal.h).
* `SIGN(1., x)` is `std::copysign(1., x)` (so `SIGN(1., -0.)` = −1, D-067 T10). `NINT` rounds half away from zero.
  Integer division truncates toward zero. Implicit typing: variables starting with `i`–`n` are integers.
* `-fno-automatic`: every local variable is SAVEd. A local that is not assigned on some path keeps the value of the
  previous call (`MonteCarlo::tmpfact_`, `EepReactionBase::krel_`): it becomes a member of the object.
* `-ff2c`: REAL functions return double under the f2c convention; no effect on the values with `-fdefault-real-8`.

## 4.5 libm and toolchain dependence (D-044)

The bitwise references depend on the libm of the machine: the SIMC objects link libm and libgfortran dynamically, so a
glibc update changes their results by a few ulp (2026-10-09: glibc 2.42 → 2.43 changed `tgamma` by ≤ 2 ulp in 49 of 274
table points; every radiative weight downstream changed). Rules:

* The reference toolchain is recorded in `reference/TOOLCHAIN.txt` (glibc/libm md5, libgfortran, compilers); every C++
  baseline README.txt has a `toolchain:` line; `tools/check_all.sh` records the current fingerprint
  (`tools/toolchain_fingerprint.sh`).
* After an OS compiler or libc update, regenerate the references (Fortran side: `fortran_harness/P10F2/regen_all.sh`
  by the fortran-reference agent; C++ baselines: `validation/scripts/regen_baselines.sh`) **before** trusting a bitwise
  failure. `tests/toolchain_guard.sh` prints SAME / DIFFERENT next to a failing L1 test.
* Use the same libm functions as gfortran: `std::tgamma` for the intrinsic `gamma`, `std::exp`, `std::log`, `std::pow`,
  `std::atan2`... The `tgamma` memo of `RadcState` (P15.4) is bitwise neutral because `tgamma` is pure.
* No `-march=native` (FMA contraction), no `-ffast-math`, no `-ffp-contract=fast` (developer manual, chapter 2).

## 4.6 Random-number order invariants

* Every random number is drawn by a stage, through `IRandom`, in the **same order and number as the Fortran**: one
  `uniform()` per `grnd()`, one `gauss1(nsig)` per `gauss1(nsig)`. Draws that the Fortran makes even when their value is
  unused must stay (D-024: the beam-energy grnd is drawn with `dEbeam = 0`; the target energy-loss gaussians with
  `using_Eloss = 0`). Example budget (D-024): H(e,e'p) norad generate = 6 grnd + 12 gauss1; montecarlo success
  = 1 grnd (+410 gauss1 with smearing).
* `gauss1` is the polar method of gauss1.f (`gauss1Polar`, core/random/Gauss1.h): 2 uniforms per attempt, only v1 used,
  rejection if s > 1 or s = 0 or |g| > nsigmax; RANLUX quirks (u1 = 0.5 exactly gives 0) are kept.
* `ranlux` reproduces call_ranlux.f: `sgrnd(seed)` = RLUXGO(3, seed, 0, 0), grnd hands out numbers 1..999 of each
  block of 1000 (the 1000th is never used); every value is a multiple of 2⁻⁴⁸ (exact in double).
* A correction (`FORTRAN-BUG`) that changes the number of draws makes the corrected mode non-replayable against
  Fortran streams; that is accepted (D-058), but `--fortran-compat` must keep the Fortran draws exactly. State in the
  comment whether a correction changes the stream.
* `mt64` (default) gives 53-bit uniforms strictly inside (0,1); the worker streams of `--threads N` are
  `seedStream(seed, k)` (developer manual, chapter 5).

## 4.7 Checking a numerical question

1. Read the Fortran line and the flags (`-fdefault-real-8`, `-fno-automatic`).
2. If precision or operand order is in doubt, disassemble the Fortran object (`objdump -d
   /home/ayerbe/NewDisk/simc_gfortran/<file>.o`, read-only) and look for `__*tf3` (quad), `cvtss2sd` (float),
   `mulsd`/`xorpd` order, calls to `__powidf2`, `pow`, `tgamma`.
3. Compare against the L1 table with the runner's report (max |rel|, bitwise-equal rows, NaN mismatches, worst row).
