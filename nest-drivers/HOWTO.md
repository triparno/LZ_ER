# Building and running the NEST drivers, step by step

NEST is not redistributed with this archive. You clone it, build it, and then
compile the three drivers of this archive against the library it produces. The
drivers stay here; **nothing from this archive goes into the NEST tree**, and
no NEST build file is edited. The whole procedure is below; it assumes a
Unix-like system with git, cmake (>= 3.10), a C++17 compiler and GSL.

If you only want to check the paper's arithmetic, you do not need any of this:
the scan output is already in `../scan/scan24_results.txt`, and every script in
`../analysis/` runs against it. This is needed only to regenerate that file or
to query the response model at new energies.

---

## 1. Get NEST at the right commit

    git clone https://github.com/NESTCollaboration/nest.git
    cd nest
    git checkout 541b643

**Do not use a tagged release.** The Letter's main text cites v2.4.5 and its
supplement discusses v2.4.0, but neither reproduces these numbers. The EFT2024
parameter vector has a component centre above unity, which selects the
log10(Nq) two-component branch of `RecombOmegaER`; that branch was introduced
in v2.4.5beta (December 2025) and the configuration itself was added in commit
307fb27 (27 August 2026). Commit 541b643 is what the paper used.

Confirm the configuration is present before going further:

    grep -n "SetEFT2024Configuration" include/Detectors/LZ_WS2024.hh

You should see the member function declared on class `LZ_Detector_2024`. If
grep finds nothing, the checkout is wrong.

## 2. Build NEST itself

The checkout must be unmodified at this point, and in particular the drivers
of this archive must **not** be inside it. They are compiled separately in
section 3, against the library built here. If anything of theirs is in the
NEST tree, `make` will try to build it, fail on `NEST.hh`, and stop before
producing the library -- see the note at the end of this section.

    cd /path/to/NEST
    mkdir -p build && cd build
    cmake .. && make

This produces `build/libNESTCore.a`. The drivers link against it.

Keep the NEST tree clean: nothing from this archive goes into it. If a
previous attempt added driver sources or `add_executable` targets there, plain
`make` will build them too and die on `NEST.hh` before the library is
produced, which makes the failure look like a problem with NEST. Restore
`CMakeLists.txt`, remove the added sources, and delete `build/` before
re-running cmake -- the cache remembers the targets until you do.

## 3. Build the drivers

The drivers stay in this archive; nothing of theirs goes into the NEST tree.
Use the script, which finds the include directories and the library wherever
your checkout puts them:

    cd /path/to/repo/nest-drivers
    ./build.sh /path/to/NEST

It prints what it found and leaves `lz_tail_test`, `lz_scan` and
`lz_scan_s1s2` beside the sources. Pass the build directory as a second
argument if it is not `/path/to/NEST/build`.

By hand, if you would rather. Do not assemble the include list from the NEST
headers alone: `RandomGen.hh` includes `gcem.hpp`, a header-only dependency
that NEST's CMake fetches into the build tree, and there may be others. Take
the flags CMake itself used:

    NESTSRC=/path/to/NEST
    FLAGS=$(find "$NESTSRC/build" -path '*.dir/flags.make' | grep -i nest | head -1)
    INCS=$(sed -n 's/^CXX_INCLUDES *= *//p' "$FLAGS" | head -1)
    LIB=$(find "$NESTSRC/build" -name 'libNEST*' -print -quit)   # libNESTCore.a
    g++ -std=c++17 -O2 lz_scan_s1s2.cpp -o lz_scan_s1s2 $INCS "$LIB"

`flags.make` is written by CMake's Makefile generator and holds the exact
include flags NEST compiled itself with, vendored dependencies included.

### Do not add them as CMake targets

Do not copy the drivers into `examples/`, add `add_executable` targets and
link `NEST`. It fails twice over.

The headers are not found, because `target_link_libraries(... NEST)` brings in
include directories only if the target exports them as `PUBLIC` or
`INTERFACE`, and at this commit it does not:

    fatal error: NEST.hh: No such file or directory

And the target name is wrong anyway: at this commit the library is
`NESTCore`, not `NEST`.

Worse, once those targets are in NEST's own `CMakeLists.txt`, a plain `make`
intended to build NEST tries to build them too, fails on `NEST.hh`, and stops
before producing the library -- so section 2 breaks as well, and the error
points at a driver source file rather than at the added target. Section 2 has
the recovery steps.

Compiling the drivers separately, as above, avoids all of this. They are three
single translation units with one dependency; they do not need a build
system.

## 4. Check it works

    ./lz_tail_test 1000000

prints the response at the **simulated** point, E = 70.95 keVee, which is
where the scan is anchored. That is not the event-conditioned energy of the
paper, E = 69.84 keVee: the event has N_q = 5161, this point has N_q = 5243,
and the quantities below are correspondingly a little larger. At 10^6 trials
it takes a few seconds. To within sampling error:

    recombProb  0.7157        Variance  43440   (binomial part 903)
    sd(Ne)      208.4         mean(Ne)  1262
    skew(Ne)    0.42
    omega  extracted 0.04648  code form 0.04646  Eq. 4 as printed 0.06853
           extracted/code 1.0003   printed/code 1.4750

The three omegas are the point of the test. *Extracted* comes from the
variance NEST returns, *code form* from Eq. 4 with the component-wise mode
normalisation `RecombOmegaER` applies, and *as printed* from Eq. 4 evaluated
literally with Table S4. The first two must agree to well under a per cent --
the driver says so explicitly if they do not -- and the third is larger by
the factor of 1.475 that Sec. 2 of the paper is about. That last ratio is
the paper's central discrepancy, arrived at here in C++ from the variance
NEST returns, independently of anything in `analysis/`. If the extracted value
comes out near 0.0685 rather than 0.0465, `SetEFT2024Configuration()` did not
take.

`P(Ne<409)` and `P(Ne<269)` both print as zero at 10^6 trials. That is
correct: the leakage is of order 10^-7, so 10^6 trials expect 0.13 counts.
Seeing the tail at all is what the 10^10-trial scan of section 5 is for.

The event-conditioned numbers of the paper -- mean 1238, sd 204.1,
gamma_1 0.423, omega 0.0462 -- are reproduced by `analysis/response.py`,
which evaluates the same model at N_q = 5161 without simulating.

## 5. Regenerate the scan

From this directory, where section 3 left the binaries:

    cd /path/to/repo/nest-drivers
    ./run_scan24.sh 8 ./lz_scan

The first argument is the number of parallel processes, the second the driver.
About two hours on eight cores for 1.084e10 trials across 24 energies. The
script writes `../scan/scan24_results.txt`, anchoring that path on its own
location rather than the working directory.

It replaces the archived copy only if the run validates. It checks the exit
status of every parallel job, and then that the output carries exactly the 24
energies of the job list with none missing or repeated, that every row has its
full thirteen fields, that every field is a finite number, and that the trial
count each row reports is the one that row was given. On any failure the
archived file is left as it was and the candidate is kept beside it, named in
the error. `validate_scan.sh` holds those checks and `run_scan24_smeared.sh`
uses the same ones against its twenty columns.

To satisfy yourself the guards work, `./test_runner_guards.sh` drives the
runner with stub drivers that fail in each of those ways and asserts it
refuses them all. It needs no NEST build and runs in a second.

**The run is deterministic.** run_scan24.sh fixes seeds 201-224 and NEST's
`RandomGen` is a singleton with a deterministic initial state, so a rerun
reproduces the archived numbers exactly, not merely statistically. If your
output differs, something has changed in the configuration or the commit.

## 6. The scan carried through to S1c and S2c

`lz_scan.cpp` stops at the quanta level: it cuts on \(N_e\) and \(N_{ph}\),
so the photon-detection and electron-extraction fluctuations are absent and
the edges of the S1c window and the S2c ceiling are sharp. `lz_scan_s1s2.cpp`
carries the same events through `GetS1` and `GetS2` and applies the cuts in
phd, which is what the detector does.

### 6.1 Find the two vector indices first

`GetS1` and `GetS2` each return a vector whose layout differs between NEST
versions, and the driver will not guess: in scan mode both indices are
required arguments. Run the probe:

    ./lz_scan_s1s2 probe 70

It prints every element of both vectors with its mean, its spread, and its
ratio to the expected corrected pulse area, and marks the candidates:

    <Nph> = 4892.3   ->  g1*<Nph> = 539.91 phd   (expected S1c)
    <Ne>  =  268.6   ->  g2*<Ne>  = 9266.6 phd   (expected S2c)

    idx           mean          sd   /g1<Nph>  integer
      ...
      3       540.1234     23.8011     1.0004       no   <-- candidate s1c_idx

Choose the element that has a non-zero spread, is not integer-valued, and
whose ratio is **closest to one**: integer entries are hit or spike counts
rather than areas, and \(g_1\) and \(g_2\) are defined so that the
*corrected* area equals \(g N\) exactly, leaving the raw area a per cent or
so off. Do not just take the later index -- in `Parametric` S1 mode the spike
counts are set equal to the phd values, so a later index can be a copy of the
right one here and a different quantity under `Full`.

The two vectors do not share a layout. `GetS2` carries the extracted-electron
count and the S2 photon count ahead of the areas, so its indices run two
higher than `GetS1`'s:

    quantity                GetS1   GetS2
    raw area, phd               4       6
    corrected area, phd         5       7

**At commit 541b643 the answer is `s1c_idx = 5`, `s2c_idx = 7`.** Both land
at 0.9992 of \(g N\), which is the signature to look for. Run the probe
anyway on a different checkout: getting this wrong would not crash anything,
it would quietly scan the wrong quantity, which is why there is no default.

The probe also reports the extracted-electron count as element 0 of `GetS2`,
and that is worth reading: at commit 541b643 it gives 868.4 for
\(\langle N_e\rangle=1240\), so the extraction efficiency is 0.70 and the
single-electron size \(g_2/\varepsilon = 49.3\) phd. `analysis/smearing.py`
takes both from here rather than assuming them.

### 6.2 Run it

    ./run_scan24_smeared.sh 8 ./lz_scan_s1s2 <s1c_idx> <s2c_idx>

This writes `scan/scan24_smeared.txt` and never touches
`scan/scan24_results.txt`, so the two can be compared.

Runtime is close to the quanta-level scan despite the extra work. `GetS1` and
`GetS2` are called only for events that lie within eight standard deviations
of a cut boundary; at a leakage of order \(10^{-7}\) almost none do, and the
`n_slow` column reports how many took the slow path. Expect well under a per
cent.

The quanta-level columns in the new file will agree with
`scan24_results.txt` statistically but not bit-for-bit: `GetS1` and `GetS2`
draw from the same `RandomGen` singleton, so the sequence of `GetQuanta`
draws differs. That does not matter, because the comparison that the file is
for is paired *within* it -- the same events classified both ways -- and that
one is exact.

### 6.3 Read it

    cd ../analysis && python read_smeared_scan.py

prints the per-energy ratios, the window integral both ways, and the pooled
paired ratio with its error. Compare it against `smearing.py`, which
estimates the same quantities analytically from the quanta-level scan: the
analytic estimate is that the leakage rises by a factor of about 1.2 and the
between-family ratios are compressed by 9 to 15 per cent, with the
thresholds moving by under 0.01 in kappa. If the measured factor differs,
the measured one is what belongs in the text, and `smearing.py`'s knob scan
indicates which assumption to revisit.

## Command-line interfaces

    build.sh <NEST source dir> [build dir]

    lz_scan <E_keV> <N_trials> <seed>

Writes one line per invocation with the columns documented in
`../scan/COLUMNS.md`. The seed argument is mandatory: because `RandomGen` is a
singleton with a fixed initial state, parallel processes without distinct
seeds return byte-identical output.

    lz_tail_test [N_trials] [skewness_override]

With no arguments it prints the response parameters analytically. With
N_trials it samples. The second argument overrides `SkewnessER` and was used
for the sensitivity scan of Sec. 4; omit it for the nominal model.

    lz_scan_s1s2 probe <E_keV> [N]
    lz_scan_s1s2 <E_keV> <N_trials> <seed> <s1c_idx> <s2c_idx>

Probe mode prints the layout of the vectors returned by `GetS1` and `GetS2`
so the two indices can be read off; scan mode requires them, and has no
default, because a wrong index would not fail -- it would scan the wrong
quantity. The output carries the quanta-level and the S1c/S2c counters for
the same events, plus `n_slow`, the number that needed the full detector
treatment.

## If something goes wrong

**cmake cannot find GSL.** NEST needs it. `apt install libgsl-dev` or
`brew install gsl`, then re-run cmake.

**`fatal error: NEST.hh: No such file or directory`.** Two causes, and the
file named in the error tells you which. If it is one of this archive's
drivers and the command was `make` in NEST's build directory, then driver
targets have been added to NEST's `CMakeLists.txt` and `make` is failing on
them before it builds the library; `ls build/libNESTCore.a` confirms it is
missing. Undo the edits as section 2 describes. If instead you invoked the
compiler yourself, the include directories are simply not on the command
line: use `build.sh`, or name them as in section 3.

**Undefined reference to `LZ_Detector_2024` or `NESTcalc::`.** The library
was not found or not linked. At this commit it is `build/libNESTCore.a` and
the CMake target is `NESTCore`, not `NEST`. `build.sh` prints the path it
used; check that it is the library you built, and that `make` finished.

**Undefined references to `gsl_*`.** Append `$(gsl-config --libs)` to the
compile line, or run `EXTRA_LIBS="$(gsl-config --libs)" ./build.sh ...`.

**`fatal error: gcem.hpp: No such file or directory`**, or any other header
you have never heard of. NEST's headers include header-only dependencies that
its CMake fetches into the build tree rather than the source tree, so an
include list built from the NEST headers alone is incomplete. `build.sh`
avoids this by reading the include flags out of CMake's own `flags.make`; if
it reports that it fell back to discovery and still fails, configure NEST
first -- `cmake ..` alone is enough to unpack the dependencies.

**The yields look nothing like the paper's.** Check that your driver calls
`SetEFT2024Configuration()` after instantiating the detector. Without it the
detector carries the WS2024 defaults, which differ.

**Different numbers from the archived scan.** Check the commit first
(`git rev-parse HEAD` should give 541b643), then the seeds.

**`no matching function for call to NESTcalc::GetS1` / `GetS2` /
`SetDriftVelocity`.** These signatures changed across NEST versions. At
commit 541b643 they are, from `include/NEST/NEST.hh`:

    line 335  const vector<double>& GetS1(const QuantaResult&,
                  truthX, truthY, truthZ, smearX, smearY, smearZ,
                  driftVelocity, dV_mid, INTERACTION_TYPE, evtNum,
                  dfield, energy, S1CalculationMode, useTiming,
                  wf_time, wf_amp)                              -- 17 arguments
    line 351  const vector<double>& GetS2(Ne,
                  truthX, truthY, truthZ, smearX, smearY, smearZ,
                  dt, driftVelocity, evtNum, dfield,
                  S2CalculationMode, useTiming,
                  wf_time, wf_amp, g2_params)                   -- 16 arguments
    line 368  double SetDriftVelocity(T, density, field, pressure)

Neither getter takes an `outputTiming` flag, both return a const reference
rather than a value, and `SetDriftVelocity` takes the pressure. If your
checkout differs, the compiler prints the candidate it found, which is
enough to adjust the call.
