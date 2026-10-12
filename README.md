# Reproduction archive

Code and derived data for "The Electron-Recoil Response Tail and the LZ
High-Energy Event Background".

Everything here is derived from public material: the LZ preprint and its
supplement, the public NEST repository, and published figures of
Phys. Rev. C 102, 065501. No collaboration-internal data is used.

The extracted inputs are checked against the sources. Eq. 4 of the Supplement
is reproduced verbatim and carries no normalisation factor, which is what the
factor of 1.475 rests on; Table S4's eight parameters and Table S3's ten match
response.py term for term, including the m6 = 0 that the released
configuration contradicts; the WS RoI bounds and the event's S1c and S2c
match; and Fig. 5's bottom panel is 500 < S1c < 600 phd, so the 0.0106 the
thresholds are expressed against is the total background of the same region
the paper computes in.

## Dependencies

NEST, public repository, **commit 541b643**. Not redistributed here; clone it
separately, build it, then run nest-drivers/build.sh against the checkout --
see nest-drivers/HOWTO.md for the full procedure. Do not add the drivers as
CMake targets linking the NEST library and expect the headers to be found:
NEST does not export its include directories in this commit, and the result
is "NEST.hh: No such file or directory". build.sh discovers them.

The EFT2024 configuration selects a branch of `RecombOmegaER` that tagged
releases before this commit do not contain, so a tagged release will not
reproduce these numbers. See nest-drivers/BUILD.md.

Python 3.10+ with numpy, scipy and matplotlib for everything in analysis/.
The digitisation scripts in spectra/ additionally need Pillow, and are only
needed to re-extract the CSVs from the source figures.

## Layout

    nest-drivers/   C++ drivers against the LZ detector configuration
                    (build.sh compiles them; HOWTO.md is the step-by-step
                    build and run guide)
    scan/           the 24-point leakage scan: the one expensive artefact
    analysis/       verification, family construction, thresholds, figures
    spectra/        digitised beta spectra and the digitisation code
    notes/          material not in the paper
    paper/          manuscript source, figures, LZ.bib and JHEP.bst

## Checking it

    ./verify.sh              # about five minutes
    ./verify.sh --figures    # and regenerate the figures

Runs every consistency check and prints a summary; exits non-zero if anything
disagrees. Nothing here needs NEST -- the scan output is committed. Nineteen
checks (sixteen plus three figures), of which four can fail rather than merely
compute:

  check_manuscript.py reads paper/main.tex and paper/supplementary.tex, pulls
  every table and the headline scalars out of them, recomputes each from the
  code, and tests whether the printed number is the correctly rounded form of
  the computed one. 129 values: the family table including the event-density
  column, the matched alpha values, the response at the event, the threshold
  ladder, the empty-interval table, the measured smearing factor, the
  calibration-size table of Sec. 5 with the two entries the summary repeats,
  the cross-family comparison table of Sec. 3, the window-weighted ratios, the
  Table I accounting, the appendix's OldW13eV figures, and the double-vacancy
  table of Sec. S4. A row edited by hand and not regenerated is caught here,
  which is what it exists for. It is the only check that reads the manuscript
  rather than keeping its own copy of the numbers, which is why it is the
  authority when it and leakage.py disagree.

  leakage.py recomputes the 61 values quoted in the prose of Sec. 4 and
  prints each beside its counterpart in the text.

  discriminate.py --selftest evaluates the Sec. 5 divergence at kappa = 2,
  where it must vanish because that member IS the implemented response. What
  it returns instead is the numerical floor of that calculation, and the check
  fails if the floor is not far above everything the calculation reports.

  window_ratio.py rebuilds every family member at every energy across the
  parent window, which takes about two minutes and is the slowest check.

verify.sh begins with a preflight that refuses to run on a partial or
inconsistent tree. It requires every script to import response.py, and it
refuses any private numeric definition of K, the ceiling, NLEAK, DNDE, TOTBKG
or WININT -- including inside a script that does import response.py, since an
import does not prevent a stale constant sitting beside it.

## Reproducing the paper

    # 0. build the drivers against your NEST checkout
    cd nest-drivers && ./build.sh /path/to/NEST

    # 1. the scan (~2 h on 8 cores) -- regenerates scan/scan24_results.txt
    bash run_scan24.sh 8 ./lz_scan

    # 1b. optional: the same scan carried through to S1c and S2c, so that the
    #     detection fluctuations are included. Find the two vector indices
    #     first -- see HOWTO.md section 8 -- then
    ./lz_scan_s1s2 probe 70
    bash run_scan24_smeared.sh 8 ./lz_scan_s1s2 <s1c> <s2c>

    # 2. the checks: every number quoted in the text, recomputed
    cd analysis
    python check_manuscript.py  # every table in the paper, against the code
    python response.py          # the one response convention, printed
    python leakage.py           # the Sec. 4 prose values; non-zero exit on any
                                #   mismatch, so it works as a regression test
    python fit_K.py             # -> K, its uncertainty, deviance, residuals
    python printed_omega.py     # printed vs implemented Eq. 4, and what the
                                #   printed form would do to the background
    python gap.py               # the empty-interval test of Sec. 4
    python wvalue.py            # OldW13eV acts through gamma_ER alone
    python family_table.py      # regenerates the family table of Sec. 3
    python families.py          # the three families, matched-deviation table
    python window_ratio.py      # window-weighted vs single-energy ratios (~2 min)
    python double_vacancy.py    # the Supplement's K/L/M mode table
    python smearing.py          # what detector smearing does to it (~3 min)
    python read_smeared_scan.py # ... and what the S1c/S2c scan measures, if run
    python discriminate.py --selftest   # the numerical floor of the next one
    python discriminate.py      # calibration events needed, Sec. 5
    python discriminate.py --grid       # ... and its grid convergence

    # 3. the figures (a minute or two each)
    python figfam.py            # -> lz_tail_family.pdf
    python fig_ratio.py         # -> lz_ratio_vs_deviation.pdf
    python fig_leakage_scan.py  # -> lz_leakage_scan.pdf         (reads ../scan/)

Runtimes on 8 cores: leakage.py and fig_ratio.py about a minute each,
window_ratio.py and smearing.py two to three, response.py, fit_K.py, gap.py,
printed_omega.py and double_vacancy.py a few seconds, figfam.py several
minutes (it solves a moment-matching problem on a fine grid for every family
member). Only the scan in step 1 requires hours.

## Expected results

    scan             1.084e10 trials, 24 energies, 935 leakage counts
    ceiling          P(Ne <= 409), the CDF at the half-integer edge 409.5
    normalisation    K = 0.937 +/- 0.031, Poisson deviance/dof = 0.91 (p = 0.59)
    leakage          N_leak^ER = 1.845e-4 counts, 1.74% of the quoted total
    printed Eq. 4    omega 0.0682 against 0.0462, leakage ratio 3.2e3,
                     N_leak^ER = 0.60 counts, 56x the quoted total
    reference        gamma_ER = 1.8923, sigma_Ne = 204.12, P = 1.318e-7
    families         ratio 0.0032 (kappa 2.3) to 521 (kappa 1.4), 5.2 orders
    window weighting integrated over the parent window: 12.8, 36.9, 96.3, 521
                     against single-energy 13, 37, 97, 521 -- within 2%
    thresholds       0.14% -> 10% of quoted background; 0.47% -> the whole of it
    empty interval   269 < Ne < 578: 0.02 counts at kappa = 2, 0.69 at kappa = 1.4
    width needed     +6.2% of sigma_Ne reaches a tenth of the quoted total,
                     +16.6% reaches the whole of it
    OldW13eV         acts through gamma_ER alone, 1.89 -> 2.29, leakage / 10.3
    event density    R_ev = 49 (kappa 1.8) to 7.0e3 (kappa 1.4), ~10x steeper
                     than the integrated ratio R
    Table I          the six rows cover 1646.9 of the 1706.8 expected counts;
                     the residual is 59.9, of which 59.7 is ER and 0.2 NR
    smearing         MEASURED N_leak x1.172 (counts 1.18 +- 0.10); modelled
                     ratio compression 7-13%; thresholds move < 0.01 in kappa
    S1c/S2c indices  s1c_idx = 5, s2c_idx = 7 at commit 541b643
    smeared scan     its quanta-level window integral reproduces the archived
                     1.497e-6 keV to four figures; smeared 1.755e-6 keV
    calibration size N for 3 sigma with the three moments fitted: 2.9e4 at
                     kappa 1.87, 3.2e3 at 1.66, 5.4e2 at 1.39; x2.8 for 5
                     sigma; floating the moments costs 1.1-1.2, the charge
                     resolution 1.07; bulk-only fraction 0.38 / 0.78 / 0.91
    double vacancy   K/L/M modes at 2.61, 1.54, 1.24 sigma above the ceiling,
                     P 8.3e-4 / 2.8e-2 / 7.7e-2, Gaussian-equivalent 3.1/1.9/1.4
    numerical floor  kappa = 2 against the response it should equal gives
                     D = 2.1e-10, i.e. N = 2.1e10, 2e4 times the ceiling

The scan is deterministic: run_scan24.sh fixes seeds 201-224 and NEST's
generator is a singleton, so a rerun reproduces these numbers exactly rather
than only statistically. run_scan24_smeared.sh uses the same seeds but will
not reproduce them bit-for-bit, because GetS1 and GetS2 draw from the same
singleton and shift the GetQuanta sequence; its own quanta-level and S1c/S2c
columns are a paired comparison on identical events, and that one is exact.

## The response, and the one place it is defined

response.py holds the single definition of the executable response at the
event-conditioned point. Every script that needs it imports it; none
redefines it. The convention is the implementation's, evaluated rather than
rounded: omega = 0.046218 from RecombOmegaER, gamma_ER = 1.8923 from the ER
skewness parametrisation at E = 69.8366 keVee rather than the 1.9 quoted in
the text, and N_q = 5160.85 from the event's S1c and S2c. It also publishes
the constants that scale results -- K, NLEAK, DNDE, TOTBKG and WININT --
which fit_K.py, leakage.py and read_smeared_scan.py re-derive and which fail
if they no longer match, so a constant cannot drift from the measurement it
came from.

The ceiling needs care. The scan cuts on an INTEGER electron count,
`(double)q.electrons < pow(10,4.15)/34.5`, which admits every Ne up to and
including 409: what it measures is P(Ne <= 409). A continuous density
reproduces that at the half-integer edge, so the analytic representation is
evaluated at 409.5, not at 10^4.15/g2 = 409.43 and not at 409. The choice
matters far less than making it once. K is a normalisation fit at the scan
energies and absorbs most of it: evaluated consistently at 409.0, 409.43 or
409.5 the event-conditioned K*P comes to 1.3204, 1.3187 and 1.3184e-7, a
spread of 0.15%. Fitting at one ceiling and applying at another does not
cancel and costs 1.75% on every absolute number, which is why verify.sh
refuses a private definition of it.

## What each script establishes

printed_omega.py evaluates Eq. 4 of the Supplement twice -- once with the
component-wise mode normalisation that RecombOmegaER applies, and once with
the printed expression taken literally -- and folds both through the parent
window of Sec. 4. It is the source of the two omega values, the width and
leakage ratios, and the count on which the paper's argument rests that the
implementation value, not the printed one, is what produced the published
background expectation. NEST has no printed-form branch, so no scan is
available for that case and the window integral there is necessarily
analytic; the script prints the analytic integral for the implementation
alongside the measured 24-point value, so that the size of the omitted Nq
(Fano) smearing can be read off directly. It is 0.2%, which bounds the error
the analytic treatment of the printed form introduces.

fit_K.py regenerates the normalisation from scan/scan24_results.txt by Poisson
maximum likelihood, using the moments NEST returned at each scan point. It
prints K, its uncertainty, the deviance and the per-point residuals, and
asserts that the K it fits is the K response.py publishes.

leakage.py builds dN/dE in the parent window from the digitised spectra and
Table I of the Letter, folds it through the response to give N_leak, and
recomputes the sensitivities and the three thresholds. It also holds all
thirteen rows of Table I and asserts that the components Sec. 4 folds in, the
components it names as excluded, and the nuclear-recoil entries close exactly
on the expected total -- so the prose cannot drift from the table. Table I
prints two columns and they must not be mixed: the 1713 +/- 39 at its foot is
the FIT total and carries the fitted signal, while the expected column, which
is what a reconstruction is compared against, sums to 1706.8.

family_table.py regenerates the family table of Sec. 3. Its ratio column is
normalised on the grid-built kappa = 2 member rather than on the exact
skew-normal, so the +0.097% tail error of the grid cancels between numerator
and denominator instead of biasing every row.

families.py builds the three families of Sec. 3 -- generalised normal, skew-t
and Johnson S_U, the last via the moment-matching fit in johnson.py, each
matched to the executable model's first three moments --
and prints the matched-deviation table. fig_ratio.py plots these lists and
check_manuscript.py tests the table against them, so the figure and the table
cannot disagree. One caveat the script states: the generalised-normal and
skew-t members are grid-built and carry about +0.1% of tail error, the Johnson
S_U members are analytic and carry none, so the three are mutually consistent
only to that order. That is far below the factor of about two the section
claims, and it is the reason the entries are not pushed past three significant
figures.

window_ratio.py integrates the family ratios over the parent window instead of
evaluating them at the event-conditioned energy, which is the comparison
Sec. 3 makes when it states what the single-energy approximation costs. It is
2%. The construction is the point: the target moments move with energy -- the
mean and width through the yield and recombination model, the skewness through
gamma_ER(E) -- so each member is rebuilt at every energy rather than
translated across the window. The integration limits move too, and over most
of the window it is the S1c window rather than the S2c ceiling that binds.

gap.py computes the expected electron-recoil counts in the interval of Fig. 4
of the Letter observed to be empty, for each member of the family, and
tabulates the sensitivity to the one input read off that figure -- the lower
edge of the observed population at log10 S2c = 4.30.

wvalue.py checks the appendix claim that the OldW13eV convention acts on the
tail through the skewness alone. The yields and the width are evaluated from
N_q, which the flag does not touch, so recombProb, omega and sigma_Ne are
identical between the two settings; W enters only where quanta are converted
to an energy to evaluate the ER skewness parametrisation, moving gamma_ER from
1.89 to 2.29. At fixed mean and width that shift gives 10.0 linearised from
the quoted derivative and 12.5 with the curvature, against 10.3 measured in
NEST. The script also reproduces the factor of 57 that Sec. 2 quotes for
gamma_ER in [1.55, 2.25], which is the same calculation over a wider range.

smearing.py replaces the hard quanta-level cuts of the scan with the
probability that a given (N_ph, N_e) satisfies them, so that the photon
detection and electron extraction fluctuations the scan omits can be sized.
At the event's N_ph the S1c fluctuation is 24 phd, which is 216 photons or a
quarter of the window width; the S2c ceiling blurs by 7 electrons in 409. The
leakage of every family rises about 32% against a measured 17.2%, so the
estimate is 13% high and 1.4 sigma from the measurement. Calibrated to the
measurement, the between-family ratios are compressed by 7% at kappa = 1.8
rising to 13% at kappa = 1.4, and the two effects very nearly cancel in the
thresholds, each kappa moving by under 0.01. The compression is the one
quantity only this script can supply, the scan existing only for the response
the configuration implements. Two of the three detector parameters are
measured rather than assumed: lz_scan_s1s2's probe mode reports the
extracted-electron count, which gives an extraction efficiency of 0.70 and a
single-electron size of 49.3 phd, and separating the observed spread of the
corrected S2 area from the spread of Ne leaves a single-electron resolution of
0.20. The single-photoelectron resolution remains a knob and the last block of
the script scans all three; the dominant counting term N_ph g1 (1-g1) has no
free parameter.

lz_scan_s1s2.cpp is the scan driver carried through GetS1 and GetS2, with the
cuts applied in phd. It is the measured answer to the question smearing.py
estimates. Two notes. The layout of the vectors GetS1 and GetS2 return varies
between NEST versions, so the driver takes the two indices as required
arguments and ships a probe mode that prints them with their means, spreads
and expected values -- a wrong index would not crash, it would quietly scan
the wrong quantity. And GetS1/GetS2 are called only for events within eight
standard deviations of a cut boundary, which at a leakage of 1e-7 is a tiny
fraction, so the run costs about what the quanta-level scan costs; the n_slow
column reports the fraction. read_smeared_scan.py reduces the output and
compares it with the analytic estimate.

discriminate.py answers the objection to Sec. 5 that the families differ in
the bulk as well as the tail, so that a shape fit to an ordinary beta
calibration might serve instead of measuring the tail to 1e-5. It takes the
heavier member as truth, fits the implemented response to it with the mean,
width and skewness free, and reads the separation off the Asimov form
Z^2 = 2 N D(q||p_hat). Three things come out of it. The nuisances are nearly
free: fitting the three moments rather than fixing them raises N by only
10-20%, because moments the two distributions share exactly cannot absorb a
difference of shape, and what the fit does instead is move the skewness by
about one standard error of the skewness at that same N. The detector is
nearly irrelevant: folding in the measured charge resolution raises every
entry by 7%. And the premise holds only in part -- restricting the divergence
to the central 99.98% of the truth leaves 38% of it at kappa = 1.87 and 78% at
1.66, so for the member nearest the implemented response the shape fit reaches
into the same tail by a less direct route.

Two guards, because a divergence calculation on a grid can report a small
number for numerical reasons rather than physical ones. --selftest evaluates
kappa = 2, which IS the implemented response, so the divergence must vanish;
what it returns instead is the numerical floor, N = 2.1e10, four orders above
the 1e6 ceiling the main run asserts nothing crosses. It also checks that each
constructed truth carries its three target moments to 1e-4. --grid reruns at
4, 2 and 1 e- steps; the answers agree to four figures, so the default 2 e-
step is not doing any of the work. Both are in verify.sh.

double_vacancy.py reproduces the K/L/M mode table of Sec. S4: how far the
charge-suppressed 125-I modes reach into the WS RoI. Its point is the gap
between its last two columns. The L and M modes sit 1.5 and 1.2 sigma_Ne above
the ceiling but 1.9 and 1.4 GAUSSIAN equivalents, because the response is
positively skewed, and reading the first as the second is the error the table
forestalls.

band.py, band_check.py, trace_band.py and band_published.csv belong to an
attempted reconstruction of Fig. 2 of the Letter which did not succeed;
notes/band-reconstruction.md records what was verified, what was not, and why
no claim in the paper rests on it.

## Licences

Code MIT. Tabulated data CC BY 4.0. The digitised spectra are digitisations
of published figures, not the original authors' data -- see spectra/PROVENANCE.md.
