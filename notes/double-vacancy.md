# The double-vacancy populations

The published treatment is Sec. S4 of the supplementary material. This note
is the fuller working behind it, including the bounds and the placement
arithmetic the Supplement does not carry. The manuscript keeps the 124Xe rate
cross-check, as a validation of the reconstruction, and a short paragraph
motivating the request for the double-vacancy distributions.

The question is whether the
calibrated ordinary ER response describes the 124Xe and 125I populations at
all. The Letter states that these are modelled following PRD 112 012024 and
that "there is a systematic uncertainty in projecting into the tail of the
distribution that we do not include in the statistical inference."

## 124Xe two-neutrino double electron capture

T_1/2 = 1.09e22 yr, KK fraction 64.8% (J. Phys. G 52, 015103). The exposure
contains 787 2nuDEC decays, of which 510 are KK captures at 64.3 keV.

Rate cross-check: that paper states 2.3% of the signal lies below the region
of interest of that analysis; that fraction gives 18 events against the
21.0 +/- 6.3 of Table I of the Letter.

The KK mode is not strongly charge-suppressed. In 88% of K-captures a 28 keV
X-ray carries energy O(100) um from the primary site, so the deposit comprises
three independent recombination sites. Only the both-Auger channel,
(1-f_K)^2 ~ 1.6% or ~7 events, is a single dense deposit; with the measured
box-size ratio a_L/a_beta = 0.84 +/- 0.03 in the Thomas-Imel framework this
gives Qy/Qy_beta >= 0.774 and a contribution of 5e-5 counts. Negligible.

Why a_KK >= a_L is the conservative choice:
 - scaling the box size with the electron range would give much larger
   suppression and is not supported by the data: a_L/a_beta = 0.84 corresponds
   to a range ratio of 0.10, so the box size is far less sensitive to track
   structure than the range;
 - the KK double-Auger cascade is dominated by ~23 keV KLL electrons, more
   beta-like than the sub-keV L cascade.

## 125I

The mechanism is the Letter's own: its Fig. 4 caption records that the
features above 25 keVee come predominantly from 125I and 133Xe and "appear
slightly below the beta-decay ER band due to increased recombination from
Auger cascades".

125I decays by electron capture to 125Te, followed by a 35.49 keV nuclear
transition internally converted in 93% of decays. That energy lies just above
the xenon K edge at 34.56 keV, where the mean free path is ~0.14 mm, so it
deposits locally in essentially every decay and produces a second vacancy;
the capture-shell label refers to the first vacancy only.

K, L and M modes deposit 67.3, 40.4 and 36.5 keV, all above the 20.02 keVee
at which a mean-recombination ER meets the S2c ceiling. The 8.9 +/- 2.7 events
of Table I therefore cannot be ordinary electron recoils of those energies:
they are a population whose MEAN charge yield is suppressed.

  shell  Qy/Qy_b  <Ne>   sigma_Ne  (<Ne>-ceil)/sigma  P(Ne<=409)  z_G
  K      0.773    916    194       2.61               8.3e-4      3.1
  L      0.783    562     99       1.54               2.8e-2      1.9
  M      0.785    517     87       1.24               7.7e-2      1.4

These are computed by analysis/double_vacancy.py, against the ceiling
response.py defines. Do not edit the figures here: change the script and
re-run verify.sh, which checks them against the Supplement.

<Ne> = (Qy/Qy_beta) Qy(E) E at each mode's own energy. sigma_Ne is the full
width including the binomial term, with omega from Supplement Eq. 4 at that
N_q -- not a constant: it falls from 0.0456 at N_q = 4973 to 0.0369 at 2697.
z_G = Phi^-1(1-P) exceeds (<Ne>-ceiling)/sigma because the distribution is
positively skewed, so its low-side tail is thinner than Gaussian.

Implied recombination fractions, 1 - <Ne>/N_i with N_i at each mode's own
energy: 0.783, 0.778, 0.774. These are NOT the Qy/Qy_beta ratios, which they
resemble numerically and order in the opposite sense; each exceeds the nominal
r-bar at that energy by the double-vacancy enhancement.

### What this does and does not show

The L and M modes reach the WS RoI through 1.9 and 1.4 Gaussian-equivalent
sigma, probabilities of order 1e-2 to 1e-1. They are NOT a probe of the
P ~ 1e-7 regime. What places them in the search region is the shift of the
mean, not a rare fluctuation.

Only the K mode, at 3.1 sigma, approaches the displacement at issue, and only
its S1c interval brackets the event: S1c > 504 phd from Ne < 409, and
S1c < 549 phd from Ne >= 0, i.e. S1c <= N_q g1.

Mode placement: L and M, at N_q = 2988 and 2697, land at S1c ~ 285-330 and
253-298 phd -- inside the WS RoI but well outside the comparison region.
Weighted by capture fraction they supply ~92% of the 125I population in the
WS RoI and K the remaining 8%; not because K is absent from its S1c range,
but because reaching Ne < 409 from a mean of 916 demands a much larger
downward fluctuation than L and M require from 562 and 517. Implied total
~1e3 decays.

### Non-identifiability

Assigning the L mode the single-vacancy L-capture ratio measured in PRD 112
012024 (0.877) rather than the double-vacancy treatment moves the fraction
above S1c = 504 phd from 8% to 82%, and the inferred total 125I decay count
by an order of magnitude. Neither assignment is a measurement of the 125I
response. This is the regime in which that paper states the Thomas-Imel
framework "cannot predict the size of the effect a priori".
