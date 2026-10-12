"""The electron-recoil contribution the *printed* form of Eq. 4 of the
Supplement would give, against the one the released implementation gives.

Eq. 4 of the Supplement specifies the recombination-fluctuation width as a sum
of two skew-Gaussians in x = log10 Nq, with the amplitudes, centres, widths and
skew parameters of Table S4.  The released code (RecombOmegaER) evaluates the
same sum with each component divided by its own value at its own mode.  The two
therefore disagree at every Nq, by a constant factor; this script quantifies
what that disagreement does to the background expectation.

The quantity of interest is

    N_leak  =  dN/dE  x  INT P(E) dE   over the parent window,

with dN/dE = 123.2 per keVee from leakage.py and P(E) the per-decay probability
of landing in the comparison region.  The paper takes the window integral for
the implementation from the 24-point NEST scan (WININT below), because the
analytic skew-normal representation omits Nq (Fano) smearing.  That correction
is small: with the ceiling applied consistently (see response.py) the analytic
integral sits 0.2% from the measured one, not the 2% reported before the
ceiling was made consistent -- most of what was attributed to Fano smearing
was in fact a fit-versus-apply mismatch in the ceiling.  NEST contains no printed-form branch, so no scan is available for that
case and the integral must be analytic.  Two numbers are therefore printed:

  (a) the bare analytic integral, which uses the analytic representation for
      both forms and so is self-consistent but understates both;
  (b) the same rescaled by WININT / I_analytic(implementation), which assumes
      the Fano correction is common to the two forms.

The two now differ by 0.2%, so the choice is immaterial; it mattered more when
the correction was thought to be 2%.  The argument of the paper does not
depend on which is used, and the text should quote (a) as the conservative
choice.  Note that (b) is an upper bound on the correction: the Fano term
matters less at larger omega, since the ceiling then sits at a smaller
displacement in units of the total width.

    cd analysis && python printed_omega.py
"""
import numpy as np
from math import erf, sqrt, exp, log10, pi
from scipy.stats import skewnorm
from scipy.integrate import quad

# ---- the same constants leakage.py uses -------------------------------------
import sys; sys.path.insert(0, '.')
import response as _R
g1, g2, W, RATIO, C, K = _R.g1, _R.g2, _R.W, _R.RATIO, _R.C, _R.K
WININT = _R.WININT   # measured, 24-point scan; response.py publishes it
DNDE   = _R.DNDE         # parent rate in the window [per keVee]; leakage.py derives it
TOTBKG = _R.TOTBKG       # LZ total for the comparison region
WLO, WHI = 61.3, 79.1    # parent window [keVee]
inv, s2 = 1/sqrt(2*pi), sqrt(2)

M = dict(m1=12.31, m2=84.91, m3=0.5707, m4=2.804, m5=34.04,
         m7=82.57, m8=4.557, m9=0.2207, m10=0.1272)
Qy = lambda E: C*(M['m1'] + (M['m2']-M['m1'])/(1+(E/M['m3'])**M['m4'])**M['m9']
                + M['m5'] + (0.0-M['m5'])/(1+(E/M['m7'])**M['m8'])**M['m10'])
S4 = [(0.03174, 2.588, 0.3627, -0.3583), (0.06211, 5.200, 1.359, -2.599)]

def omega(Nq, printed=False):
    """Eq. 4 of the Supplement.  printed=False applies the component-wise mode
    normalisation carried in RecombOmegaER; printed=True does not."""
    x, t = log10(Nq), 0.
    for A, c, w, sk in S4:
        n = 1.0
        if not printed:
            mo = c + 2*inv*sk*w/sqrt(1+sk*sk)          # mode of component i
            n  = 1/(exp(-.5*(mo-c)**2/w**2)*(1+erf(sk*(mo-c)/(w*s2))))
        t += n*A*exp(-.5*(x-c)**2/w**2)*(1+erf(sk*(x-c)/(w*s2)))
    return t

F = 96.5
a0,c0,c1,E0,E1,E2,E3,F0,F1 = 1.39,4.,22.1,7.7,54.,26.7,6.4,225.,71.
alpha = lambda E: (1/(1+exp((E-E2)/E3))*(a0+c0*exp(-F/F0)*(1-exp(-E/E0)))
                 + 1/(1+exp(-(E-E2)/E3))*c1*exp(-E/E1)*exp(-sqrt(F)/sqrt(F1)))

NPH_MIN, NPH_MAX = 500/g1, 600/g1
NE_CEIL = _R.NE_CEIL

def state(E, printed=False):
    Nq = E/W; Ni = Nq/(1+RATIO); mu = Qy(E)*E; r = 1-mu/Ni
    sd = sqrt(r*(1-r)*Ni + (omega(Nq, printed)*Ni)**2)
    return Nq, Ni, mu, sd, alpha(E)

def P_window(E, printed=False):
    """per-decay probability of falling in the comparison region"""
    Nq, Ni, mu, sd, a = state(E, printed)
    d = a/sqrt(1+a*a); w = sd/sqrt(1-(2/pi)*d*d); xi = mu - w*d*sqrt(2/pi)
    hi = min(NE_CEIL, Nq-NPH_MIN); lo = max(0.0, Nq-NPH_MAX)
    if hi <= lo: return 0.0
    return K*(skewnorm.cdf(hi,a,loc=xi,scale=w) - skewnorm.cdf(lo,a,loc=xi,scale=w))

# ---- at the event-conditioned energy ----------------------------------------
NQ_EV = 540.1/g1 + 9268/g2
E_EV  = NQ_EV*W
print('\n--- at the event-conditioned Nq = %.0f (E = %.2f keVee) ---' % (NQ_EV, E_EV))
oi, op = omega(NQ_EV, False), omega(NQ_EV, True)
si = state(E_EV, False)[3]; sp = state(E_EV, True)[3]
Pi, Pp = P_window(E_EV, False), P_window(E_EV, True)
print('  omega                implementation %7.4f    printed %7.4f    ratio %6.3f'
      % (oi, op, op/oi))
print('  sigma_Ne [e-]        implementation %7.1f    printed %7.1f    ratio %6.3f'
      % (si, sp, sp/si))
print('  P(comparison region) implementation %7.3g    printed %7.3g    ratio %6.0f'
      % (Pi, Pp, Pp/Pi))
print('  (both evaluated at the same Nq; the omega ratio is Nq-independent,')
print('   the component normalisations being constants.  At Nq = 5243:')
print('   implementation %.4f, printed %.4f, ratio %.3f)'
      % (omega(5243.,False), omega(5243.,True), omega(5243.,True)/omega(5243.,False)))

# ---- folded through the parent window ---------------------------------------
Ii = quad(lambda E: P_window(E, False), WLO, WHI, limit=400)[0]
Ip = quad(lambda E: P_window(E, True),  WLO, WHI, limit=400)[0]
print('\n--- window integral INT P(E) dE over %.1f--%.1f keVee [keV] ---' % (WLO, WHI))
print('  analytic, implementation  %10.4g' % Ii)
print('  analytic, printed         %10.4g      ratio %6.0f' % (Ip, Ip/Ii))
print('  measured (24-point scan)  %10.4g      Fano correction x%.3f'
      % (WININT, WININT/Ii))

print('\n--- electron-recoil contribution in the comparison region [counts] ---')
nom = WININT*DNDE
print('  implementation, as quoted in the paper        %10.4g' % nom)
for lbl, val in (('(a) printed, bare analytic integral', DNDE*Ip),
                 ('(b) printed, Fano-corrected as above', DNDE*Ip*WININT/Ii)):
    print('  %-44s %10.4g   = %5.1f x the quoted total of %.4f'
          % (lbl, val, val/TOTBKG, TOTBKG))
print('\n  The printed form overshoots the Letter\'s own total background for the')
print('  region by one to two orders of magnitude on either construction, so it')
print('  cannot be the specification that produced that total.  The paper quotes')
print('  (a).\n')
