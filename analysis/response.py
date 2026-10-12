"""The one definition of the executable response at the event-conditioned point.

Every script that needs the ER response imports it from here, so that the
family construction, the figures and the leakage integrals cannot drift apart.
The convention is the released implementation's, evaluated rather than rounded:

  * omega from RecombOmegaER, with the component-wise mode normalisation
    (Eq. 4 of the Supplement divided by each component's value at its mode);
  * gamma_ER from the energy-dependent ER-skewness parametrisation evaluated at
    the event-conditioned energy, NOT the rounded value 1.9 quoted in the text;
  * N_q from the event's S1c and S2c, not rounded to 5161.

Exported: W, RATIO, C, g1, g2, K, Qy, omega, alpha, and the reference block
NQ, E, NI, MU, R, SD, GAM, SK, PREF, QREF.
"""
from math import sqrt, pi, exp, erf, log10
from scipy.stats import skewnorm

# --- released configuration --------------------------------------------------
W, RATIO, C, g1, g2 = 13.532e-3, 0.1812, 1.0051, 0.1104, 34.5

# The S2c ceiling of the WS RoI. The scan cuts on an integer electron count,
#   (double)q.electrons < pow(10,4.15)/34.5
# which admits every N_e up to and including 409, so what the scan measures is
# P(N_e <= NE_CUT). A continuous density approximating a count distribution
# reproduces that at the half-integer edge, so the analytic representation is
# evaluated at NE_CUT + 0.5 and NOT at 10^4.15/g2 = 409.4312.
#
# The choice matters far less than using one value everywhere. K is a
# normalisation fit at the scan energies and absorbs most of it: evaluated
# consistently at 409.0, 409.4312 or 409.5 the event-conditioned K*P comes to
# 1.3204, 1.3187 and 1.3184e-7, a spread of 0.15%. Fitting K at one ceiling and
# applying it at another does not cancel and costs 1.75%, which is what earlier
# versions of this archive did; see the note in ../README.md.
NE_CUT  = 409                  # the largest electron count inside the RoI
NE_CEIL = NE_CUT + 0.5         # 409.5, the continuity-corrected edge

# Analytic normalisation, fit to the scan in fit_K.py AT NE_CEIL. fit_K.py
# re-derives it and fails if it no longer matches, so this constant cannot
# drift away from the scan it came from.
K, DK = 0.9374, 0.0307

# The electron-recoil leakage at kappa = 2, folded through the parent window
# and acceptance. leakage.py derives it and fails if this constant no longer
# matches, the same contract K has. It scales the whole N_leak column of the
# family table and the threshold ladder, so it is published here and nowhere
# else. TOTBKG is the Letter's total for the comparison region.
NLEAK, TOTBKG = 1.845e-4, 0.0106
DNDE = 123.2                   # parent rate in the window [per keVee]

# The measured window integral over the parent window, from the 24-point scan.
# Published here rather than repeated as a literal in leakage.py and
# printed_omega.py; read_smeared_scan.py recomputes it from the scan file and
# fails if the two disagree.
WININT = 1.497e-6              # [keV]

M = dict(m1=12.31, m2=84.91, m3=0.5707, m4=2.804, m5=34.04,
         m7=82.57, m8=4.557, m9=0.2207, m10=0.1272)
Qy = lambda E: C*(M['m1'] + (M['m2']-M['m1'])/(1+(E/M['m3'])**M['m4'])**M['m9']
                + M['m5'] + (0.0-M['m5'])/(1+(E/M['m7'])**M['m8'])**M['m10'])

_inv, _s2 = 1/sqrt(2*pi), sqrt(2)
_S4 = [(0.03174, 2.588, 0.3627, -0.3583), (0.06211, 5.200, 1.359, -2.599)]
def omega(Nq, printed=False):
    """Eq. 4 of the Supplement; printed=False applies RecombOmegaER's
    component-wise mode normalisation (see printed_omega.py)."""
    x, t = log10(Nq), 0.
    for A, c, w, sk in _S4:
        n = 1.0
        if not printed:
            mo = c + 2*_inv*sk*w/sqrt(1+sk*sk)
            n  = 1/(exp(-.5*(mo-c)**2/w**2)*(1+erf(sk*(mo-c)/(w*_s2))))
        t += n*A*exp(-.5*(x-c)**2/w**2)*(1+erf(sk*(x-c)/(w*_s2)))
    return t

_F = 96.5
_a0,_c0,_c1,_E0,_E1,_E2,_E3,_F0,_F1 = 1.39,4.,22.1,7.7,54.,26.7,6.4,225.,71.
alpha = lambda E: (1/(1+exp((E-_E2)/_E3))*(_a0+_c0*exp(-_F/_F0)*(1-exp(-E/_E0)))
                 + 1/(1+exp(-(E-_E2)/_E3))*_c1*exp(-E/_E1)*exp(-sqrt(_F)/sqrt(_F1)))

# --- the event-conditioned reference ----------------------------------------
S1C_EV, S2C_EV = 540.1, 9268.
NQ = S1C_EV/g1 + S2C_EV/g2          # 5160.85
E  = NQ*W                           # 69.8366 keVee
NI = NQ/(1+RATIO)
MU = Qy(E)*E                        # mean charge
R  = 1 - MU/NI                      # recombProb
SD = sqrt(R*(1-R)*NI + (omega(NQ)*NI)**2)
GAM = alpha(E)                      # gamma_ER, evaluated not rounded
_d  = GAM/sqrt(1+GAM*GAM)
SK  = (4-pi)/2*(_d*sqrt(2/pi))**3/(1-2*_d*_d/pi)**1.5
_sc = SD/sqrt(1-2*_d*_d/pi)
_lo = MU - _sc*_d*sqrt(2/pi)
PREF = K*skewnorm.cdf(NE_CEIL, GAM, loc=_lo, scale=_sc)
QREF = {p: skewnorm.ppf(p, GAM, loc=_lo, scale=_sc) for p in (0.10, 0.50, 0.90)}

if __name__ == "__main__":
    print('N_q        = %.2f' % NQ)
    print('E          = %.4f keVee' % E)
    print('N_i        = %.1f' % NI)
    print('mu         = %.2f e-' % MU)
    print('recombProb = %.5f' % R)
    print('omega      = %.6f' % omega(NQ))
    print('sigma_Ne   = %.3f e-' % SD)
    print('gamma_ER   = %.4f   (rounded value quoted in the text: 1.9)' % GAM)
    print('K          = %.4f +/- %.4f   (fit at the ceiling, fit_K.py)' % (K, DK))
    print('skewness   = %.6f' % SK)
    print('P(Ne<=%d) = %.4g   (CDF at %.1f)' % (NE_CUT, PREF, NE_CEIL))
    print('quantiles  = ' + ', '.join('%.1f' % QREF[p] for p in (0.10,0.50,0.90)))
    print('displacement (mu-ceiling)/sigma = %.2f' % ((MU-NE_CEIL)/SD))
