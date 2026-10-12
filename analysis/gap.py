"""Expected electron-recoil counts in the interval of Fig. 4 of the Letter that
is observed to be empty, for each member of the generalised-normal family.

The Letter notes "the lack of events between the event of interest and the
bottom of the ER band in Fig. 4". Every point in that figure is search data
plotted across the whole plane, so the remark concerns a genuinely empty
region of the observable space and not a selection boundary. This script
quantifies what the emptiness constrains.

Within 500 < S1c < 600 phd the observed population extends down to
log10 S2c = 4.30 (Ne = 578); the lowest point nearby, at 4.25, sits at
S1c = 495 and so outside the region. The empty interval is therefore between
the event at log10 S2c = 3.98 (Ne = 269) and that edge.

Both bounds are read off the published figure. EDGE below is the one input
that is not computed, and the sensitivity to it is tabulated at the end.

    python gap.py
"""
import numpy as np
from math import sqrt, pi, erf, exp, log10
from scipy.optimize import fsolve
from scipy.integrate import trapezoid
from scipy.stats import poisson

# ---- inputs read off Fig. 4 -------------------------------------------------
EDGE   = 4.30        # log10 S2c of the lowest search event within 500<S1c<600
S2C_EV = 9268.       # S2c of the event of interest [phd]; log10 = 3.967
S1_LO, S1_HI = 500., 600.

# ---- response model, as in Sec. 2 of the paper ------------------------------
import sys; sys.path.insert(0, '.')
import response as _R0
W, RATIO, C, g1, g2, KNORM = (_R0.W, _R0.RATIO, _R0.C, _R0.g1, _R0.g2, _R0.K)
DNDE = _R0.DNDE      # parent rate per keVee; leakage.py derives it
M = dict(m1=12.31,m2=84.91,m3=0.5707,m4=2.804,m5=34.04,
         m7=82.57,m8=4.557,m9=0.2207,m10=0.1272)
Qy = lambda E: C*(M['m1']+(M['m2']-M['m1'])/(1+(E/M['m3'])**M['m4'])**M['m9']
                + M['m5']+(0.0-M['m5'])/(1+(E/M['m7'])**M['m8'])**M['m10'])
S4 = [(0.03174,2.588,0.3627,-0.3583),(0.06211,5.200,1.359,-2.599)]
inv, s2 = 1/sqrt(2*pi), sqrt(2)
def om(Nq):
    x, t = log10(Nq), 0.
    for A,c,w,sk in S4:
        mo = c + 2*inv*sk*w/sqrt(1+sk*sk)
        n  = 1/(exp(-.5*(mo-c)**2/w**2)*(1+erf(sk*(mo-c)/(w*s2))))
        t += n*A*exp(-.5*(x-c)**2/w**2)*(1+erf(sk*(x-c)/(w*s2)))
    return t

z = np.linspace(-40,40,200001); dz = z[1]-z[0]
import sys; sys.path.insert(0,'.')
import response as _R
SK = _R.SK           # the single convention: see response.py
_c = {}
def family(k):
    """generalised normal of shape k, matched in the first three moments"""
    if k in _c: return _c[k]
    g = np.exp(-np.abs(z)**k/k); g /= g.sum()*dz
    G = np.cumsum(g)*dz; G /= G[-1]
    def mo(a):
        f = g*np.interp(a*z,z,G); n = f.sum()*dz
        m1 = (z*f).sum()*dz/n
        m2 = ((z-m1)**2*f).sum()*dz/n
        m3 = ((z-m1)**3*f).sum()*dz/n
        return m1, sqrt(m2), m3/m2**1.5, np.cumsum(f)*dz/n
    a = fsolve(lambda x:[mo(x[0])[2]-SK],[1.9])[0]
    m1, s1, _, cdf = mo(a); _c[k] = (m1,s1,cdf); return _c[k]

NPH_LO, NPH_HI = S1_LO/g1, S1_HI/g1
def counts(k, ne_lo, ne_hi):
    """expected ER counts with ne_lo < Ne < ne_hi and S1_LO < S1c < S1_HI"""
    m1, s1, cdf = family(k)
    Eg = np.linspace(45,170,500); v = []
    for E in Eg:
        Nq = E/W; Ni = Nq/(1+RATIO); mu = Qy(E)*E; r = 1-mu/Ni
        sd = sqrt(max(r*(1-r),1e-9)*Ni+(om(Nq)*Ni)**2)
        s = sd/s1; xi = mu - m1*s
        hi = min(ne_hi, Nq-NPH_LO); lo = max(ne_lo, Nq-NPH_HI)
        v.append(0. if hi<=lo else
                 KNORM*(np.interp((hi-xi)/s,z,cdf)-np.interp((lo-xi)/s,z,cdf)))
    return DNDE*trapezoid(v,Eg)

KS = (2.0,1.8,1.7,1.6,1.5,1.4)
NE_EV, NE_ED = S2C_EV/g2, 10**EDGE/g2
print('Empty interval of Fig. 4: %.3f < log10 S2c < %.2f, i.e. %.0f < Ne < %.0f'
      % (log10(S2C_EV), EDGE, NE_EV, NE_ED))
print('within %.0f < S1c < %.0f phd.\n' % (S1_LO, S1_HI))
print('   kappa   expected   P(0 obs)   exposure x N to reach 3 expected')
for k in KS:
    n = counts(k, NE_EV, NE_ED)
    print('   %.1f     %8.3f     %.3f          %6.1f' % (k, n, poisson.pmf(0,n), 3.0/n))
print('\n  Across the whole band in this S1c slice the same calculation gives')
print('  %.0f electron recoils, which is why the population there is a continuum.'
      % counts(2.0, 0., 4000.))
print('\nSensitivity to the one input read off the figure:\n')
print('   EDGE    Ne      ' + '  '.join('k=%.1f' % k for k in KS))
for ed in (4.20,4.25,4.30,4.35,4.40):
    row = '   %.2f  %5.0f   ' % (ed, 10**ed/g2)
    print(row + '  '.join('%6.3f' % counts(k,NE_EV,10**ed/g2) for k in KS))
print('\n  The ordering is unchanged over this range; no family is excluded at')
print('  any edge value consistent with the figure.')
