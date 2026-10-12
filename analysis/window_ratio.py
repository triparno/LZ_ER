"""The family ratios of Sec. 3, integrated over the parent window.

The table of Sec. 3 quotes tail-probability ratios at the event-conditioned
energy. The alternative is to integrate numerator and denominator over the
parent window with the parent spectrum as weight, which is what folds into the
background expectation. Sec. 3 states how much that changes the ratios, and
that figure was quoted in two places with two different values, neither of
them produced by any code in this archive. This script produces it.

The construction differs from family_table.py in one essential way. At the
event-conditioned energy the family is matched to one set of moments. Over the
window the target moments move with energy -- the mean and width through the
yield and recombination model, and the skewness through gamma_ER(E) -- so the
family member has to be rebuilt at every energy rather than translated. The
integration limits move too, and over most of the window it is the S1c window
and not the S2c ceiling that binds; at the ends of the window the interval
closes entirely.

The parent rate is taken flat across the window, as Sec. 4 does: the window is
61.3--79.1 keVee and every component's endpoint is far above it, so dN/dE is
constant to better than the twenty per cent assigned to it.

    cd analysis && python window_ratio.py
"""
import sys
sys.path.insert(0, '.')
import numpy as np
from math import sqrt, pi
from scipy.optimize import fsolve
from scipy.stats import skewnorm
from scipy.integrate import trapezoid
import response as R

KAPPAS = (1.8, 1.7, 1.6, 1.4)
NPH_MIN, NPH_MAX = 500/R.g1, 600/R.g1
EGRID = np.linspace(61.29, 79.08, 121)      # the kinematic support of Sec. 4

z = np.linspace(-40, 40, 200001); dz = z[1]-z[0]
_kern = {}


def _kernel(k):
    if k not in _kern:
        g = np.exp(-np.abs(z)**k/k); g /= g.sum()*dz
        Gc = np.cumsum(g)*dz; Gc /= Gc[-1]
        _kern[k] = (g, Gc)
    return _kern[k]


def _moments(k, a):
    g, Gc = _kernel(k)
    f = g*np.interp(a*z, z, Gc); n = f.sum()*dz
    m1 = (z*f).sum()*dz/n
    m2 = ((z-m1)**2*f).sum()*dz/n
    m3 = ((z-m1)**3*f).sum()*dz/n
    return m1, sqrt(m2), m3/m2**1.5, np.cumsum(f)*dz/n


# The shape parameter that matches a given skewness depends only on kappa and
# the target, so build a lookup once per kappa and interpolate. Rebuilding it
# at each of the 121 energies would cost a hundredfold for no accuracy.
_tab = {}
def _a_for(k, target):
    if k not in _tab:
        aa = np.linspace(0.05, 6.0, 60)
        sk = [_moments(k, a)[2] for a in aa]
        _tab[k] = (np.array(sk), aa)
    sk, aa = _tab[k]
    a0 = float(np.interp(target, sk, aa))
    return float(fsolve(lambda x: [_moments(k, float(x[0]))[2]-target], [a0])[0])


def state(E):
    """NEST's moments at energy E, as Sec. 2 defines them"""
    Nq = E/R.W
    Ni = Nq/(1+R.RATIO)
    mu = R.Qy(E)*E
    r = 1 - mu/Ni
    sd = sqrt(r*(1-r)*Ni + (R.omega(Nq)*Ni)**2)
    a = R.alpha(E)
    d = a/sqrt(1+a*a)
    return Nq, mu, sd, a, (4-pi)/2*(d*sqrt(2/pi))**3/(1-2*d*d/pi)**1.5


def limits(Nq):
    hi = min(R.NE_CEIL, Nq-NPH_MIN)
    lo = max(0.0, Nq-NPH_MAX)
    return lo, hi


def P_at(E, k=None):
    """the cut probability at energy E, for the family member k (None = NEST)"""
    Nq, mu, sd, a, skw = state(E)
    lo, hi = limits(Nq)
    if hi <= lo:
        return 0.0
    if k is None:
        d = a/sqrt(1+a*a); w = sd/sqrt(1-2*d*d/pi); xi = mu - w*d*sqrt(2/pi)
        return R.K*(skewnorm.cdf(hi, a, loc=xi, scale=w)
                    - skewnorm.cdf(lo, a, loc=xi, scale=w))
    af = _a_for(k, skw)
    m1, s1, _, cdf = _moments(k, af)
    s = sd/s1; xi = mu - m1*s
    return R.K*(np.interp((hi-xi)/s, z, cdf) - np.interp((lo-xi)/s, z, cdf))


if __name__ == "__main__":
    print(__doc__.split('\n\n')[0])
    print('\n  window %.2f--%.2f keVee, %d points, ceiling %.1f e-'
          % (EGRID[0], EGRID[-1], len(EGRID), R.NE_CEIL))

    # the kappa = 2 member built on the same grid, so the grid error cancels
    # exactly as it does in family_table.py
    I0 = trapezoid([P_at(e, 2.0) for e in EGRID], EGRID)
    Inest = trapezoid([P_at(e) for e in EGRID], EGRID)
    print('  grid-built kappa = 2 window integral %.4g keV, against the exact'
          % I0)
    print('  skew-normal %.4g keV: %+.3f%%, which cancels in the ratios below.'
          % (Inest, 100*(I0/Inest-1)))

    print('\n  kappa   window-weighted   single-energy   difference')
    rows = []
    for k in KAPPAS:
        Ik = trapezoid([P_at(e, k) for e in EGRID], EGRID)
        rw = Ik/I0
        # the single-energy value, from the same construction at the event
        Nq, mu, sd, a, skw = state(R.E)
        af = _a_for(k, skw); m1, s1, _, cdf = _moments(k, af)
        s = sd/s1; xi = mu - m1*s
        p = np.interp((R.NE_CEIL-xi)/s, z, cdf)
        af0 = _a_for(2.0, skw); m10, s10, _, cdf0 = _moments(2.0, af0)
        s0 = sd/s10; xi0 = mu - m10*s0
        p0 = np.interp((R.NE_CEIL-xi0)/s0, z, cdf0)
        rs = p/p0
        rows.append((k, rw, rs))
        print('   %4.2f   %13.1f   %13.0f   %+9.1f%%'
              % (k, rw, rs, 100*(rw/rs-1)))

    worst = max(abs(100*(rw/rs-1)) for _, rw, rs in rows)
    print('\n  largest difference over the tabulated members: %.0f%%' % worst)
    print('  LaTeX: window-weighted \\(%s\\) against single-energy \\(%s\\)'
          % (', '.join('%.1f' % r[1] for r in rows),
             ', '.join('%.0f' % r[2] for r in rows)))
    print("""
  Sec. 3 and Sec. 4 both state this figure and stated two different values,
  neither computed. The window-weighted ratio is the larger because the
  heavier members gain most where the cut interval is widest, in the middle of
  the window, while the single-energy value is evaluated at one point near its
  upper end. The difference is small against the factors of 10 to 10^3 the
  argument turns on, which is why the single-energy values are used; but the
  figure quoted for it has to be the computed one.
""")
