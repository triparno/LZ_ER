"""The three distribution families of Sec. 3, constructed once.

Each family reproduces the executable model's mean, standard deviation and
standardised skewness exactly; the remaining shape parameter (kappa for the
generalised normal, nu for the skew-t, delta for the Johnson S_U) is scanned.
For each member this returns the largest fractional shift in the 10/50/90%
quantiles of N_e relative to the executable model, and the tail probability
relative to it.

This lives in its own module so that one construction serves three users:
fig_ratio.py plots these lists, check_manuscript.py tests the Sec. 3 table
against them, and neither builds its own copy. A table produced inside a
figure script is checked by nothing and drifts silently.

One caveat on the normalisation. The generalised-normal and skew-t members are
built on a grid and their tail probabilities carry its error, about +0.1% at
this displacement; the Johnson S_U members use an analytic distribution
function and carry none. The three are therefore mutually consistent only to
that order. It is far below the factor of about two the section claims and
below the precision of the printed table, but it is the reason the entries
cannot be pushed past three significant figures.

    cd analysis && python families.py
"""
import sys
sys.path.insert(0, '.')
import numpy as np
from math import sqrt, pi
from scipy.stats import t as tdist
from scipy.optimize import fsolve
import response as _R

# johnson.py supplies the S_U fit and the reference quantiles; it imports
# response.py, so the convention is the same one.
exec(open('johnson.py').read().split('if __name__')[0])

KAPPAS = (2.4, 2.3, 2.2, 2.1, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.35)
NUS    = (80, 40, 25, 18, 14, 12, 9, 7, 6)
DELTAS = (6, 5, 4.5, 4, 3.5, 3, 2.5, 2.2, 1.8, 1.5, 1.3)
DEVS   = (0.6, 0.8, 1.0, 1.2)        # the matched deviations of the Sec. 3 table

z = np.linspace(-40, 40, 300001); dz = z[1]-z[0]


def kernel_family(kern):
    """skew a symmetric kernel: p ~ g(z) G(alpha z); match the three moments"""
    g = kern(z); g /= g.sum()*dz
    Gc = np.cumsum(g)*dz; Gc /= Gc[-1]
    def mo(a):
        f = g*np.interp(a*z, z, Gc); n = f.sum()*dz
        m1 = (z*f).sum()*dz/n
        m2 = ((z-m1)**2*f).sum()*dz/n
        m3 = ((z-m1)**3*f).sum()*dz/n
        return m1, sqrt(m2), m3/m2**1.5, np.cumsum(f)*dz/n
    a = fsolve(lambda x: [mo(x[0])[2]-SK], [1.9])[0]
    m1, s1, _, cdf_ = mo(a)
    return SD/s1, Ne0-m1*(SD/s1), cdf_


P  = lambda f: K*np.interp((_R.NE_CEIL-f[1])/f[0], z, f[2])
qn = lambda f, p: f[1]+f[0]*np.interp(p, f[2], z)
qd = lambda f: max(abs(qn(f, p)-q)/q for p, q in QREF.items())

gn, tt, su = [], [], []
for b in KAPPAS:
    f = kernel_family(lambda u, b=b: np.exp(-np.abs(u)**b/b))
    gn.append((100*qd(f), P(f)/PREF, b))
for nu in NUS:
    f = kernel_family(lambda u, nu=nu: tdist.pdf(u, nu))
    tt.append((100*qd(f), P(f)/PREF, nu))
for dl in DELTAS:
    g_, xi_, lam_ = fit(dl)
    p = K*cdf(_R.NE_CEIL, g_, dl, xi_, lam_)
    q = max(abs(ppf(pp, g_, dl, xi_, lam_)-qq)/qq for pp, qq in QREF.items())
    su.append((100*q, p/PREF, dl))
gn.sort(); tt.sort(); su.sort()

# the heavier branch of the generalised normal is the one the section uses;
# kappa > 2 is the lighter side and is plotted dashed
gh = [p for p in gn if p[2] < 2]
gl = [p for p in gn if p[2] > 2]

at = lambda s, t: float(np.exp(np.interp(t, [p[0] for p in s],
                                         [np.log(p[1]) for p in s])))


def table():
    """the matched-deviation table of Sec. 3: one row per deviation"""
    out = []
    for d in DEVS:
        v = [at(s, d) for s in (gh, tt, su)]
        out.append((d, v[0], v[1], v[2], max(v)/min(v)))
    return out


if __name__ == "__main__":
    print('reference: P(N_e <= %d) = %.4g at the ceiling %.1f e-'
          % (_R.NE_CUT, PREF, _R.NE_CEIL))
    print('\n  dev     gen-normal   skew-t   Johnson S_U   spread')
    for d, a, b, c, sp in table():
        print('  %.1f%%   %9.0f %8.0f %13.0f    x%.1f' % (d, a, b, c, sp))
    print('\nLaTeX rows for Sec. 3:\n')
    for d, a, b, c, sp in table():
        print(r'\(%.1f\%%\) & %.0f & %.0f & %.0f & \(\times%.1f\)\\'
              % (d, a, b, c, sp))
    print('\n  The spread column is what Sec. 3 claims: three families spanning')
    print('  two construction mechanisms agree within a factor of about two at')
    print('  matched quantile deviation. It is insensitive to the ceiling; the')
    print('  individual entries are not, at the third significant figure.')
