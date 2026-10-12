"""Regenerates the family table of Sec. 3 from the single convention in
response.py: moments matched to the executable model at the event-conditioned
point, kappa scanned, everything else derived.

Columns: the largest fractional shift in the 10/50/90% quantiles of Ne; the
largest fractional deviation of the density across the 10-90% interval; the
tail probability and its ratio to the executable model; and the resulting ER
contribution, ratio x N_leak(kappa=2).

    cd analysis && python family_table.py
"""
import sys; sys.path.insert(0, '.')
import numpy as np
from math import sqrt, pi
from scipy.optimize import fsolve
from scipy.stats import skewnorm
import response as R

NLEAK0 = R.NLEAK                      # N_leak at kappa = 2; leakage.py derives it
KAPPAS = (2.3, 2.2, 2.0, 1.8, 1.7, 1.6, 1.4)
z = np.linspace(-40, 40, 400001); dz = z[1]-z[0]
_sc = R.SD/sqrt(1-2*(R.GAM/sqrt(1+R.GAM**2))**2/pi)
_lo = R.MU - _sc*(R.GAM/sqrt(1+R.GAM**2))*sqrt(2/pi)
ref_pdf = lambda x: skewnorm.pdf(x, R.GAM, loc=_lo, scale=_sc)

def member(k):
    """generalised normal of shape k, matched in the first three moments"""
    g = np.exp(-np.abs(z)**k/k); g /= g.sum()*dz
    Gc = np.cumsum(g)*dz; Gc /= Gc[-1]
    def mo(a):
        f = g*np.interp(a*z, z, Gc); n = f.sum()*dz
        m1 = (z*f).sum()*dz/n
        m2 = ((z-m1)**2*f).sum()*dz/n
        m3 = ((z-m1)**3*f).sum()*dz/n
        return m1, sqrt(m2), m3/m2**1.5, np.cumsum(f)*dz/n, f/n
    a = fsolve(lambda x: [mo(x[0])[2]-R.SK], [1.9])[0]
    m1, s1, sk, cdf, pdf = mo(a)
    s = R.SD/s1; xi = R.MU - m1*s
    return a, s, xi, cdf, pdf, sk

print('\nreference: gamma_ER = %.4f, sigma_Ne = %.3f, ceiling = %.2f e-, P = %.4g'
      % (R.GAM, R.SD, R.NE_CEIL, R.PREF))
print('target standardised skewness = %.6f\n' % R.SK)
# Normalise the ratio on the GRID-CONSTRUCTED kappa = 2 member, not on the
# exact skew-normal R.PREF. The construction carries a grid error of about
# 0.1% in the tail, and dividing by the same construction cancels it from
# every row instead of biasing all of them by it. The kappa = 2 row printing
# exactly 1 is the sign that the cancellation is in place.
assert 2.0 in KAPPAS, 'kappa = 2 must be in KAPPAS to normalise on'


def quantities(k):
    a, s, xi, cdf, pdf, sk = member(k)
    assert abs(sk-R.SK) < 1e-6, (k, sk)
    P  = R.K*np.interp((R.NE_CEIL-xi)/s, z, cdf)
    qs = max(abs((xi+s*np.interp(q, cdf, z))-v)/v for q, v in R.QREF.items())
    xs = np.linspace(R.QREF[0.10], R.QREF[0.90], 601)
    dd = max(abs(np.interp((x-xi)/s, z, pdf)/s - ref_pdf(x))/ref_pdf(x) for x in xs)
    return a, P, qs, dd


P_REF_GRID = quantities(2.0)[1]
print('  ratios normalised on the grid-built kappa = 2 member, which sits')
print('  %+.3f%% from the exact skew-normal; dividing by it cancels that.\n'
      % (100*(P_REF_GRID/R.PREF - 1)))

print('  kappa  alpha   q.shift  dens.dev   P(Ne<409)      ratio    N_leak^ER')
rows = []
for k in KAPPAS:
    a, P, qs, dd = quantities(k)
    rows.append((k, a, 100*qs, 100*dd, P, P/P_REF_GRID, NLEAK0*P/P_REF_GRID))
    print('  %5.1f  %6.3f  %6.2f%%   %6.1f%%  %10.3g  %9.4g  %11.3g'
          % rows[-1])

print('\nLaTeX rows for Sec. 3:\n')
for k, a, qs, dd, P, rat, nl in rows:
    q = '---' if k == 2.0 else '%.2f\\%%' % qs
    d = '---' if k == 2.0 else '%.1f\\%%' % dd
    r = '1'   if k == 2.0 else ('%.4f' % rat if rat < 0.1 else '%.0f' % rat)
    print('%.1f & %s & %s & \\(%.2f\\times10^{%d}\\) & %s & \\(%.1f\\times10^{%d}\\)\\\\'
          % (k, q, d, P/10**int(np.floor(np.log10(P))), int(np.floor(np.log10(P))),
             r, nl/10**int(np.floor(np.log10(nl))), int(np.floor(np.log10(nl)))))

print('\nspread across the tabulated range: %.0f / %.4f = %.3g  (%.2f orders)'
      % (rows[-1][5], rows[0][5], rows[-1][5]/rows[0][5],
         np.log10(rows[-1][5]/rows[0][5])))
alpha_list = ', '.join('%.2f' % r[1] for r in rows)
print('alpha values for the text: {%s}' % alpha_list)
