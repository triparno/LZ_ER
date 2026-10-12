"""How many calibration events would separate the family members of Sec. 3.

Sec. 5 asks for the low-side tail measured to the 10^-5 level. That is not the
only route: the members differ in the bulk as well as the tail -- the density
deviation column of the family table runs from 3% to 18% -- so a shape fit to
ordinary calibration statistics might separate them without reaching the tail
at all. This script asks how many events that would take.

The question has an optimistic answer and a defensible one, and the difference
is the whole point.

  Optimistic. Treat the mean, width and skewness as known and ask how well the
  remaining shape difference can be told apart.

  Defensible. A real calibration FITS the mean, width and skewness from the
  same events. So the test is: take the heavier member as truth, fit the
  implemented response to it with its three moments free, and ask what
  divergence survives. By construction the members share those moments
  exactly, so the fit starts at the truth -- but it does not have to stay
  there, and it will move to wherever it best absorbs the shape difference.

Both are computed by the Asimov form

    Z^2 = 2 N D( q || p_hat ),   D the Kullback-Leibler divergence,

p_hat being the implemented response, at the reference moments in the first
case and at its best-fitting moments in the second. Each is evaluated at fixed
N_q, as a calibration slice would be, and with and without the measured charge
resolution.

The third number is the one that decides whether the answer means anything.
The divergence is reported twice: over the whole response, and restricted to
the central 99.98% of the truth. If almost all of it comes from outside that
range, then the "bulk" test is the tail test wearing a different hat, and the
event count is not an alternative to Sec. 5 but a restatement of it. Printing
the split is the point; the number alone would be misleading.

Two things this cannot settle, and the paper must say so: a real beta
calibration is a broad spectrum, so the number of events in a slice at
N_q ~ 5e3 is a fraction of the total that is not published in usable form; and
whether such data exist at that N_q at all is a question for the collaboration.

    cd analysis && python discriminate.py
    cd analysis && python discriminate.py --grid   # grid convergence only
"""
import sys
sys.path.insert(0, '.')
import numpy as np
from math import sqrt, pi
from scipy.optimize import fsolve, minimize
from scipy.stats import skewnorm, norm
import response as R

KAPPAS = (1.90, 1.87, 1.80, 1.70, 1.66, 1.60, 1.40, 1.39)  # 1.87, 1.66 and 1.39
                                                           # are the Sec. 4 thresholds
EPS, SE_RES = 0.700, 0.20          # from the probe; see smearing.py
SE = R.g2/EPS
TRUNC = 1e-4                       # the bulk is the central 1-2*TRUNC of q
REPORTED_MAX = 1e6                 # nothing this script reports comes near this;
                                   # the numerical floor measured by --selftest
                                   # must sit well above it, and the main run
                                   # checks that no reported value crosses it

# The charge grid. Wide enough that both tails of every family member fall
# below double precision inside it, and fine enough that the 204 e- width is
# resolved several hundred times over; --grid checks that it is.
LO, HI, STEP = 0., 8000., 2.0
NE = None

# The standardised grid the generalised-normal family is built on. +/-150 is
# past the point where exp(-|z|^k/k) underflows for every k used here, so
# interpolation never has to extrapolate and the clamping np.interp would
# otherwise do cannot invent a flat tail.
z = np.linspace(-150., 150., 150001); dz = z[1]-z[0]
_mem, _Kmat = {}, None


def setgrid(step=STEP, lo=LO, hi=HI):
    global NE, _Kmat
    NE = np.arange(lo, hi+step/2, step)
    _mem.clear(); _Kmat = None


def gn_pdf(k, mu, sg, gam):
    """generalised normal of shape k, skewed, with the given first three moments

    Azzalini's construction, 2 g(z) G(a z), with a set by the skewness and the
    location and scale by the mean and width.  The shape solve is cached on k.
    """
    if k not in _mem:
        g = np.exp(-np.abs(z)**k/k); g /= g.sum()*dz
        Gc = np.cumsum(g)*dz; Gc /= Gc[-1]
        _mem[k] = (g, Gc, {})
    g, Gc, acache = _mem[k]

    def mo(a):
        f = g*np.interp(a*z, z, Gc); n = f.sum()*dz
        m1 = (z*f).sum()*dz/n
        m2 = ((z-m1)**2*f).sum()*dz/n
        m3 = ((z-m1)**3*f).sum()*dz/n
        return m1, sqrt(m2), m3/m2**1.5, f/n

    if gam not in acache:
        acache[gam] = float(fsolve(lambda x: [mo(float(x[0]))[2]-gam], [1.9])[0])
    m1, s1, _, pdf = mo(acache[gam])
    s = sg/s1; xi = mu - m1*s
    return np.interp((NE-xi)/s, z, pdf)/s


SKMAX = (4-pi)/2*(sqrt(2/pi))**3/(1-2/pi)**1.5   # the skew-normal's limit
_at_bound = [False]


def sn_pdf(mu, sg, gam):
    """the implemented response: the skew-normal, kappa = 2, moments given

    Inverted in closed form rather than numerically: with b = delta sqrt(2/pi),
    gamma = ((4-pi)/2) (b^2/(1-b^2))^{3/2}, so delta^2 = (pi/2) m/(m+c) with
    m = |gamma|^{2/3} and c = ((4-pi)/2)^{2/3}.  The skew-normal cannot reach
    |gamma| > %.4f; a fit that wants more is clipped just inside, and the clip
    is recorded so that a minimiser sitting on the boundary is not read as a
    converged interior fit.
    """ % SKMAX
    s = 1. if gam >= 0 else -1.
    if abs(gam) >= SKMAX:
        _at_bound[0] = True
    m = min(abs(gam), 0.9999*SKMAX)**(2./3.)
    c = ((4-pi)/2)**(2./3.)
    d = s*sqrt(pi/2*m/(m+c))
    a = d/sqrt(1-d*d)
    w = sg/sqrt(1-2*d*d/pi); xi = mu - w*d*sqrt(2/pi)
    return skewnorm.pdf(NE, a, loc=xi, scale=w)


def sd_e(n):
    """the measured charge resolution: extraction binomial + single-electron"""
    return np.sqrt(np.maximum(n, 1.)*(EPS*(SE_RES*SE)**2
                                      + EPS*(1-EPS)*SE**2))/R.g2


def smear(p):
    """fold in that resolution, which varies with Ne, so the kernel is not a
    convolution and has to be carried as a matrix.  Banded would be cheaper;
    dense is 2501^2 at the default step, which is small enough not to bother."""
    global _Kmat
    if _Kmat is None:
        step = NE[1]-NE[0]
        _Kmat = np.empty((len(NE), len(NE)))
        for j, n0 in enumerate(NE):
            _Kmat[:, j] = norm.pdf(NE, n0, sd_e(n0))*step
    out = _Kmat @ p
    return out/np.trapezoid(out, NE)


def norml(p):
    return p/np.trapezoid(p, NE)


_floored = [False]


def kl(q, p, mask=None):
    """D(q||p) on the support of q, optionally restricted"""
    m = q > 0
    if mask is not None:
        m = m & mask
    pp = np.maximum(p, 1e-300)
    if (p[m] < 1e-300).any():
        _floored[0] = True
    return float(np.trapezoid(q[m]*np.log(q[m]/pp[m]), NE[m]))


def bulk(q):
    """the central 1-2*TRUNC of q, which is what a calibration of ordinary
    size actually populates"""
    c = np.concatenate(([0.], np.cumsum((q[1:]+q[:-1])/2*np.diff(NE))))
    c /= c[-1]
    return (c > TRUNC) & (c < 1-TRUNC)


def n_for(Z, D):
    return Z*Z/(2*D) if D > 0 else float('inf')


def run(sm, kappas=KAPPAS, quiet=False):
    rows = []
    p0_cache = {}
    for k in kappas:
        q = norml(gn_pdf(k, R.MU, R.SD, R.SK))
        if sm:
            q = smear(q)
        if 'p0' not in p0_cache:
            p = norml(sn_pdf(R.MU, R.SD, R.SK))
            p0_cache['p0'] = smear(p) if sm else p
        p0 = p0_cache['p0']
        b = bulk(q)

        d_fix, d_fix_b = kl(q, p0), kl(q, p0, b)

        def nll(t):
            p = norml(sn_pdf(R.MU*t[0], R.SD*t[1], R.SK*t[2]))
            if sm:
                p = smear(p)
            return kl(q, p)

        res = minimize(nll, [1., 1., 1.], method='Nelder-Mead',
                       options=dict(xatol=1e-6, fatol=1e-14, maxiter=6000))
        phat = norml(sn_pdf(R.MU*res.x[0], R.SD*res.x[1], R.SK*res.x[2]))
        if sm:
            phat = smear(phat)
        d_fit, d_fit_b = res.fun, kl(q, phat, b)
        rows.append((k, d_fix, d_fit, d_fix_b, d_fit_b, res.x))

    if quiet:
        return rows
    tag = 'with the measured charge resolution' if sm else 'at the quanta level'
    print('\n--- %s ---' % tag)
    print('            N for 3 sigma           N for 5 sigma       bulk only   cost of   fitted moments        skew shift')
    print('  kappa     fixed      fitted      fixed      fitted    fitted 3s  floating  mu/mu0 sd/sd0  g/g0   in sigma_g')
    for k, dx, df, dxb, dfb, x in rows:
        n3 = n_for(3, df)
        # the skewness a calibration of that size can itself determine
        sg = sqrt(6.0/n3) if n3 < float('inf') else float('inf')
        print('  %5.2f  %9.3g  %10.3g  %9.3g  %10.3g  %9.3g  %8.1f  %.4f %.4f %.4f  %9.2f'
              % (k, n_for(3, dx), n3, n_for(5, dx), n_for(5, df),
                 n_for(3, dfb), n3/n_for(3, dx), x[0], x[1], x[2],
                 abs(R.SK*(x[2]-1))/sg))
    print('  bulk fraction of the fitted divergence: ' +
          ', '.join('%.2f at %.2f' % (r[4]/r[2] if r[2] > 0 else 0., r[0])
                    for r in rows))
    return rows


def selftest(strict=True):
    """kappa = 2 IS the implemented response, so the divergence must vanish; what
    it comes to instead is the numerical floor of the whole calculation.

    With strict=True this exits non-zero if the floor is not far above every
    event count the script reports, or if the constructed truths do not carry
    the moments they were built to carry, so that it works as a check and not
    merely as a printout.
    """
    setgrid()
    q = norml(gn_pdf(2.0, R.MU, R.SD, R.SK))
    p = norml(sn_pdf(R.MU, R.SD, R.SK))
    b = bulk(q)
    d = kl(q, p)
    print('kappa = 2 against the skew-normal it should equal:')
    print('  max |q/p - 1| over the bulk   %.3g' % np.abs(q[b]/p[b]-1).max())
    print('  D(q||p)                       %.3g' % d)
    print('  -> N for 3 sigma              %.4g   (the noise floor; every number'
          % n_for(3, d))
    print('                                        below must sit far under it)')
    bad = []
    floor = n_for(3, d)
    if floor < 100*REPORTED_MAX:
        bad.append('the noise floor is N = %.3g, within a hundredfold of the '
                   'ceiling %.3g on what the script reports' % (floor, REPORTED_MAX))
    print('\nmoments of the constructed truths against their targets:')
    for k in KAPPAS:
        qq = norml(gn_pdf(k, R.MU, R.SD, R.SK))
        m1 = np.trapezoid(NE*qq, NE)
        m2 = np.trapezoid((NE-m1)**2*qq, NE)
        m3 = np.trapezoid((NE-m1)**3*qq, NE)
        e = (m1/R.MU-1, sqrt(m2)/R.SD-1, m3/m2**1.5/R.SK-1)
        print('  kappa %4.2f   mu %+.1e   sigma %+.1e   gamma_1 %+.1e  (relative)'
              % (k, *e))
        if max(abs(x) for x in e) > 1e-4:
            bad.append('kappa %.2f does not carry its target moments' % k)
    if bad:
        print('\nFAILED:\n  ' + '\n  '.join(bad))
        if strict:
            sys.exit(1)
    elif strict:
        print('\nok: the floor is %.3g, a factor %.0f above the %.3g ceiling on '
              'anything' % (floor, floor/REPORTED_MAX, REPORTED_MAX))
        print('    reported, and every constructed truth carries its three '
              'moments to 1e-4.')
    return floor


if __name__ == "__main__":
    setgrid()
    print('reference: mu = %.2f e-, sigma = %.2f e-, gamma_1 = %.4f at N_q = %.0f'
          % (R.MU, R.SD, R.SK, R.NQ))
    print('grid: %.0f to %.0f e- in steps of %.1f (%d points); '
          'bulk = central %.4g%% of the truth'
          % (NE[0], NE[-1], NE[1]-NE[0], len(NE), 100*(1-2*TRUNC)))

    if '--selftest' in sys.argv:
        selftest()
        sys.exit(0)

    if '--grid' in sys.argv:
        print('\ngrid convergence, quanta level, N for 3 sigma with moments fitted:')
        base = None
        for st in (4.0, 2.0, 1.0):
            setgrid(step=st)
            rows = run(False, kappas=(1.66, 1.40), quiet=True)
            vals = [n_for(3, r[2]) for r in rows]
            print('  step %4.1f e-   kappa 1.66: %9.4g    kappa 1.40: %9.4g%s'
                  % (st, vals[0], vals[1],
                     '' if base is None else
                     '   (%.2f%%, %.2f%% from the step above)'
                     % (100*abs(vals[0]/base[0]-1), 100*abs(vals[1]/base[1]-1))))
            base = vals
        sys.exit(0)

    rows_q = run(False)
    rows_s = run(True)

    print('\n  ratio of the smeared to the quanta-level event count, moments fitted: ' +
          ', '.join('%.2f' % (n_for(3, s[2])/n_for(3, q[2]))
                    for q, s in zip(rows_q, rows_s)))
    bad = []
    if _at_bound[0]:
        bad.append('the skewness fit reached the skew-normal boundary %.4f in at '
                   'least one case, so that row is a bound and not a fit' % SKMAX)
    if _floored[0]:
        bad.append('the implemented response underflowed inside the support of the '
                   'truth, so the fitted divergence is a lower bound there')
    worst = max(n_for(3, r[2]) for r in rows_q + rows_s)
    if worst > REPORTED_MAX:
        bad.append('a reported count, %.3g, exceeds the ceiling %.3g that '
                   '--selftest validates the numerical floor against'
                   % (worst, REPORTED_MAX))

    print("""
  The fitted column is the number to quote. The floating column is what the
  three nuisances cost: a calibration that must determine the mean, width and
  skewness from the same events pays that factor over one that knows them.

  The last three columns are where the fit goes. A fitted moment differing from
  unity by more than the calibration itself could determine means the shape
  difference is being absorbed into the moments rather than resolved, and the
  separation is correspondingly harder than the fixed column suggests.

  The bulk column is the honest version of the headline. It is the same fitted
  test with the outer %.2g of the truth at each end removed -- the part of the
  response a calibration of this size would not populate. Where the bulk
  fraction printed above is small, the apparent bulk discrimination is in fact
  tail discrimination, and does not substitute for Sec. 5.

  N is per slice at fixed N_q. A real beta calibration is a spectrum, so the
  events available near N_q = 5e3 are a fraction of the total that the
  released data do not give; this is a scaling statement, not a limit.
""" % TRUNC)

    if bad:
        print('FAILED:\n  ' + '\n  '.join(bad))
        sys.exit(1)
