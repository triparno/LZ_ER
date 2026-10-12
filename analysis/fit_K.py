"""Normalisation of the analytic representation against the NEST scan.

Fits the single factor K in

    P_analytic(E) = K * Phi_sn(NE_CEIL ; mu(E), sigma(E), gamma_1(E))

The ceiling comes from response.py and is NOT 10^4.15/g2. The scan cuts on an
integer electron count, so it measures P(N_e <= 409); a continuous density
reproduces that at the half-integer edge, 409.5. Fitting K at one ceiling and
applying it at another does not cancel and costs 1.75% on every absolute
number -- see the comment in response.py. This script therefore takes the
ceiling from there rather than naming its own, and asserts at the end that the
K it fits is the K response.py publishes, so the two cannot drift apart.

to the measured leakage counts, using the mean, standard deviation and
standardised skewness that NEST itself returned at each scan point (columns
3-5 of scan24_results.txt). Nothing about the yield or width model is
re-derived here: the fit tests only whether the skew-normal form reproduces
what the sampler produced.

Likelihood: the counts are Poisson, n_i ~ Pois(K * p_i * N_i), with p_i the
unit-normalisation analytic probability and N_i the trials. The MLE is then
closed-form,

    K_hat = sum(n_i) / sum(p_i N_i),     var(K_hat) = K_hat / sum(p_i N_i)

Goodness of fit is the Poisson deviance, 2*sum[ n ln(n/mu) - (n-mu) ],
which is the correct statistic for counts this small; the Gaussian chi^2 is
printed alongside for comparison but is not used.

    python fit_K.py [path-to-scan24_results.txt]
"""
import sys, numpy as np
from math import sqrt, pi
from scipy.stats import skewnorm, chi2
from scipy.optimize import brentq

sys.path.insert(0, '.')
import response as R

PATH = sys.argv[1] if len(sys.argv) > 1 else '../scan/scan24_results.txt'
NE_CEIL = R.NE_CEIL

d = np.genfromtxt(PATH)
E, mu, sd, g1s, N, n = d[:,0], d[:,3], d[:,4], d[:,5], d[:,6], d[:,7]

def alpha_from_skew(g):
    """invert gamma_1 = (4-pi)/2 * (delta sqrt(2/pi))^3 / (1-2 delta^2/pi)^1.5"""
    def f(a):
        dl = a/sqrt(1+a*a)
        return (4-pi)/2*(dl*sqrt(2/pi))**3/(1-2*dl*dl/pi)**1.5 - g
    return brentq(f, 1e-6, 30.0)

a = np.array([alpha_from_skew(g) for g in g1s])
dl = a/np.sqrt(1+a*a)
w  = sd/np.sqrt(1-(2/pi)*dl**2)
xi = mu - w*dl*sqrt(2/pi)
p  = skewnorm.cdf(NE_CEIL, a, loc=xi, scale=w)

den = (p*N).sum()
K   = n.sum()/den
dK  = sqrt(K/den)
exp_ = K*p*N

with np.errstate(divide='ignore', invalid='ignore'):
    term = np.where(n > 0, n*np.log(n/exp_), 0.0)
dev = 2*(term - (n - exp_)).sum()
chi2g = (((n - exp_)**2)/exp_).sum()
dof = len(n) - 1

print('Normalisation of the analytic representation against the scan\n')
print('  points %d, trials %.3e, leakage counts %d' % (len(n), N.sum(), n.sum()))
print('  K      = %.4f +/- %.4f   (Poisson MLE, closed form)' % (K, dK))
print('  deviance/dof = %.3f  (%.2f / %d),  p = %.3f'
      % (dev/dof, dev, dof, chi2.sf(dev, dof)))
print('  Gaussian chi2/dof = %.3f  (for comparison only)\n' % (chi2g/dof))
print('   E[keV]   trials      n   expected   (n-exp)/sqrt(exp)')
for i in range(len(n)):
    print('   %6.1f  %.2e  %5d   %8.2f        %+6.2f'
          % (E[i], N[i], n[i], exp_[i], (n[i]-exp_[i])/sqrt(exp_[i])))
print('\n  All 24 points enter; none has zero counts in the ceiling column.')
print('  Ceiling %.1f e- (P(N_e <= %d), continuity corrected), from response.py.'
      % (NE_CEIL, R.NE_CUT))

# The constant response.py publishes must be the one fitted here. Rounding to
# four places is what response.py stores, so compare at that precision.
if abs(round(K, 4) - R.K) > 1e-9 or abs(round(dK, 4) - R.DK) > 1e-9:
    print('\nFAILED: response.py publishes K = %.4f +/- %.4f, this fit gives '
          '%.4f +/- %.4f.\n        Update response.py, then re-run every check: '
          'K scales every absolute number.' % (R.K, R.DK, K, dK))
    sys.exit(1)
print('  Matches response.py (K = %.4f +/- %.4f).' % (R.K, R.DK))
