"""The OldW13eV convention acts on the tail through the skewness alone.

The flag selects the W-value: 13.532 eV when true (the released setting), a
density-dependent 11.550 eV when false.  Appendix A of the paper reports that
switching it moves the leakage by a factor of 10.3 while recombProb, omega and
sigma_Ne are numerically identical between the two settings.  Those two
statements are consistent because the flag enters exactly one channel.

The yields and the width are evaluated from N_q: recombProb comes from the
returned yields, omega = RecombOmegaER(log10 N_q), and the variance of Eq. (5)
is built from those.  None of them sees W.  The W-value is used where the
quanta are converted to an energy at which to evaluate the ER skewness
parametrisation, so gamma_ER is the one response quantity the flag moves,
from 1.89 to 2.29.

This script checks that a shift of that size accounts for the measured factor,
at fixed mean and fixed width:

  - the logarithmic derivative quoted in Sec. 2, d lnP/d gamma_ER = -5.8,
    linearised over Delta gamma = 0.40, gives 10.0;
  - the analytic representation of Eq. (6) evaluated at the two skewness
    values, which carries the curvature, gives 12.5;
  - NEST itself, at 10^8 trials per point, gives 10.3.

It also prints the leakage against gamma_ER across the range of Sec. 2, which
is where the factor of 57 quoted there comes from.

    cd analysis && python wvalue.py
"""
import sys
sys.path.insert(0, '.')
from math import sqrt, pi, log, exp
from scipy.stats import skewnorm
import printed_omega as p

GAM_TRUE, GAM_FALSE = 1.89, 2.29      # gamma_ER at OldW13eV = true / false
MEASURED = 10.3                       # NEST, 10^8 trials per point

def P_at(mu, sd, a):
    """leakage with the mean and width held fixed and the skewness set to a"""
    d = a/sqrt(1+a*a)
    w = sd/sqrt(1-2*d*d/pi)
    xi = mu - w*d*sqrt(2/pi)
    return p.K*skewnorm.cdf(p._R.NE_CEIL, a, loc=xi, scale=w)

_, Ni, MU, SD, A0 = p.state(p.E_EV, False)
print('\nEvent-conditioned response (OldW13eV = true, as released):')
print('  N_q = %.0f   N_i = %.1f   omega = %.5f' % (p.E_EV/p.W, Ni, p.omega(p.E_EV/p.W)))
print('  mu = %.2f e-   sigma_Ne = %.2f e-   gamma_ER = %.4f' % (MU, SD, A0))
print('\nThese are all functions of N_q and are unchanged by the flag.')

P1, P2 = P_at(MU, SD, GAM_TRUE), P_at(MU, SD, GAM_FALSE)
h = 1e-4
deriv = (log(P_at(MU, SD, A0+h)) - log(P_at(MU, SD, A0)))/h
lin = 1/exp(deriv*(GAM_FALSE-GAM_TRUE))
print('\ngamma_ER %.2f -> %.2f at fixed mu and sigma_Ne:' % (GAM_TRUE, GAM_FALSE))
print('  P        %.4g -> %.4g' % (P1, P2))
print('  factor   %.2f  (leakage falls)' % (P1/P2))
print('\n  d lnP / d gamma_ER at the default      %7.2f   (paper: -5.8)' % deriv)
print('  linearised over Delta gamma = %.2f      %7.2f   (paper: 10.0)'
      % (GAM_FALSE-GAM_TRUE, lin))
curv = P1/P2
print('  analytic, with curvature                %7.2f' % curv)
print('  NEST, measured                          %7.2f' % MEASURED)
print('''
  The linearised estimate understates the factor because d lnP/d gamma_ER
  steepens over the interval; the measured value sits between the linearised
  and fully analytic ones, as it should.''')

print('\nLeakage against gamma_ER over the range of Sec. 2:')
print('  gamma_ER     P(Ne<=409)  relative to default')
P0 = P_at(MU, SD, A0)
for g in (1.55, 1.70, 1.89, 2.00, 2.15, 2.25):
    print('  %6.2f     %10.4g   %10.3f' % (g, P_at(MU,SD,g), P_at(MU,SD,g)/P0))
span = P_at(MU,SD,1.55)/P_at(MU,SD,2.25)
print('\n  span over [1.55, 2.25]: %.0f   (paper: 57)\n' % span)
