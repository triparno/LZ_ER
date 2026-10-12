"""The double-vacancy modes of Sec. S4: do they reach the WS RoI ceiling?

125-I decays by electron capture, followed by a 35.49 keV transition that is
internally converted in 93% of decays. The K, L and M modes deposit 67.3, 40.4
and 36.5 keV, all well above the 20.0 keVee at which an ordinary electron
recoil of mean yield would reach the ceiling. So the features Table I of the
Letter places above 25 keVee cannot be ordinary recoils of those energies, and
the question is how far into the RoI the charge-suppressed modes reach.

This reproduces the table in the Supplement. Each mode carries a measured
charge suppression Qy/Qy_beta, applied to the beta yield at the mode's own
energy; sigma_Ne is the full width including the binomial term, with omega from
the same Eq. 4 branch the rest of this archive uses, which is not a constant
across these energies. The displacement is quoted two ways: in units of
sigma_Ne, and as the Gaussian-equivalent Phi^-1(1-P), which is larger because
the response is positively skewed and its low-side tail correspondingly thinner
than Gaussian. Reading the first as if it were the second is the error this
table exists to forestall.

The ceiling comes from response.py, so this table moves with the rest of the
archive rather than carrying its own. check_manuscript.py tests the Supplement's
printed table against what this script computes.

    cd analysis && python double_vacancy.py
"""
import sys
sys.path.insert(0, '.')
from math import sqrt, pi
from scipy.stats import skewnorm, norm
import response as R

# mode, deposited energy [keV], measured charge suppression Qy/Qy_beta
MODES = (('K', 67.3, 0.773), ('L', 40.4, 0.783), ('M', 36.5, 0.785))

print(__doc__.split('\n\n')[0])
print('\n  ceiling %.1f e- (P(N_e <= %d)), from response.py\n'
      % (R.NE_CEIL, R.NE_CUT))
print('  mode   Qy/Qyb   <Ne>   sigma_Ne   (<Ne>-ceil)/sigma   P(Ne<=ceil)   '
      'Phi^-1(1-P)   1-<Ne>/Ni')
rows = []
for name, E, sup in MODES:
    Nq = E/R.W
    Ni = Nq/(1+R.RATIO)
    Ne = sup*R.Qy(E)*E
    r = 1 - Ne/Ni
    om = R.omega(Nq)
    sd = sqrt(r*(1-r)*Ni + (om*Ni)**2)
    gam = R.alpha(E)
    d = gam/sqrt(1+gam*gam)
    sc = sd/sqrt(1-2*d*d/pi)
    xi = Ne - sc*d*sqrt(2/pi)
    P = skewnorm.cdf(R.NE_CEIL, gam, loc=xi, scale=sc)
    disp = (Ne-R.NE_CEIL)/sd
    rows.append((name, E, Nq, sup, Ne, sd, disp, P, -norm.ppf(P), r, om))
    print('   %s     %.3f   %5.0f   %6.0f       %6.2f         %9.1e      %6.1f'
          '        %.3f' % (name, sup, Ne, sd, disp, P, -norm.ppf(P), r))

print('\n  omega is not a constant across these modes: %.4f at N_q = %.0f '
      'down to' % (rows[0][10], rows[0][2]))
print('  %.4f at N_q = %.0f.' % (rows[-1][10], rows[-1][2]))
print('\n  LaTeX rows for the Supplement:\n')
for name, E, Nq, sup, Ne, sd, disp, P, ge, r, om in rows:
    e = int(('%e' % P).split('e')[1])
    print('%s & %.3f & %.0f & %3.0f & %.2f & \\(%.1f\\times10^{%d}\\)\\\\'
          % (name, sup, Ne, sd, disp, P/10**e, e))
print('\n  recombination fractions 1-<Ne>/Ni: '
      + ', '.join('%.3f' % r[9] for r in rows))
print('  Gaussian-equivalent displacements:  '
      + ', '.join('%.1f' % r[8] for r in rows))
print("""
  The point of the table is the gap between the last two columns. The L and M
  modes sit 1.5 and 1.2 sigma_Ne above the ceiling but 1.9 and 1.4 Gaussian
  equivalents, so they are not excluded from the RoI by the response; what
  excludes them is the S1c window, which Sec. S4 treats separately.
""")
