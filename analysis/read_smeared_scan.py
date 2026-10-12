"""Reads scan/scan24_smeared.txt and reports what the detection fluctuations do.

The file carries, for the same simulated events, the quanta-level counters of
the original scan and the S1c/S2c counters. The comparison is therefore paired:
P_win_s/P_win measures the smearing directly, and its error follows from the
overlap of the two accepted sets rather than from the two rates separately.

Compare the output of this script with smearing.py, which estimates the same
quantities analytically from the quanta-level scan. If they agree, the
analytic treatment in the paper is vindicated; if not, these numbers are the
ones to quote, and smearing.py says which assumption to look at.

    cd analysis && python read_smeared_scan.py [../scan/scan24_smeared.txt]
"""
import sys
sys.path.insert(0, '.')
import numpy as np
from scipy.integrate import trapezoid
import response as R

path = sys.argv[1] if len(sys.argv) > 1 else '../scan/scan24_smeared.txt'
try:
    d = np.genfromtxt(path, names=True, comments='#')
except OSError:
    sys.exit('not found: %s\n'
             'Run nest-drivers/run_scan24_smeared.sh first; see HOWTO.md.' % path)

need = ('E', 'N', 'n_ceil', 'P_ceil', 'n_win', 'P_win',
        'n_ceil_s', 'P_ceil_s', 'n_win_s', 'P_win_s', 'n_slow')
missing = [c for c in need if c not in d.dtype.names]
if missing:
    sys.exit('columns missing from %s: %s' % (path, ', '.join(missing)))

E = d['E']
print('\n%d energy points, %.3g trials, %.3g on the slow S1/S2 path (%.2g%%)'
      % (len(E), d['N'].sum(), d['n_slow'].sum(),
         100*d['n_slow'].sum()/d['N'].sum()))

print('\n  E      P_ceil quanta   P_ceil S1c/S2c   x      P_win quanta    P_win S1c/S2c   x')
for i in range(len(E)):
    rc = d['P_ceil_s'][i]/d['P_ceil'][i] if d['P_ceil'][i] > 0 else float('nan')
    rw = d['P_win_s'][i]/d['P_win'][i] if d['P_win'][i] > 0 else float('nan')
    print('  %5.1f   %13.4e   %14.4e %5.2f   %13.4e   %14.4e %5.2f'
          % (E[i], d['P_ceil'][i], d['P_ceil_s'][i], rc,
             d['P_win'][i], d['P_win_s'][i], rw))

# Integrate over EVERY sampled energy, zeros included. Masking to the points
# where the quanta-level count is non-zero would discard precisely the
# energies the smearing opens up -- the parent support is wider once the cuts
# are soft, and that widening is part of the effect being measured.
Ih = trapezoid(d['P_win'],   E)
Is = trapezoid(d['P_win_s'], E)
sup_h = E[d['n_win']   > 0]
sup_s = E[d['n_win_s'] > 0]
print('\nparent support where P_win is resolved')
if len(sup_h): print('  quanta level   %.0f -- %.0f keVee' % (sup_h[0], sup_h[-1]))
if len(sup_s): print('  S1c / S2c      %.0f -- %.0f keVee' % (sup_s[0], sup_s[-1]))
print('\nwindow integral over the full sampled range')
print('  quanta level   %.4e keV' % Ih)
print('  S1c / S2c      %.4e keV   (x%.3f)' % (Is, Is/Ih))
DNDE = R.DNDE        # leakage.py derives it
print('\n  the archived quanta-level value is 1.497e-06 keV; this column should')
print('  reproduce it, since the same cuts are applied to the same model.')
if abs(Ih/R.WININT - 1) > 0.002:
    print('\nFAILED: response.py publishes WININT = %.4g keV, this scan gives %.4g.'
          % (R.WININT, Ih))
    sys.exit(1)
print('\n  quanta-level window integral %.4g keV matches response.py to %.2f%%.'
      % (Ih, 100*abs(Ih/R.WININT-1)))
print('\nN_leak^ER = dN/dE x INT P_win dE, with dN/dE = %.1f per keVee:' % DNDE)
print('  quanta level   %.4g counts' % (DNDE*Ih))
print('  S1c / S2c      %.4g counts   (x%.3f)' % (DNDE*Is, Is/Ih))
tot = d['n_win'].sum(); tots = d['n_win_s'].sum()
if tot and tots:
    print('\npooled window counts %d -> %d; the paired ratio is %.3f +- %.3f'
          % (tot, tots, tots/tot,
             (tots/tot)*np.sqrt(1./tot + 1./tots)))
print("""
  The paper quotes the quanta-level numbers and states that they understate
  the expectation. If the ratio above is greater than one, that statement
  holds and the factor replaces the analytic estimate of smearing.py in the
  text. The between-family ratios need family_table.py rerun against these
  acceptances to be updated; the thresholds were shown analytically to be
  insensitive to the smearing, which this file can confirm.
""")
