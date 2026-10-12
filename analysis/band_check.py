"""External check of the band construction in band.py.

The construction is only trustworthy if, at kappa = 2 (the executable model
itself), it reproduces the electron-recoil band the Letter publishes. This
script prints the reconstructed band in the Letter's own coordinates so the
three percentiles can be read directly off Fig. 2 and compared by hand.

It cannot do that comparison for you: the band is published as a figure, not
as a table, and no HEPData entry exists for it. What it gives you is the
reconstruction in a form you can lay against the figure.

    python band_check.py            # table across the WS RoI
    python band_check.py --digitise # reads band_published.csv if you make one
"""
import sys, numpy as np
from band import band_quantiles, g2, W, g1, Qy, RATIO

PCT = (0.10, 0.50, 0.90)

def report(S1c_values):
    print('  reconstructed ER band at kappa = 2 (the executable model)\n')
    print('   S1c     ------------- N_e -------------   ---------- log10 S2c ----------   parent E')
    print('   [phd]      10%      50%      90%            10%      50%      90%      [keVee] at median')
    for S1c in S1c_values:
        q = band_quantiles(S1c, 2.0, ps=PCT)
        l = np.log10(q*g2)
        Emed = (S1c/g1 + q[1])*W
        print('   %5.0f   %7.1f  %7.1f  %7.1f        %6.3f   %6.3f   %6.3f     %6.1f'
              % (S1c, q[0], q[1], q[2], l[0], l[1], l[2], Emed))

def compare(path):
    """band_published.csv: S1c_phd,log10S2c_10,log10S2c_50,log10S2c_90
    digitised from Fig. 2 of the Letter."""
    d = np.genfromtxt(path, delimiter=',', names=True)
    print('\n  against %s\n' % path)
    print('   S1c    percentile   published   reconstructed   difference')
    worst = 0.
    for row in np.atleast_1d(d):
        S1c = float(row['S1c_phd'])
        q = np.log10(band_quantiles(S1c, 2.0, ps=PCT)*g2)
        for p, name, col in zip(q, ('10%','50%','90%'),
                                ('log10S2c_10','log10S2c_50','log10S2c_90')):
            pub = float(row[col]); diff = p - pub
            worst = max(worst, abs(diff))
            print('   %5.0f     %4s       %7.3f      %7.3f       %+6.3f'
                  % (S1c, name, pub, p, diff))
    print('\n   largest discrepancy: %.3f dex (%.1f%% in S2c)'
          % (worst, 100*(10**worst - 1)))
    print('   A faithful reconstruction should sit within the line width of')
    print('   Fig. 2, i.e. a few hundredths of a dex. Larger than ~0.05 dex')
    print('   means the construction or its inputs need revisiting before any')
    print('   band deviation computed from it is used.')

if __name__ == '__main__':
    report([100., 200., 300., 400., 500., 540.1, 600.])
    if '--digitise' in sys.argv:
        try:
            compare('band_published.csv')
        except OSError:
            print('\n  band_published.csv not found. To run the comparison,')
            print('  digitise the 10/50/90 percentile curves of Fig. 2 of the')
            print('  Letter into a CSV with columns')
            print('      S1c_phd,log10S2c_10,log10S2c_50,log10S2c_90')
            print('  one row per S1c you can read off cleanly.')
