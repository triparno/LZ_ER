"""Figure: tail probability against the shift in the event-conditioned Ne
quantiles, for three distribution families. Produces lz_ratio_vs_deviation.pdf / .png.

Self-contained: no external data. Needs families.py alongside it, which in
turn needs johnson.py and response.py.

Each family is constrained to reproduce the executable model's mean, standard
deviation and standardised skewness; the remaining shape parameter
(kappa, nu or delta) is scanned. The x axis is the largest fractional shift in the 10/50/90% quantiles of Ne
relative to NEST, all evaluated at the event-conditioned energy. It is NOT a
measure of agreement with the band of Fig. 2 of the Letter, which is built at
fixed S1c from a flat energy spectrum and mixes parent energies; see
notes/band-reconstruction.md. The y axis is P(Ne<=409) relative to NEST.

Knobs worth editing: the annotated members in the ANNOT dicts and XLIM/YLIM.
The scanned grids KAPPAS / NUS / DELTAS and the family construction live in
families.py, which the Sec. 3 table check also uses; NLEAK and the ceiling come
from response.py. Nothing numeric should be defined in this file.
"""
import os
_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper")
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
from math import sqrt, pi
from scipy.stats import t as tdist
from scipy.optimize import fsolve
import sys; sys.path.insert(0, '.')
import response as _R   # also reached via johnson.py below; named here so
                        # the verify.sh preflight can see it
# The families come from families.py, which check_manuscript.py also imports,
# so the figure and the Sec. 3 table cannot disagree. NLEAK and TOTBKG come
# from response.py; they set the "equals the quoted total" line.
import families as _F
from families import gn, gl, gh, tt, su, at, KAPPAS, NUS, DELTAS, SD, Ne0, SK, K, QREF, PREF
NLEAK, TOTBKG = _R.NLEAK, _R.TOTBKG
XLIM, YLIM = (0, 2.0), (2e-3, 3e4)
ANNOT_B, ANNOT_N, ANNOT_D = (1.8, 1.5), (25, 9), (5, 3.5, 2.2)

fig, ax = plt.subplots(figsize=(5,3))
ax.plot([p[0] for p in gh],[p[1] for p in gh],'o-', color='lightcoral', lw=1.9, ms=5,
        label=r'generalised normal ($\kappa$)')
ax.plot([p[0] for p in gl],[p[1] for p in gl],'o--',color='lightcoral', lw=1.4, ms=5, mfc='none')
ax.plot([p[0] for p in tt],[p[1] for p in tt],'^-', color='steelblue',lw=1.9, ms=6,
        label=r'skew-$t$ ($\nu$)')
ax.plot([p[0] for p in su],[p[1] for p in su],'D-',
        color='mediumseagreen', lw=1.9, ms=5,
        label=r'Johnson $S_U$ ($\delta$)')
ax.axhline(1, color='k', lw=1.2); ax.text(0.04, 1.35, 'NEST', fontsize=11)
ax.axhline(TOTBKG/NLEAK, color='#777777', lw=1.1, ls=':')
ax.text(0.8, TOTBKG/NLEAK*0.10, 'ER contribution reaches the\nquoted total background',
        fontsize=11, color='#555555')
ax.set_yscale('log'); ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
ax.set_xlabel('max. shift in the 10, 50, 90% quantiles of $N_e$  [%]',
              fontsize=11)
ax.set_ylabel(r'$P(N_e\leq409)\,/\,P_{\rm NEST}$', fontsize=11)
ax.legend(fontsize=11, loc='lower right', framealpha=.55)
ax.tick_params(axis='both', which='major', direction='in', length=7.5, width=1)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=1)
plt.tight_layout()
plt.savefig(os.path.join(_OUT,'lz_ratio_vs_deviation.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(_OUT,'lz_ratio_vs_deviation.png'), dpi=185, bbox_inches='tight')
print('written. matched-deviation ratios (families.py prints the LaTeX rows):')
for d, a, b, c, sp in _F.table():
    print('  %.1f%%  gen-normal %6.0f   skew-t %6.0f   Johnson %6.0f   spread x%.1f'
          % (d, a, b, c, sp))
