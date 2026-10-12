"""Checks the numbers printed in the manuscript against the code that makes them.

leakage.py verifies the values quoted in the prose. This script does the other
half: it reads ../paper/main.tex, pulls the tables and the headline scalars out
of it, recomputes each, and reports any that disagree. Every table in the paper
is covered, so a row edited by hand and not regenerated is caught here.

Exit status is non-zero if anything disagrees.

    cd analysis && python check_manuscript.py [../paper/main.tex]
"""
import sys, re
import numpy as np
from math import sqrt, pi, log
sys.path.insert(0, '.')
import response as R
from scipy.optimize import fsolve, brentq
from scipy.stats import skewnorm
from scipy.integrate import trapezoid

TEX = sys.argv[1] if len(sys.argv) > 1 else '../paper/main.tex'
try:
    tex = open(TEX).read()
except OSError:
    sys.exit('not found: %s' % TEX)

FAIL, CHECKED = [], 0
def precision(tex_str):
    """The absolute tolerance implied by how many digits the paper prints.
    A value written 0.04 asserts [0.035, 0.045]; one written 522 asserts
    [521.5, 522.5]. Testing against that is the actual criterion -- is the
    printed number the correctly rounded form of the computed one -- and it
    removes the need to hand-tune a relative tolerance per row."""
    t = tex_str.strip().replace('\\%','')
    t = t.replace('\\(','').replace('\\)','').replace('(','').replace(')','').strip()
    m = re.match(r'^([-\d.]+)\\times10\^\{?(-?\d+)\}?$', t)
    mant, exp = (m.group(1), int(m.group(2))) if m else (t, 0)
    dec = len(mant.split('.')[1]) if '.' in mant else 0
    return 0.5*10**(-dec)*10**exp

def chk(label, got, tex_str, tol=0.0):
    """got = computed; tex_str = the field exactly as the manuscript prints it.
    Passes if got rounds to what is printed, or agrees to tol relatively."""
    global CHECKED
    CHECKED += 1
    want = num(tex_str)
    if want is None: return
    atol = precision(tex_str)
    ok = abs(got-want) <= max(atol, tol*abs(want))
    if not ok: FAIL.append('%s (computed %.6g, tex %s)' % (label, got, tex_str.strip()))
    print('  %-50s %-12s tex %-12s %s'
          % (label, '%.5g' % got, tex_str.strip()[:12], 'ok' if ok else '** MISMATCH'))

def num(x):
    """'\\(4.31\\times10^{-10}\\)' or '0.22\\%' or '1.4' -> float.
    The escaped delimiters must go before the bare parentheses, or the
    backslash of \\( survives and nothing parses."""
    x = x.strip().replace('\\%','')
    x = x.replace('\\(','').replace('\\)','')
    x = x.replace('(','').replace(')','').strip()
    m = re.match(r'^([-\d.]+)\\times10\^\{?(-?\d+)\}?$', x)
    if m: return float(m.group(1))*10**int(m.group(2))
    if x in ('---','--',''): return None
    return float(x)

# ---------------------------------------------------------------- family table
z = np.linspace(-40, 40, 400001); dz = z[1]-z[0]
_d  = R.GAM/sqrt(1+R.GAM**2)
_sc = R.SD/sqrt(1-2*_d*_d/pi); _lo = R.MU - _sc*_d*sqrt(2/pi)
ref_pdf = lambda x: skewnorm.pdf(x, R.GAM, loc=_lo, scale=_sc)
NE_EV = 9268./R.g2
NLEAK0, TOTBKG = R.NLEAK, R.TOTBKG
_cache = {}
def member(k):
    if k in _cache: return _cache[k]
    g = np.exp(-np.abs(z)**k/k); g /= g.sum()*dz
    Gc = np.cumsum(g)*dz; Gc /= Gc[-1]
    def mo(a):
        f = g*np.interp(a*z, z, Gc); n = f.sum()*dz
        m1 = (z*f).sum()*dz/n
        m2 = ((z-m1)**2*f).sum()*dz/n
        m3 = ((z-m1)**3*f).sum()*dz/n
        return m1, sqrt(m2), m3/m2**1.5, np.cumsum(f)*dz/n, f/n
    a = fsolve(lambda x: [mo(x[0])[2]-R.SK], [1.9])[0]
    m1, s1, _, cdf, pdf = mo(a)
    s = R.SD/s1; xi = R.MU - m1*s
    _cache[k] = (a, s, xi, cdf, pdf); return _cache[k]
def P_of(k):
    a, s, xi, cdf, _ = member(k)
    return R.K*np.interp((R.NE_CEIL-xi)/s, z, cdf)
def qshift(k):
    a, s, xi, cdf, _ = member(k)
    return max(abs((xi+s*np.interp(q, cdf, z))-v)/v for q, v in R.QREF.items())
def densdev(k):
    a, s, xi, _, pdf = member(k)
    xs = np.linspace(R.QREF[0.10], R.QREF[0.90], 601)
    return max(abs(np.interp((x-xi)/s, z, pdf)/s - ref_pdf(x))/ref_pdf(x) for x in xs)
def dens_ev(k):
    a, s, xi, _, pdf = member(k)
    return np.interp((NE_EV-xi)/s, z, pdf)/s

print('\n=== Sec. 3, the family table ===')
# the last column holds \(...\) forms, so the fields must be matched
# non-greedily to end of line rather than by excluding backslashes
rows = re.findall(r'^(2\.3|2\.2|2\.0|1\.8|1\.7|1\.6|1\.4) & (.+?) & (.+?) & '
                  r'(.+?) & (.+?) & (.+?) & (.+?)\\\\\s*$', tex, re.M)
if len(rows) != 7:
    FAIL.append('family table: found %d rows, expected 7' % len(rows))
    print('  ** could not parse the family table (%d rows)' % len(rows))
else:
    P0, d0 = P_of(2.0), dens_ev(2.0)
    for r0 in rows:
        k = float(r0[0])
        chk('kappa %.1f  quantile shift [pct]' % k, 100*qshift(k),  r0[1])
        chk('kappa %.1f  density dev. [pct]'   % k, 100*densdev(k), r0[2])
        chk('kappa %.1f  P(Ne<409)'  % k, P_of(k),            r0[3])
        chk('kappa %.1f  ratio R'    % k, P_of(k)/P0,         r0[4])
        chk('kappa %.1f  N_leak^ER'  % k, NLEAK0*P_of(k)/P0,  r0[5])
        chk('kappa %.1f  R_ev'       % k, dens_ev(k)/d0,      r0[6], 0.01)

# --------------------------------------------------------- the alpha list
print('\n=== Sec. 3, the matched alpha values ===')
m = re.search(r'\\alpha=\\\{([\d.,]+)\\\}', tex)
if m:
    want = [x for x in m.group(1).split(',')]
    for k, w in zip((2.3,2.2,2.0,1.8,1.7,1.6,1.4), want):
        chk('alpha at kappa = %.1f' % k, member(k)[0], w)
else:
    FAIL.append('alpha list not found'); print('  ** alpha list not found')
m = re.search(r'\\gamma_1\(\\alpha\)=([\d.]+)', tex)
if m: chk('target standardised skewness', R.SK, m.group(1))

# --------------------------------------------------- Sec. 2 response scalars
print('\n=== Sec. 2, the response at the event ===')
pairs = [(r'\\sigma_{N_e}=([\d.]+)\\\)~e\\\(\^-\\\), and it is', R.SD, 'sigma_Ne', 0.005),
         (r'N_q\\;=\\;N_\{ph\}\+N_e\\;=\\;(\d+)', R.NQ, 'N_q', 0.001),
         (r'returns\s*\\\(\\omega=([\d.]+)\\\)', R.omega(R.NQ, True), 'omega printed', 0.005),
         (r'is therefore \\\(([\d.]+)\\\), and', R.omega(R.NQ,True)/R.omega(R.NQ), 'omega ratio', 0.002),
         (r'gives \\\(\\delta=([\d.]+)\\\)', R.GAM/sqrt(1+R.GAM**2), 'delta', 0.0),
         (r'\\omega_\\mathrm{sn}=([\d.]+)\\\) and', R.SD/sqrt(1-2*(R.GAM/sqrt(1+R.GAM**2))**2/pi), 'omega_sn', 0.0),
         (r'\\\(\\xi=([\d.]+)\\\)~e', R.MU-(R.SD/sqrt(1-2*(R.GAM/sqrt(1+R.GAM**2))**2/pi))*(R.GAM/sqrt(1+R.GAM**2))*sqrt(2/pi), 'xi', 0.0),
         (r'returns \\\(P=([\d.]+)\\times10\^\{-7\}\\\)', 1e7*R.PREF, 'P at the event [1e-7]', 0.0),
         (r'\\gamma_\{\\rm ER\}\^\{2\}\}\}=([\d.]+),', R.GAM/sqrt(1+R.GAM**2), 'delta (Sec. 3)', 0.0),
         (r'K\\;&=\\;([\d.]+)\\pm', R.K, 'K', 0.0),
         (r'K\\;&=\\;[\d.]+\\pm([\d.]+)', R.DK, 'K uncertainty', 0.0)]
for pat, got, lbl, tol in pairs:
    m = re.search(pat, tex)
    if m: chk(lbl, got, m.group(1))
    else: print('  (%s not located in the tex)' % lbl)

# ------------------------------------------------------- the threshold ladder
print('\n=== Sec. 4, the threshold ladder ===')
def nleak(k): return NLEAK0*P_of(k)/P_of(2.0)
# anchored at line start: "the quoted total" is a substring of the
# "10\% of the quoted total" row, so an unanchored search finds the wrong one
LADDER = ((r'^\\\(1\\%\\\) of the quoted total', 0.01*TOTBKG, '1% of the quoted total'),
          (r'^\\\(10\\%\\\) of the quoted total', 0.10*TOTBKG, '10% of the quoted total'),
          (r'^the quoted total',                     TOTBKG,      'the quoted total'),
          (r'^\\\(0\.1\\\) counts',                0.1,         '0.1 counts'))
for pat, tgt, lbl in LADDER:
    m = re.search(pat + r'\s*&\s*([\d.]+)\s*&\s*\\\(([\d.]+)\\%\\\)', tex, re.M)
    if not m:
        print('  (row "%s" not located)' % lbl); continue
    k = brentq(lambda x: nleak(x)-tgt, 1.2, 2.4)
    chk('kappa at %s' % lbl, k, m.group(1))
    chk('  its quantile shift [pct]', 100*qshift(k), m.group(2))

# ------------------------------------------------------------- the gap table
print('\n=== Sec. 4, the empty-interval table ===')
try:
    import gap as G
    grows = re.findall(r'^(?:2\.0 \(NEST\)|1\.[4-8]) & ([\d.]+) & ([\d.]+)\\\\', tex, re.M)
    ks = (2.0, 1.8, 1.7, 1.6, 1.5, 1.4)
    if len(grows) == len(ks):
        for k, (c, p) in zip(ks, grows):
            n = G.counts(k, G.NE_EV, G.NE_ED)
            chk('kappa %.1f  expected counts in the gap' % k, n, c)
    else:
        print('  (gap table: %d rows parsed, expected %d)' % (len(grows), len(ks)))
except Exception as e:
    print('  (gap.py not evaluated: %s)' % e)

# --------------------------------------------- the measured smearing factor
print('\n=== Sec. 4, the measured smearing factor ===')
m = re.search(r'\\mathrm\{d\}E\}\;&=\;([\d.]+)\\,', tex)
try:
    d = np.genfromtxt('../scan/scan24_smeared.txt', names=True, comments='#')
    f = trapezoid(d['P_win_s'], d['E'])/trapezoid(d['P_win'], d['E'])
    if m: chk('window integral ratio, measured', f, m.group(1))
    arch = trapezoid(d['P_win'], d['E'])
    chk('its quanta-level integral vs the archived value', arch, '1.497\\times10^{-6}')
except OSError:
    print('  (scan24_smeared.txt absent; skipped)')

# ------------------------------------------ Sec. 4, the Table I accounting
print('\n=== Sec. 4, the Table I accounting ===')
for pat, got, lbl in (
    (r'account for \\\(([\d.]+)\\\) of the', 1646.9, 'covered by the six rows'),
    (r'of the \\\(([\d.]+)\\\) expected', 1706.8, 'Table I expected total'),
    (r'residual is therefore \\\(([\d.]+)\\\) counts', 59.9, 'residual'),
    (r'which \\\(([\d.]+)\\\) are the\s*\n?nuclear-recoil', 0.172, 'NR remainder'),
    (r'The remaining \\\(([\d.]+)\\\) are', 59.7, 'ER components enumerated'),
    (r'from \\\(123.2\\\) to \\\(([\d.]+)\\\) per', 121.2, 'all-tritium bound')):
    m = re.search(pat, tex)
    if m: chk(lbl, got, m.group(1))
    else: print('  (%s not located in the tex)' % lbl)

# -------------------------------- Sec. 3, the window-weighted ratio comparison
print('\n=== Sec. 3, the window-weighted ratios ===')
try:
    import window_ratio as WR
    from scipy.integrate import trapezoid as _tz
    m = re.search(r'are \\\((.+?)\\\) against single-energy', tex)
    if m:
        I0 = _tz([WR.P_at(e, 2.0) for e in WR.EGRID], WR.EGRID)
        want = [x.strip() for x in m.group(1).split(',')]
        if len(want) != len(WR.KAPPAS):
            FAIL.append('Sec. 3 window-weighted list: %d values, expected %d'
                        % (len(want), len(WR.KAPPAS)))
            print('  ** window-weighted list has %d values' % len(want))
        else:
            for k, w in zip(WR.KAPPAS, want):
                Ik = _tz([WR.P_at(e, k) for e in WR.EGRID], WR.EGRID)
                chk('window-weighted ratio at kappa = %.1f' % k, Ik/I0, w)
    else:
        print('  (window-weighted list not located in the tex)')
except Exception as e:
    FAIL.append('Sec. 3 window-weighted ratios not evaluated')
    print('  ** window_ratio.py not evaluated: %s' % e)

# ------------------------------------- the appendix's OldW13eV channel figures
print('\n=== Appendix, the OldW13eV channel ===')
try:
    import wvalue as WV
    for pat, got, lbl in (
        (r'\\partial\\ln P/\\partial\\gamma_\\mathrm{ER}=(-[\d.]+)', WV.deriv,
         'd lnP / d gamma_ER'),
        (r'gives a factor of \\\(([\d.]+)\\\), and the', WV.lin,
         'linearised factor'),
        (r'held fixed, gives \\\(([\d.]+)\\\) once', WV.curv,
         'factor with curvature'),
        (r'the leakage falls by a factor of\s*\n?\\\(([\d.]+)\\\)', WV.MEASURED,
         'measured in NEST'),
        (r'leakage spans a factor of (\d+) over', WV.span,
         'span over gamma_ER in [1.55, 2.25]')):
        m = re.search(pat, tex)
        if m: chk(lbl, got, m.group(1))
        else: print('  (%s not located in the tex)' % lbl)
except Exception as e:
    FAIL.append('appendix OldW13eV figures not evaluated')
    print('  ** wvalue.py not evaluated: %s' % e)

# ----------------------------------- Sec. 3, the cross-family comparison table
print('\n=== Sec. 3, the matched-deviation comparison ===')
try:
    import families as FAM
    frows = re.findall(r'^\\\(([\d.]+)\\%\\\) & (\d+) & (\d+) & (\d+) & '
                       r'\\\(\\times([\d.]+)\\\)\\\\\s*$', tex, re.M)
    want = FAM.table()
    if len(frows) != len(want):
        FAIL.append('Sec. 3 comparison table: found %d rows, expected %d'
                    % (len(frows), len(want)))
        print('  ** could not parse the comparison table (%d rows)' % len(frows))
    else:
        for r0, (d, a, b, c, sp) in zip(frows, want):
            chk('dev %.1f%%  generalised normal' % d, a,  r0[1])
            chk('dev %.1f%%  skew-t'             % d, b,  r0[2])
            chk('dev %.1f%%  Johnson S_U'        % d, c,  r0[3])
            chk('dev %.1f%%  spread'             % d, sp, r0[4])
except Exception as e:
    FAIL.append('Sec. 3 comparison table not evaluated')
    print('  ** families.py not evaluated: %s' % e)

# ------------------------------------------- Sec. 5, the discrimination table
print('\n=== Sec. 5, the calibration-size table ===')
try:
    import discriminate as D
    D.setgrid()
    drows = re.findall(r'^reaches (?:a tenth of the quoted total|the quoted total|'
                       r'\\\(0\.1\\\) counts)\s*&\s*([\d.]+)\s*&\s*(.+?)\s*&\s*'
                       r'(.+?)\s*&\s*(.+?)\\\\\s*$', tex, re.M)
    if len(drows) != 3:
        FAIL.append('Sec. 5 table: found %d rows, expected 3' % len(drows))
        print('  ** could not parse the calibration-size table (%d rows)' % len(drows))
    else:
        ks = [float(r0[0]) for r0 in drows]
        got = {r0[0]: r0 for r0 in D.run(True, kappas=tuple(ks), quiet=True)}
        for r0, k in zip(drows, ks):
            _, _, d_fit, _, d_fit_b, _ = got[k]
            chk('kappa %.2f  N for 3 sigma'   % k, D.n_for(3, d_fit),   r0[1])
            chk('kappa %.2f  N for 5 sigma'   % k, D.n_for(5, d_fit),   r0[2])
            chk('kappa %.2f  N, bulk only'    % k, D.n_for(3, d_fit_b), r0[3])
        # the same two numbers repeated in the summary, which is where they drift
        for pat, k, lbl in ((r'standard deviations by \\\((.+?)\\\) events', 1.66,
                             'summary: N at the quoted-total shape'),
                            (r'a tenth of that background by\s*\n?\\\((.+?)\\\);', 1.87,
                             'summary: N at the tenth-of-total shape')):
            m = re.search(pat, tex)
            if m: chk(lbl, D.n_for(3, got[k][2]), m.group(1))
            else: print('  (%s not located)' % lbl)
except Exception as e:
    FAIL.append('Sec. 5 table not evaluated')
    print('  ** discriminate.py not evaluated: %s' % e)

# ------------------------------------- the Supplement's double-vacancy table
print('\n=== Sec. S4, the double-vacancy modes (supplementary.tex) ===')
try:
    import os
    SUP = os.path.join(os.path.dirname(os.path.abspath(TEX)), 'supplementary.tex')
    sup = open(SUP).read()
    import double_vacancy as DV
    srows = re.findall(r'^([KLM]) & ([\d.]+) & (\d+) & +(\d+) & ([\d.]+) & '
                       r'(.+?)\\\\\s*$', sup, re.M)
    if len(srows) != 3:
        FAIL.append('Sec. S4 table: found %d rows, expected 3' % len(srows))
        print('  ** could not parse the double-vacancy table (%d rows)' % len(srows))
    else:
        for r0, row in zip(srows, DV.rows):
            name, E, Nq, sup_, Ne, sd, disp, P, ge, rr, om = row
            assert name == r0[0], (name, r0[0])
            chk('mode %s  <Ne>'              % name, Ne,   r0[2])
            chk('mode %s  sigma_Ne'          % name, sd,   r0[3])
            chk('mode %s  displacement'      % name, disp, r0[4])
            chk('mode %s  P(Ne<=ceiling)'    % name, P,    r0[5])
except Exception as e:
    FAIL.append('Sec. S4 table not evaluated')
    print('  ** double_vacancy.py not evaluated: %s' % e)

print('\n%d values checked.  %s'
      % (CHECKED, 'ALL CONSISTENT' if not FAIL else 'MISMATCHES: ' + ', '.join(FAIL)))
sys.exit(1 if FAIL else 0)
