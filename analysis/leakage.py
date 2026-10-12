"""Recomputes every number quoted in the section on the electron-recoil
contribution to the expected background, and prints it beside the value in the
paper. Exit status is non-zero if any check fails.

Self-contained except for the three digitised spectra, which it reads from
../spectra/ (set SPECTRA below), and for the measured window integral, which
comes from the NEST scan and is taken here as an input constant.

    cd analysis && python leakage.py
"""
import sys, numpy as np
from math import erf, sqrt, exp, log10, pi, log
from scipy.stats import skewnorm
from scipy.integrate import quad, trapezoid
from scipy.optimize import brentq, fsolve

import sys; sys.path.insert(0, '.')
import response as _R

SPECTRA = '../spectra'
WININT  = _R.WININT   # measured, 24-point scan; response.py publishes it
TOTBKG  = _R.TOTBKG     # LZ total for the comparison region

# ----------------------------------------------------------------- response
inv, s2 = 1/sqrt(2*pi), sqrt(2)
g1, g2, W, RATIO, C, K = _R.g1, _R.g2, _R.W, _R.RATIO, _R.C, _R.K
M = dict(m1=12.31, m2=84.91, m3=0.5707, m4=2.804, m5=34.04,
         m7=82.57, m8=4.557, m9=0.2207, m10=0.1272)
Qy = lambda E: C*(M['m1'] + (M['m2']-M['m1'])/(1+(E/M['m3'])**M['m4'])**M['m9']
                + M['m5'] + (0.0-M['m5'])/(1+(E/M['m7'])**M['m8'])**M['m10'])
S4 = [(0.03174, 2.588, 0.3627, -0.3583), (0.06211, 5.200, 1.359, -2.599)]
def omega(Nq, printed=False):
    """Supplement Eq. 4.  printed=False applies the component-wise mode
    normalisation carried in RecombOmegaER; printed=True evaluates the printed
    expression literally.  See printed_omega.py."""
    x, t = log10(Nq), 0.
    for A, c, w, sk in S4:
        n = 1.0
        if not printed:
            mo = c + 2*inv*sk*w/sqrt(1+sk*sk)
            n  = 1/(exp(-.5*(mo-c)**2/w**2)*(1+erf(sk*(mo-c)/(w*s2))))
        t += n*A*exp(-.5*(x-c)**2/w**2)*(1+erf(sk*(x-c)/(w*s2)))
    return t
F = 96.5
a0,c0,c1,E0,E1,E2,E3,F0,F1 = 1.39,4.,22.1,7.7,54.,26.7,6.4,225.,71.
alpha = lambda E: (1/(1+exp((E-E2)/E3))*(a0+c0*exp(-F/F0)*(1-exp(-E/E0)))
                 + 1/(1+exp(-(E-E2)/E3))*c1*exp(-E/E1)*exp(-sqrt(F)/sqrt(F1)))
def state(E, printed=False):
    Nq = E/W; Ni = Nq/(1+RATIO); mu = Qy(E)*E; r = 1-mu/Ni
    return Nq, Ni, mu, sqrt(r*(1-r)*Ni + (omega(Nq, printed)*Ni)**2), alpha(E)

NPH_MIN, NPH_MAX = 500/g1, 600/g1
NE_CEIL = _R.NE_CEIL
EV = (540.1*1.0)/g1 + 9268/g2          # Nq of the event
E_EV = EV*W

def P_window(E, sd_scale=1.0, dalpha=0.0, ceil_only=False, printed=False):
    Nq, Ni, mu, sd, a = state(E, printed)
    sd *= sd_scale; a += dalpha
    d = a/sqrt(1+a*a); w = sd/sqrt(1-(2/pi)*d*d); xi = mu - w*d*sqrt(2/pi)
    hi = NE_CEIL if ceil_only else min(NE_CEIL, Nq-NPH_MIN)
    lo = 0.0     if ceil_only else max(0.0, Nq-NPH_MAX)
    if hi <= lo: return 0.0
    return K*(skewnorm.cdf(hi,a,loc=xi,scale=w) - skewnorm.cdf(lo,a,loc=xi,scale=w))

FAIL = []
def chk(label, got, want, tol=0.03, fmt='%.4g'):
    ok = abs(got-want) <= tol*abs(want) if want else abs(got) <= tol
    if not ok: FAIL.append(label)
    print('  %-46s %-12s paper %-12s %s'
          % (label, fmt % got, fmt % want, 'ok' if ok else '** MISMATCH'))

print('\n--- parent window ---')
lo = brentq(lambda E: E/W - NPH_MIN, 1, 200)
hi = brentq(lambda E: (E/W) - NPH_MAX - NE_CEIL, 1, 200)
chk('support, lower edge [keVee]', lo, 61.3)
chk('support, upper edge [keVee]', hi, 79.1)
chk('S2c ceiling in electrons',    NE_CEIL, _R.NE_CUT+0.5, tol=0.01)
chk('event Nq',                    EV, 5161, tol=0.01)
chk('event E [keVee]',             E_EV, 69.84, tol=0.01)

print('\n--- ROI acceptance ---')
Ecl = brentq(lambda E: Qy(E)*E - NE_CEIL, 5, 60)
chk('ceiling, Qy(E)E = 409 [keVee]', Ecl, 19.99, tol=0.01)
Nq_min = 3.0/g1 + 645./g2
chk('floor, three cuts jointly [keVee]', Nq_min*W, 0.62, tol=0.02)
ACC = Ecl - Nq_min*W
chk('acceptance width [keVee]', ACC, 19.37, tol=0.01)
chk('Qy max = 1/[W(1+Nex/Ni)]', 1/(W*(1+RATIO)), 62.6, tol=0.01)
chk('Qy -> C m2 as E -> 0', C*M['m2'], 85.3, tol=0.01)
chk('crossover Qy = Qy_max [keVee]',
    brentq(lambda E: Qy(E)*E - E/(W*(1+RATIO)), 0.05, 50), 0.97, tol=0.03)

print('\n--- parent spectrum: window/ROI ratios ---')
LO, HI, WLO, WHI = Nq_min*W, Ecl, 61.3, 79.1
def ratio_from_csv(fn, col='intensity'):
    d = np.genfromtxt('%s/%s' % (SPECTRA, fn), delimiter=',', names=True)
    E, y = d['E_keV'], d[col]
    mean = lambda a, b: np.interp(np.linspace(a, b, 4001), E, y,
                                  left=y[0], right=0.0).mean()
    return mean(WLO, WHI)/mean(LO, HI)
specs = [('pb214_gA085_exchange.csv', '214Pb', 1.057),
         ('pb212_gA085_exchange.csv', '212Pb', 0.973),
         ('kr85_exchange_1stforb.csv', '85Kr',  0.949)]
R = {}
for fn, nm, want in specs:
    R[nm] = ratio_from_csv(fn); chk('ratio %s (digitised PRC 102 065501)' % nm, R[nm], want)

mec2 = 510.999
def fermi(Z, T):
    E_ = T+mec2; p = sqrt(max(E_*E_-mec2*mec2, 1e-12)); b = p/E_
    eta = Z/137.036/max(b, 1e-9); return 2*pi*eta/(1-exp(-2*pi*eta))
def fpE(T):
    if T <= 0: return 0.
    E_ = T+mec2; return fermi(56., T)*sqrt(E_*E_-mec2*mec2)*E_
Qbb = 2457.83
dNdT_full = lambda T: 0. if T <= 0 or T >= Qbb else \
    (Qbb-T)**5*quad(lambda u: fpE(u)*fpE(T-u), 0, T, limit=200)[0]
def dNdT_PR(T):
    if T <= 0 or T >= Qbb: return 0.
    t, q = T/mec2, Qbb/mec2
    return t*(t+2)*(q-t)**5*(1+2*t+4*t*t/3+t**3/3+t**4/30)
rat = lambda f: (quad(f, WLO, WHI, limit=250)[0]/(WHI-WLO)) / \
                (quad(f, LO, HI, limit=250)[0]/(HI-LO))
R['136Xe'] = rat(dNdT_full)
chk('ratio 136Xe, Fermi convolution', R['136Xe'], 7.49)
chk('ratio 136Xe, Primakoff-Rosen',   rat(dNdT_PR), 7.88)
R['14C'], R['nuER'] = 0.720, 1.000

print('\n--- dN/dE and the leakage ---')
comp = [('214Pb',1027.7),('212Pb',112.8),('85Kr',200.5),
        ('136Xe',110.0),('14C',55.3),('nuER',140.6)]
tot = 0.
for nm, N in comp:
    w = (N/ACC)*R[nm]; tot += w
    print('    %-7s ROI %6.1f  rate %6.2f  ratio %5.3f  window %6.2f'
          % (nm, N, N/ACC, R[nm], w))
chk('total dN/dE [per keVee]', tot, _R.DNDE, tol=0.003)

# Table I coverage. The EXPECTED column, which is the one to compare a
# reconstruction against: Table I's printed total of 1713 is the FIT column
# and carries the fitted 1.0-count WIMP besides, so comparing an expected
# subtotal against it overstates the residual by six counts. Holding the whole
# table here forces the enumeration in Sec. 4 to balance against it.
TABLE_I = {                      # expected events, Table I of the Letter
    'Internal beta': 1341.0, 'nu ER': 140.6, '136Xe': 110.0, 'CH3T & 14C': 55.3,
    '124Xe': 21.0, '83mKr': 17.1, '125I': 8.9, 'Detector ERs': 8.5,
    'Accidentals': 2.7, '127Xe+125Xe': 1.5,
    'Atmospheric nu': 0.11, '8B+hep nu': 0.057, 'MSSI': 0.0049}
ENUMERATED = ('124Xe', '83mKr', '125I', 'Detector ERs', 'Accidentals',
              '127Xe+125Xe')                      # named in Sec. 4 as excluded
NRLIKE     = ('Atmospheric nu', '8B+hep nu', 'MSSI')
COVERED    = sum(N for _, N in comp)
EXPECTED   = sum(TABLE_I.values())
chk('Table I expected total', EXPECTED, 1706.8, tol=0.001)
chk('  covered by the six rows', COVERED, 1646.9, tol=0.001)
chk('  residual against it', EXPECTED-COVERED, 59.9, tol=0.01)
chk('  the ER components Sec. 4 enumerates',
    sum(TABLE_I[k] for k in ENUMERATED), 59.7, tol=0.01)
chk('  the NR remainder', sum(TABLE_I[k] for k in NRLIKE), 0.17, tol=0.05)
# the enumeration must close: covered + enumerated + NR = expected total
_gap = EXPECTED - COVERED - sum(TABLE_I[k] for k in ENUMERATED) \
                          - sum(TABLE_I[k] for k in NRLIKE)
chk('  enumeration closes (residual gap)', _gap, 0.0, tol=1e-9)

# The CH3T/14C bound. Table I combines CH3T and 14C in one 55.3-count entry
# and gives no split; tritium's endpoint is 18.6 keV, so it contributes nothing
# in a 61-79 keVee window and treating the entry as 14C is the conservative
# choice. The bound is the window rate if the whole entry were tritium.
tot_no14C = tot - (55.3/ACC)*R['14C']
chk('dN/dE if the CH3T+14C entry were all tritium', tot_no14C, 121.2, tol=0.003)
chk('  i.e. a change of [pct]', 100*(1-tot_no14C/tot), 1.7, tol=0.03)

NLEAK = WININT*tot
# Compared against response.py, not against a literal copied from the paper.
# The tolerance is 0.3%, not the 3% that let a 1.75% normalisation error pass
# as 'ok' here while check_manuscript.py caught it.
chk('N_leak^ER [counts]', NLEAK, _R.NLEAK, tol=0.003)
chk('N_leak as fraction of the quoted total', NLEAK/TOTBKG,
    _R.NLEAK/_R.TOTBKG, tol=0.003)

print('\n--- printed vs implemented Eq. 4 (see printed_omega.py) ---')
chk('omega, implementation, at Nq = 5161', omega(5161.), 0.0462, tol=0.01)
chk('omega, printed form,   at Nq = 5161', omega(5161., True), 0.0682, tol=0.01)
chk('omega, implementation, at Nq = 5243', omega(5243.), 0.0465, tol=0.01)
chk('omega, printed form,   at Nq = 5243', omega(5243., True), 0.0685, tol=0.01)
chk('omega ratio (Nq-independent)', omega(5161.,True)/omega(5161.), 1.475, tol=0.01)
chk('  same ratio at Nq = 5243',   omega(5243.,True)/omega(5243.), 1.475, tol=0.01)
chk('sigma_Ne, implementation [e-]', state(E_EV)[3], 204.1, tol=0.01)
chk('sigma_Ne, printed form [e-]',   state(E_EV, True)[3], 299.3, tol=0.01)
P_IMP, P_PRT = P_window(E_EV), P_window(E_EV, printed=True)
chk('P, implementation', P_IMP, 1.32e-7, tol=0.02)
chk('P, printed form',   P_PRT, 4.27e-4, tol=0.02)
chk('leakage-probability ratio', P_PRT/P_IMP, 3.24e3, tol=0.02, fmt='%.3g')
I_IMP = quad(lambda E: P_window(E), WLO, WHI, limit=400)[0]
I_PRT = quad(lambda E: P_window(E, printed=True), WLO, WHI, limit=400)[0]
chk('analytic window integral, impl. [keV]', I_IMP, 1.500e-6, tol=0.005)
chk('  against the measured 24-point value', I_IMP/WININT, 1.0, tol=0.03)
chk('analytic window integral, printed [keV]', I_PRT, 4.836e-3, tol=0.02)
chk('N_leak from the printed form [counts]', I_PRT*tot, 0.60, tol=0.03)
chk('  as a multiple of the quoted total', I_PRT*tot/TOTBKG, 56.2, tol=0.01)

print('\n--- sensitivities (shape held at NEST, sigma varied) ---')
Nqe, Nie, mue, sde, ae = state(E_EV)
P0 = P_window(E_EV)
h = 1e-5
inst = lambda e: (log(P_window(E_EV, 1+e+h)) - log(P_window(E_EV, 1+e)))/h*(1+e)/(1+e)
d0 = (log(P_window(E_EV, 1+h)) - log(P0))/h
chk('dlnP/dlnsigma at zero variation', d0, 31.0, tol=0.02)
for e, want in ((0.05,4.2),(0.10,14.4),(0.15,42.0),(0.20,107.)):
    chk('  sigma +%2d%% -> factor' % (100*e), P_window(E_EV,1+e)/P0, want, tol=0.05)
dins = (log(P_window(E_EV,1.2*(1+h))) - log(P_window(E_EV,1.2)))/h
chk('instantaneous dlnP/dlnsigma at +20%', dins, 20.9, tol=0.03)
chord = log(P_window(E_EV,1.2)/P0)/log(1.2)
chk('chord in ln sigma over [0,+20%]', chord, 25.6, tol=0.02)
for frac, want in ((0.10, 6.2), (1.0, 16.6)):
    f = brentq(lambda e: NLEAK*P_window(E_EV,1+e)/P0 - frac*TOTBKG, 0., 0.5)
    chk('sigma excess reaching %4.2f x quoted bkg [%%]' % frac, 100*f, want, tol=0.02)
for da, want in ((-0.35,6.3),(0.35,0.11)):
    chk('  gamma_ER %+0.2f -> factor' % da, P_window(E_EV,1.0,da)/P0, want, tol=0.08)

print('\n--- density at the event, relative to kappa = 2 ---')
z = np.linspace(-40,40,200001); dz = z[1]-z[0]
d_ = ae/sqrt(1+ae*ae)
SK = (4-pi)/2*(d_*sqrt(2/pi))**3/(1-2*d_*d_/pi)**1.5
NE_EV = 9268./g2
_z, _dz = z, dz
def fam_dens(kap):
    g = np.exp(-np.abs(_z)**kap/kap); g /= g.sum()*_dz
    Gc = np.cumsum(g)*_dz; Gc /= Gc[-1]
    def mo(a):
        f = g*np.interp(a*_z,_z,Gc); n = f.sum()*_dz
        m1 = (_z*f).sum()*_dz/n; m2 = ((_z-m1)**2*f).sum()*_dz/n
        m3 = ((_z-m1)**3*f).sum()*_dz/n
        return m1, sqrt(m2), m3/m2**1.5, f/n
    a = fsolve(lambda x: [mo(x[0])[2]-SK], [1.9])[0]
    m1, s1, _, pdf = mo(a); s = sde/s1; xi = mue - m1*s
    return np.interp((NE_EV-xi)/s, _z, pdf)/s
_d0 = fam_dens(2.0)
for kap, want in ((1.8,49.),(1.7,220.),(1.6,814.),(1.4,7.0e3)):
    chk('density ratio R_ev at kappa = %.1f' % kap, fam_dens(kap)/_d0, want, tol=0.03)

print('\n--- thresholds (via the generalised-normal family) ---')
def fam_ratio(kap):
    g = np.exp(-np.abs(z)**kap/kap); g /= g.sum()*dz
    Gc = np.cumsum(g)*dz; Gc /= Gc[-1]
    def mo(a):
        f = g*np.interp(a*z,z,Gc); n = f.sum()*dz
        m1 = (z*f).sum()*dz/n; m2 = ((z-m1)**2*f).sum()*dz/n
        m3 = ((z-m1)**3*f).sum()*dz/n
        return m1, sqrt(m2), m3/m2**1.5, np.cumsum(f)*dz/n
    a = fsolve(lambda x: [mo(x[0])[2]-SK], [1.9])[0]
    m1, s1, _, cdf = mo(a); s = sde/s1; xi = mue - m1*s
    return K*np.interp((NE_CEIL-xi)/s, z, cdf)/P_window(E_EV, ceil_only=True)
for tgt, want_k, want_dev in ((0.10*TOTBKG,1.87,0.14),(TOTBKG,1.66,0.47),(0.1,1.40,1.10)):
    k = brentq(lambda x: NLEAK*fam_ratio(x) - tgt, 1.2, 2.3)
    chk('kappa at N_leak = %.4g counts' % tgt, k, want_k, tol=0.03)

print('\n' + ('ALL CHECKS PASSED' if not FAIL else 'FAILED: ' + ', '.join(FAIL)))
sys.exit(1 if FAIL else 0)
