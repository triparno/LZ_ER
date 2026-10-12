"""Band deviation computed the way the Letter constructs its band:
flat-in-energy spectrum, percentiles of Ne evaluated at fixed S1c.

At fixed S1c the photon number Nph is fixed, and since recombination
conserves quanta, Nq = Nph + Ne. For a flat energy spectrum the delta
function that enforces this collapses the energy integral exactly:

    p(Ne | Nph)  proportional to  f_{E*}(Ne),   E* = (Nph + Ne) W

so the density at each Ne is the response density at the energy that makes
the quanta add up. No binning or quadrature over E is involved.
"""
import numpy as np
from math import sqrt, pi, erf, exp, log10
from scipy.optimize import fsolve
from scipy.special import erf as erfv

W, RATIO, C, g1, g2 = 13.532e-3, 0.1812, 1.0051, 0.1104, 34.5
M = dict(m1=12.31,m2=84.91,m3=0.5707,m4=2.804,m5=34.04,
         m7=82.57,m8=4.557,m9=0.2207,m10=0.1272)
Qy = lambda E: C*(M['m1']+(M['m2']-M['m1'])/(1+(E/M['m3'])**M['m4'])**M['m9']
                + M['m5']+(0.0-M['m5'])/(1+(E/M['m7'])**M['m8'])**M['m10'])
S4 = [(0.03174,2.588,0.3627,-0.3583),(0.06211,5.200,1.359,-2.599)]
inv, s2 = 1/sqrt(2*pi), sqrt(2)
def om(Nq):
    """Supplement Eq. 4 with the component-wise normalisation of
    RecombOmegaER. Vectorised: Nq may be a scalar or an array."""
    x = np.log10(np.asarray(Nq, dtype=float))
    t = np.zeros_like(x)
    for A, c, w, sk in S4:
        mo = c + 2*inv*sk*w/sqrt(1+sk*sk)
        n  = 1/(exp(-.5*(mo-c)**2/w**2)*(1+erf(sk*(mo-c)/(w*s2))))
        t += n*A*np.exp(-.5*(x-c)**2/w**2)*(1+erfv(sk*(x-c)/(w*s2)))
    return t if t.ndim else float(t)
F = 96.5; a0,c0,c1,E0,E1,E2,E3,F0,F1 = 1.39,4.,22.1,7.7,54.,26.7,6.4,225.,71.
alp = lambda E: (1/(1+exp((E-E2)/E3))*(a0+c0*exp(-F/F0)*(1-exp(-E/E0)))
               + 1/(1+exp(-(E-E2)/E3))*c1*exp(-E/E1)*exp(-sqrt(F)/sqrt(F1)))

# ---- the generalised-normal family, matched in the first three moments ----
u = np.linspace(-30, 30, 200001); du = u[1]-u[0]
_d0 = 1.9/sqrt(1+1.9**2)
SK  = (4-pi)/2*(_d0*sqrt(2/pi))**3/(1-2*_d0*_d0/pi)**1.5
_cache = {}
def shape(k):
    """returns (m1, s1, pdf_on_u) of the unit-scale member with matched skewness"""
    if k in _cache: return _cache[k]
    g = np.exp(-np.abs(u)**k/k); g /= g.sum()*du
    G = np.cumsum(g)*du; G /= G[-1]
    def mo(a):
        f = g*np.interp(a*u, u, G); n = f.sum()*du
        m1 = (u*f).sum()*du/n
        m2 = ((u-m1)**2*f).sum()*du/n
        m3 = ((u-m1)**3*f).sum()*du/n
        return m1, sqrt(m2), m3/m2**1.5, f/n
    a = fsolve(lambda x: [mo(x[0])[2]-SK], [1.9])[0]
    m1, s1, _, pdf = mo(a)
    _cache[k] = (m1, s1, pdf)
    return _cache[k]

def band_quantiles_generic(S1c, dens, ps=(0.10,0.50,0.90), NE=(5.,3500.), n=400001):
    """Percentiles of Ne at fixed S1c for an arbitrary response family.

    dens(Ne, mu, sd) must return the density of Ne for a member of the family
    placed at mean mu and standard deviation sd. Quanta conservation fixes the
    parent energy of each Ne, so the flat-in-energy integral collapses and the
    density is evaluated pointwise.
    """
    Nph = S1c/g1
    Ne  = np.linspace(*NE, n)
    Nq  = Nph + Ne
    E   = Nq*W
    Ni  = Nq/(1+RATIO)
    mu  = Qy(E)*E
    r   = 1 - mu/Ni
    sd  = np.sqrt(np.clip(r*(1-r),0,None)*Ni + (om(Nq)*Ni)**2)
    d   = dens(Ne, mu, sd)
    c = np.cumsum(d); c /= c[-1]
    return np.array([np.interp(p, c, Ne) for p in ps])


def kernel_dens(m1, s1, pdf, grid=None):
    """Density callable for a family supplied as a standardised pdf on `grid`
    (defaults to this module's own grid u)."""
    gr = u if grid is None else grid
    def f(Ne, mu, sd):
        s = sd/s1; xi = mu - m1*s
        return np.interp((Ne-xi)/s, gr, pdf)/s
    return f


def band_quantiles(S1c, k, ps=(0.10,0.50,0.90), NE=(5.,3500.), n=400001):
    m1, s1, pdf = shape(k)
    Nph = S1c/g1
    Ne  = np.linspace(*NE, n)
    Nq  = Nph + Ne                       # quanta conservation at fixed S1c
    E   = Nq*W                           # the energy that makes them add up
    Ni  = Nq/(1+RATIO)
    mu  = Qy(E)*E
    r   = 1 - mu/Ni
    sd  = np.sqrt(np.clip(r*(1-r),0,None)*Ni + (om(Nq)*Ni)**2)
    s   = sd/s1
    xi  = mu - m1*s
    dens = np.interp((Ne-xi)/s, u, pdf)/s
    c = np.cumsum(dens); c /= c[-1]
    return np.array([np.interp(p, c, Ne) for p in ps])

if __name__ == '__main__':
    KS = (2.3,2.2,2.1,1.9,1.8,1.7,1.6,1.5,1.4)
    print('Band deviation at fixed S1c, flat-in-energy (the Letter\'s construction)\n')
    for S1c in (520., 540.1, 580.):
        ref = band_quantiles(S1c, 2.0)
        print('  S1c = %.1f phd   reference 10/50/90 percentiles of Ne: %s'
              % (S1c, ', '.join('%.1f'%x for x in ref)))
        print('     kappa   band dev     (paper, fixed E)')
        paper = {2.3:0.23,2.2:0.15,1.8:0.24,1.7:0.40,1.6:0.59,1.4:1.08}
        for k in KS:
            q = band_quantiles(S1c, k)
            dev = max(abs(a-b)/b for a,b in zip(q,ref))
            print('     %.1f     %6.3f%%      %s'
                  % (k, 100*dev, ('%6.2f%%'%paper[k]) if k in paper else '    --'))
        print()
