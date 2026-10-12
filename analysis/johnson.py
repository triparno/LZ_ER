"""Third family for the tail-shape construction: Johnson S_U.

A transformation family rather than a kernel-times-CDF family: if Z is standard
normal then X = xi + lam*sinh((Z-gam)/del), so
    F(x) = Phi( gam + del*asinh((x-xi)/lam) ).
Four parameters; (gam, xi, lam) are fixed by the same three moment conditions
used for the other families, leaving del as the tail parameter. Unlike the
generalised normal and the skew-t it does NOT contain the skew-normal as a
member, so it tests whether the result depends on the construction at all.
"""
import sys; sys.path.insert(0,'.')
import numpy as np
from math import sqrt,pi
from scipy.stats import norm,skewnorm
from scipy.optimize import fsolve,brentq
import response as _R
# the single convention: see response.py
W, RATIO, C, K = _R.W, _R.RATIO, _R.C, _R.K
Qy, E, Nq, Ni = _R.Qy, _R.E, _R.NQ, _R.NI
Ne0, SD, A0, SK = _R.MU, _R.SD, _R.GAM, _R.SK
CEIL, PREF, QREF = _R.NE_CEIL, _R.PREF, _R.QREF

# Johnson S_U moments, in closed form
def su_moments(gam,dl,xi,lam):
    w=np.exp(1.0/dl**2); Om=gam/dl
    mean=xi-lam*sqrt(w)*np.sinh(Om)
    var =0.5*lam**2*(w-1)*(w*np.cosh(2*Om)+1)
    m3  =-0.25*lam**3*sqrt(w)*(w-1)**2*(w*(w+2)*np.sinh(3*Om)+3*np.sinh(Om))
    return mean,sqrt(var),m3/var**1.5
def fit(dl):
    def eqs(p):
        gam,xi,lam=p
        m,s,sk=su_moments(gam,dl,xi,lam)
        return [m-Ne0,s-SD,sk-SK]
    gam,xi,lam=fsolve(eqs,[-0.4,Ne0,SD*0.9],xtol=1e-12)
    return gam,xi,abs(lam)
cdf=lambda x,g,dl,xi,lam: norm.cdf(g+dl*np.arcsinh((x-xi)/lam))
ppf=lambda p,g,dl,xi,lam: xi+lam*np.sinh((norm.ppf(p)-g)/dl)
pdf=lambda x,g,dl,xi,lam:(dl/(lam*sqrt(2*pi)))/np.sqrt(1+((x-xi)/lam)**2)\
                          *np.exp(-0.5*(g+dl*np.arcsinh((x-xi)/lam))**2)

if __name__=="__main__":
    print("targets:  mean %.1f   sd %.1f   skew %.4f"%(Ne0,SD,SK))
    print("NEST reference: P(Ne<409) = %.3e\n"%PREF)
    print("  delta    gam      xi       lam    | moments check (m, s, skew)      | qdev    P/P_NEST")
    for dl in (5.0,4.0,3.0,2.5,2.2,1.8,1.5,1.3):
        g,xi,lam=fit(dl); m,s,sk=su_moments(g,dl,xi,lam)
        P=K*cdf(CEIL,g,dl,xi,lam)
        qd=max(abs(ppf(p,g,dl,xi,lam)-q)/q for p,q in QREF.items())
        ok="ok" if abs(m-Ne0)<1e-6 and abs(s-SD)<1e-6 and abs(sk-SK)<1e-8 else "**"
        print("  %4.1f  %7.4f %8.1f %7.1f | %8.1f %7.1f %8.4f  %s | %5.2f%%  %9.1f"
              %(dl,g,xi,lam,m,s,sk,ok,100*qd,P/PREF))
