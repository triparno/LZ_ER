"""Figure: tail-shape degeneracy at fixed mean, variance and skewness.
Produces lz_tail_family.pdf / .png. Self-contained, no external data.

Members of the generalised-normal family of Eq. (9) are matched to the executable
model in the first three moments; kappa controls only the tail decay. kappa=2
reproduces NEST exactly.

Knobs: the kappa list in the two plotting loops, cols (colours), the inset limits,
and ax.set_ylim / set_xlim.
"""
import os
_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper")
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import sqrt,pi,gamma,erf,exp,log10
from scipy.integrate import quad
from scipy.optimize import fsolve,brentq

plt.rcParams["font.family"] = "serif"

import sys; sys.path.insert(0,'.')
import response as _R
# the single convention: see response.py
W, RATIO, C, Qy = _R.W, _R.RATIO, _R.C, _R.Qy
E, Nq, Ni, Ne = _R.E, _R.NQ, _R.NI, _R.MU
SD, G, SKEW, CEIL = _R.SD, _R.GAM, _R.SK, _R.NE_CEIL
d = G/sqrt(1+G*G)

def gn(z,b): return b/(2*gamma(1./b)*b**(1./b))*np.exp(-np.abs(z)**b/b)
def GN(z,b): return quad(lambda t: gn(t,b),-40,z,limit=200)[0]
def dens(z,a,b): return gn(z,b)*GN(a*z,b)


def mom(a,b):
    n=quad(lambda z:dens(z,a,b),-30,30,limit=400)[0]
    m1=quad(lambda z:z*dens(z,a,b),-30,30,limit=400)[0]/n
    m2=quad(lambda z:(z-m1)**2*dens(z,a,b),-30,30,limit=400)[0]/n
    m3=quad(lambda z:(z-m1)**3*dens(z,a,b),-30,30,limit=400)[0]/n
    return m1,sqrt(m2),m3/m2**1.5,n


params={}
for b in (2.3,2.0,1.8,1.6,1.4):
    a=fsolve(lambda x:[mom(x[0],b)[2]-SKEW],[G])[0]
    m1,s1,_,n=mom(a,b); om=SD/s1; xi=Ne-m1*om
    params[b]=(a,om,xi,n)
def pdf(x,b):
    a,om,xi,n=params[b]
    return np.array([dens((xx-xi)/om,a,b)/(om*n) for xx in np.atleast_1d(x)])
def q(b,p):
    a,om,xi,n=params[b]
    return xi+om*brentq(lambda zc: quad(lambda z:dens(z,a,b),-30,zc,limit=400)[0]/n-p,-12,12)

fig,ax=plt.subplots(figsize=(5,3))
cols={2.3:'mediumaquamarine',2.0:'k',1.8:'teal',1.6:'dodgerblue',1.4:'orchid'}
lo,hi=q(2.0,0.10),q(2.0,0.90)
ax.axvspan(lo,hi,color='azure',zorder=0)
# ax.text((lo+hi)/2,2.5e-3,r'NEST 10--90\\% quantiles at E = 69.8 keV$_{\mathrm{ee}}$', ha='center',fontsize=11)
xs=np.linspace(150,2100,900)
for b in (2.3,2.0,1.8,1.6,1.4):
    ax.plot(xs,pdf(xs,b),color=cols[b],lw=2.2 if b==2.0 else 1.5,
            ls='-' if b==2.0 else '--',
            label=('$\\kappa=2.0$ (NEST)' if b==2.0 else '$\\kappa=%.1f$'%b))
ax.axvline(CEIL,color='dimgrey',lw=1.4,ls=':')
ax.annotate('$N_e\\leq409$\n(ROI ceiling)',xy=(CEIL,4e-8),xytext=(650,5e-8),fontsize=11,
            color='dimgrey',arrowprops=dict(arrowstyle='->',color='dimgrey',lw=1))
ax.plot(268.6,2.2e-9,'*',ms=15,color='rebeccapurple',mec='none',mew=.5,zorder=9)
ax.text(296,1.7e-9,'event',fontsize=11)
ax.set_yscale('log'); ax.set_ylim(1.2e-9,7e-3); ax.set_xlim(150,2100)
ax.set_xlabel('$N_e$  [electrons]',fontsize=11)
ax.set_ylabel('probability density',fontsize=11)
ax.legend(fontsize=8.6,loc='upper left',framealpha=.75)
ax.tick_params(axis='both', which='major', direction='in', length=7.5, width=1)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=1)

ins=ax.inset_axes([0.5,0.40,0.33,0.26])
xi2=np.linspace(lo*0.97,hi*1.03,300)
for b in (2.3,2.0,1.8,1.6,1.4):
    ins.plot(xi2,pdf(xi2,b),color=cols[b],lw=1.8 if b==2.0 else 1.1,ls='-' if b==2.0 else '--')
ins.set_xlim(xi2[0],xi2[-1]); ins.set_xticks([1000,1250,1500])
ins.tick_params(labelsize=6.5); ins.set_title(r'core, linear scale',fontsize=11)
ins.tick_params(axis='both', which='major', direction='in', length=5.5, width=1)

plt.tight_layout()
plt.savefig(os.path.join(_OUT,'lz_tail_family.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(_OUT,'lz_tail_family.png'),dpi=190,bbox_inches='tight')

print("10-90%% band: %.0f - %.0f e-"%(lo,hi))
for b in (2.3,2.0,1.8,1.6,1.4):
    a,om,xi,n=params[b]; print("  kappa=%.1f  alpha=%6.3f  omega=%6.1f  xi=%7.1f"%(b,a,om,xi))
