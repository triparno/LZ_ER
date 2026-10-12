"""Figure: ER leakage probability per decay and its logarithmic sensitivities.

Produces lz_leakage_scan.pdf / .png.

DATA DEPENDENCY
---------------
Reads the 24-point NEST scan from SCAN_FILE, resolved relative to this
script rather than the working directory. Whitespace-separated,
one row per energy point, with the columns this script uses being:

    col 0  : parent energy E            [keVee]
    col 8  : P(Ne <= 409)               measured leakage probability
    col 9  : uncertainty on col 8
    col 10 : hits passing the S1c window (used only as a >0 mask)
    col 11 : P(Ne <= 409 AND 500<S1c<600)
    col 12 : uncertainty on col 11

Regenerate it with run_scan24.sh, which drives lz_scan.cpp. If the file is
absent the script draws the analytic model alone and says so.
"""
import os
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
from math import erf,sqrt,exp,log10,pi,log
from scipy.stats import skewnorm
inv = 1/sqrt(2*pi)

# Everything about the response comes from response.py: W, RATIO, g1, g2, C,
# K, omega() and alpha(), and the ceiling. Nothing numeric is defined here, so
# this figure cannot drift from the numbers in the text.
import sys; sys.path.insert(0, '.')
import response as _R
W, RATIO, g1, g2, C, K = _R.W, _R.RATIO, _R.g1, _R.g2, _R.C, _R.K
Ku, Kd = K + _R.DK, K - _R.DK
Qy, om, alp, CEIL = _R.Qy, _R.omega, _R.alpha, _R.NE_CEIL

def P(E, K=K, dsig=1., dal=0.):
    """the analytic leakage at energy E, as the scan measures it"""
    Nq = E/W; Ni = Nq/(1+RATIO); Ne = Qy(E)*E; r = 1 - Ne/Ni
    sd = dsig*sqrt(r*(1-r)*Ni + (om(Nq)*Ni)**2); a = alp(E) + dal
    d = a/sqrt(1+a*a); wc = sqrt(1-(2/pi)*d*d); sc = sd/wc
    return K*skewnorm.cdf(CEIL, a, loc=Ne-sc*d*2*inv, scale=sc)

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT  = os.path.join(_HERE, "..", "paper")
SCAN_FILE = os.path.join(_HERE, '..', 'scan', 'scan24_results.txt')
HAVE = os.path.exists(SCAN_FILE)
if HAVE:
    d=np.genfromtxt(SCAN_FILE)
    E,Pc,dPc,nw,Pw,dPw = d[:,0],d[:,8],d[:,9],d[:,10],d[:,11],d[:,12]
else:
    print('NOTE: %s not found - drawing the analytic model only.'%SCAN_FILE)
    E=Pc=dPc=nw=Pw=dPw=np.array([])
Eg=np.linspace(44,111,900)

fig,(ax,bx)=plt.subplots(1,2,figsize=(8.0,4.0),sharey=False,gridspec_kw={'width_ratios':[1,1]})
ax.axvspan(64,80,color='goldenrod',alpha=.13,zorder=0)
#ax.text(72,2.4e-7,'parents reaching $S1c\\in[500,600]$',fontsize=8,ha='center',color='darkgoldenrod')
ax.fill_between(Eg,[P(e,Ku) for e in Eg],[P(e,Kd) for e in Eg],color='navy',lw=1.9,zorder=4,alpha=0.5,label='model, $K=%.3f\\pm%.3f$' % (K, _R.DK))
if HAVE:
    ax.errorbar(E,Pc,yerr=dPc,fmt='o',ms=4.6,color='navy',mfc='w',mew=1.3,capsize=2.3,lw=1.2,
                zorder=6,label='NEST $P(N_e\\leq409)$,  $1.1\\times10^{10}$ trials')
    m=nw>0
    ax.errorbar(E[m],Pw[m],yerr=dPw[m],fmt='s',ms=4.6,color='firebrick',mfc='w',mew=1.3,
                capsize=2.3,lw=1.2,zorder=6,label='$+\\,S1c$ window')
ax.set_yscale('log'); ax.set_ylim(3e-9,3.2e-7)
ax.set_xlabel(r'parent energy  $E$  [keV$_{\mathrm{ee}}$]',fontsize=10)
ax.set_ylabel(r'leakage probability per decay',fontsize=11)
ax.legend(fontsize=10,loc='lower left',framealpha=.95)
ax.tick_params(axis='both', which='major', direction='in', length=7.5, width=1)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=1)
ax.set_xlim(45,110)

h=1e-3
bx.plot(Eg,[(np.log(P(e,dsig=1+h))-np.log(P(e)))/h for e in Eg],color='seagreen',lw=2.0,
        label=r'$\partial\ln P/\partial\ln\sigma_{N_e}$')
bx.plot(Eg,[-(np.log(P(e,dal=h))-np.log(P(e)))/h for e in Eg],color='darkorange',lw=2.0,ls='--',
        label=r'$-\,\partial\ln P/\partial\alpha$')
bx.axvspan(64,80,color='goldenrod',alpha=.13)
bx.set_xlabel('parent energy  $E$  [keV$_{ee}$]',fontsize=11)
bx.set_ylabel('logarithmic sensitivity',fontsize=11)
bx.legend(fontsize=11,loc='upper left',framealpha=.95)
bx.set_xlim(44,111); bx.set_ylim(0,44)
bx.tick_params(axis='both', which='major', direction='in', length=7.5, width=1)
bx.tick_params(axis='both', which='minor', direction='in', length=5, width=1)
plt.tight_layout()
plt.savefig(os.path.join(_OUT,'lz_leakage_scan.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(_OUT,'lz_leakage_scan.png'),dpi=190,bbox_inches='tight')
print("sensitivity at 71 keV: dlnP/dlnsigma = %.1f ; dlnP/dalpha = %.2f"
      %((np.log(P(71,dsig=1+h))-np.log(P(71)))/h,(np.log(P(71,dal=h))-np.log(P(71)))/h))
