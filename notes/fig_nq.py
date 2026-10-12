"""Figure: contours of constant Nq in the observable plane, both panels.
Produces lz_constant_nq.pdf / .png. Self-contained, no external data.

Panel (a) is linear in S2c, where the contours are straight lines of slope
-g2/g1; panel (b) is logarithmic, where the same contours are curved. The two
arrows from the event are the constant-Nq displacement and the fixed-S1c
displacement.

Knobs: the (E,label) list controlling which contours are drawn, the Eg grid for
the ER band, the axis limits in the per-panel loop, and the annotate offsets.
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import erf,sqrt,exp,log10,pi

plt.rcParams["font.family"] = "serif"


inv=1/sqrt(2*pi); s2=sqrt(2)
g1,g2,W,RATIO,C = 0.1104,34.5,13.532e-3,0.1812,1.0051
M=dict(m1=12.31,m2=84.91,m3=0.5707,m4=2.804,m5=34.04,m7=82.57,m8=4.557,m9=0.2207,m10=0.1272)
Qy=lambda E: C*(M['m1']+(M['m2']-M['m1'])/(1+(E/M['m3'])**M['m4'])**M['m9']
              + M['m5']+(0.0-M['m5'])/(1+(E/M['m7'])**M['m8'])**M['m10'])
S4=[(0.03174,2.588,0.3627,-0.3583),(0.06211,5.200,1.359,-2.599)]



def om(Nq):
    x=log10(Nq);t=0
    for A,c,w,sk in S4:
        mo=c+2*inv*sk*w/sqrt(1+sk*sk)
        n=1/(exp(-.5*(mo-c)**2/w**2)*(1+erf(sk*(mo-c)/(w*s2))))
        t+=n*A*exp(-.5*(x-c)**2/w**2)*(1+erf(sk*(x-c)/(w*s2)))
    return t


def st(E):
    Nq=E/W; Ni=Nq/(1+RATIO); Ne=Qy(E)*E; r=1-Ne/Ni
    return Nq,Ni,Ne,sqrt(r*(1-r)*Ni+(om(Nq)*Ni)**2)


EV=(540.1,9268.0); NQ_EV=EV[0]/g1+EV[1]/g2
Nq0,Ni0,Ne0,sd0=st(NQ_EV*W)
NOM=((Nq0-Ne0)*g1, Ne0*g2)
# band median at fixed S1c = 540
from scipy.optimize import brentq
E91=brentq(lambda E:(E/W-Qy(E)*E)*g1-540.1,5.,300.)
BM=(540.1, Qy(E91)*E91*g2)
Eg=np.arange(6,136,0.5)
band=np.array([[ (st(E)[0]-(st(E)[2]+k*st(E)[3]))*g1 for k in (-2,-1,0,1,2)]+
               [ (st(E)[2]+k*st(E)[3])*g2 for k in (-2,-1,0,1,2)] for E in Eg])

fig,axes=plt.subplots(1,2,figsize=(10,4))
for ax,logy in zip(axes,(False,True)):
    ax.tick_params(axis='both', which='both', direction='in')
    # constant-Nq loci
    for E,lab in ((40,'40'),(55,'55'),(NQ_EV*W,None),(85,'85'),(105,'105')):
        Nq=E/W; ne=np.linspace(0.5,Nq*0.85,600); s1=(Nq-ne)*g1; s2v=ne*g2
        ev=abs(E-NQ_EV*W)<1e-6
        ax.plot(s1,s2v,color='darkorange' if ev else '0.72',
                lw=2.4 if ev else 1.0, ls='-' if ev else '--',
                zorder=3 if ev else 1,
                label=r'constant $N_q=5161$ (event)' if ev else None)
        if lab:
            j=np.argmin(abs(s1-((Nq-Qy(E)*E)*g1)))
            if 20<s1[j]<640 and (4e2<s2v[j]<1.0e5 if logy else 2e3<s2v[j]<7.0e4):
                ax.text(s1[j]+8,s2v[j]*(0.80 if logy else 1.0)-(0 if logy else 3200),
                        lab,fontsize=11,color='0.45')
    # ROI box
    ax.add_patch(plt.Rectangle((3,10**2.75),597,10**4.15-10**2.75,fill=False,
                 ec='seagreen',lw=1.4,ls='-',zorder=2))
    ax.text(250,10**4.*(1.0 if logy else 1.0)+(0 if logy else 600),'WS RoI',
            fontsize=11,color='seagreen',ha='center')
    ax.axvspan(500,600,color='goldenrod',alpha=.13,zorder=0)
    # Response at FIXED PARENT ENERGY. At each E, Ne is displaced by +/- sigma_Ne
    # with Nph = Nq - Ne compensating, so every +/-1 sigma point lies on the
    # constant-Nq contour through that energy's mean. This is NOT the band of the
    # Letter, which is built at fixed S1c from a flat energy spectrum.
    ax.fill(np.concatenate([band[:,1],band[::-1,3]]),
            np.concatenate([band[:,6],band[::-1,8]]),
            color='navy',alpha=.08,lw=0,zorder=1)
    ax.plot(band[:,2],band[:,7],color='navy',lw=2.0,zorder=4,
            label=r'mean response, fixed $E$')
    ax.plot(band[:,1],band[:,6],color='navy',lw=0.9,ls=':',zorder=4)
    ax.plot(band[:,3],band[:,8],color='navy',lw=0.9,ls=':',zorder=4,
            label=r'$\pm1\sigma_{N_e}$ at fixed $E$')
    # displacements
    ax.annotate('',xy=NOM,xytext=EV,arrowprops=dict(arrowstyle='<->',color='darkorange',lw=2.0,mutation_scale=20),zorder=6)
    ax.annotate('',xy=BM,xytext=EV,arrowprops=dict(arrowstyle='<->',color='dimgrey',lw=1.8,mutation_scale=20),zorder=6)
    ax.plot(*EV,'*',ms=17,color='rebeccapurple',mec='none',mew=.5,zorder=8)
    ax.plot(*NOM,'o',ms=7,color='darkorange',mec='k',mew=.5,zorder=8)
    ax.plot(*BM,'s',ms=6,color='dimgrey',mec='k',mew=.5,zorder=8)
    ax.set_xlim(0,680); ax.set_xlabel('$S1c$  [phd]',fontsize=10.5)
    if logy:
        ax.set_yscale('log'); ax.set_ylim(4e2,1.2e5)
        ax.set_ylabel(r'S2c  [phd]  (log)',fontsize=11)
    else:
        ax.set_ylim(0,7.4e4); ax.set_ylabel(r'S2c  [phd]',fontsize=11)


for ax,logy in zip(axes,(False,True)):
    ax.annotate('event',xy=EV,xytext=(-46,-26),textcoords='offset points',fontsize=11,
                color='saddlebrown',arrowprops=dict(arrowstyle='-',color='saddlebrown',lw=.7))
    ax.annotate('nominal, same $N_q$',xy=NOM,xytext=(-118,18),textcoords='offset points',
                fontsize=11,color='darkorange',arrowprops=dict(arrowstyle='-',color='darkorange',lw=.7))
    ax.annotate('band median at fixed $S1c$',xy=BM,xytext=((-150 if logy else -200),(25 if logy else -8)),textcoords='offset points',
                fontsize=11,color='dimgrey',arrowprops=dict(arrowstyle='-',color='dimgrey',lw=.7))
axes[1].set_ylim(4e2,2.2e5)
axes[1].legend(fontsize=9,loc='lower left',framealpha=.95,handlelength=1.8)


plt.tight_layout()
plt.savefig('lz_constant_nq.pdf',bbox_inches='tight')
plt.savefig('lz_constant_nq.png',dpi=185,bbox_inches='tight')

print("event  (%.1f, %.0f)   Nq = %.0f"%(EV[0],EV[1],NQ_EV))
print("nominal(%.1f, %.0f)   same Nq, mean recombination"%NOM)
print("band median at S1c=540 -> (%.1f, %.0f), parent %.1f keVee"%(BM[0],BM[1],E91))
print("slope of the locus, linear plane: %.1f"%(-g2/g1))
