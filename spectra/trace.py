"""Digitise the g_A = 0.85 + Exchange curve and its uncertainty band from the
214Pb beta-spectrum figure.

Three reds are present: the solid curve (255,0,0), the band fill (254,223,219)
and the dashed g_A = 0.85 curve (249,127,114). The dashed curve runs through the
band and would otherwise split the fill into disconnected pieces, so it is
treated as band interior and plays no part in locating the edges.

Where the band is narrower than the drawn line its edges cannot be separated
from the line itself; there the line's own extent is reported, which bounds the
true band from above. band_resolved flags which regime each row is in.

    python trace.py [source.png] [output.csv]

The source image is a crop of the published figure and is NOT redistributed
here; see PROVENANCE.md. An overlay <output>_check.png is written beside the
CSV so the trace can be checked against the figure.
"""
import sys, os
_H = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_H, 'Selection_151.png')
DST = sys.argv[2] if len(sys.argv) > 2 else os.path.join(_H, 'pb214_gA085_exchange.csv')
OVL = os.path.splitext(DST)[0] + "_check.png"
from PIL import Image
import numpy as np, csv
a=np.array(Image.open(SRC).convert('RGB')).astype(int)
X0,X1,Y0,Y1=162.5,975.0,509.0,122.5          # 0 and 1000 keV ; 0 and 1.50e-3
kx=1000./(X1-X0); ky=1.50e-3/(Y0-Y1)
R,G,B=a[:,:,0],a[:,:,1],a[:,:,2]
solid =(R>200)&(G<80)&(B<80)                  # the curve we want
anyred=(R>200)&(R-G>8)&(R-B>8)                # fill + dashed + both antialiasings
region=anyred|solid                           # dashed treated as band interior
def groups(m,gap=3):
    i=np.where(m)[0]
    if not len(i): return []
    o=[];c=[i[0]]
    for j in i[1:]:
        if j-c[-1]<=gap: c.append(j)
        else: o.append(c); c=[j]
    o.append(c); return o
xs=range(163,998)
cen={};lt={};lb={};prev=None
for x in xs:                                   # continuity rejects the legend line
    g=groups(solid[:,x])
    if not g: continue
    k=max(g,key=len) if prev is None else min(g,key=lambda z:abs(np.mean(z)-prev))
    cen[x]=np.mean(k); lt[x]=min(k); lb[x]=max(k); prev=cen[x]
bt={};bb={}
for x in xs:
    if x not in cen: continue
    g=[z for z in groups(region[:,x],gap=4) if z[0]-5<=lb[x] and z[-1]+5>=lt[x]]
    if not g: continue
    k=min(g,key=lambda z:abs(np.mean(z)-cen[x]))
    if abs(np.mean(k)-cen[x])>60: continue
    bt[x]=min(k); bb[x]=max(k)
rec=[]
for x in xs:
    if x not in cen: continue
    E=(x-X0)*kx
    if not(-1<=E<=1002): continue
    t=min(bt.get(x,lt[x]),lt[x]); b=max(bb.get(x,lb[x]),lb[x])
    res=(lt[x]-t)>1 or (b-lb[x])>1
    rec.append((E,(Y0-cen[x])*ky,(Y0-t)*ky,(Y0-b)*ky,int(res)))
A=np.array(rec); out=[]
for e in np.arange(0,1001,2.):
    m=(A[:,0]>=e-1.5)&(A[:,0]<=e+1.5)
    if m.sum(): out.append((e,A[m,1].mean(),A[m,2].mean(),A[m,3].mean(),int(round(A[m,4].mean()))))
with open(DST,'w',newline='') as f:
    w=csv.writer(f); w.writerow(["E_keV","intensity","intensity_upper","intensity_lower","band_resolved"])
    for e,c,u,d,r in out: w.writerow([f"{e:.0f}",f"{c:.6e}",f"{u:.6e}",f"{d:.6e}",r])
print("rows %d ; band resolved on %d"%(len(out),sum(r[4] for r in out)))
print("\n   E      central      upper       lower    res   halfwidth")
for e,c,u,d,r in out:
    if int(e)%100==0 or e in (2,10,30,60,150):
        print("  %4.0f  %.4e  %.4e  %.4e    %d   %5.2f%%"%(e,c,u,d,r,100*(u-d)/2/c if c>0 else 0))
