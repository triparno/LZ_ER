"""Digitise the g_A = 0.85 + Exchange curve and its band from the 212Pb figure.

Same three reds as the 214Pb case: solid curve (255,0,0), band fill (254,223,219)
and dashed g_A = 0.85 (249,127,114). The dashed curve lies inside the band and is
treated as band interior so that it never acts as an edge. Where the band is
narrower than the drawn line, the line's own extent is reported instead, which
bounds the band from above; band_resolved flags which regime each row is in.

    python trace212.py [source.png] [output.csv]

The source image is a crop of the published figure and is NOT redistributed
here; see PROVENANCE.md. An overlay <output>_check.png is written beside the
CSV so the trace can be checked against the figure.
"""
import sys, os
_H = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_H, 'Selection_152.png')
DST = sys.argv[2] if len(sys.argv) > 2 else os.path.join(_H, 'pb212_gA085_exchange.csv')
OVL = os.path.splitext(DST)[0] + "_check.png"
from PIL import Image
import numpy as np, csv
a=np.array(Image.open(SRC).convert('RGB')).astype(int)
X0,X1=112.5,878.5          # 0 and 500 keV
Y0,YREF=513.25,163.0       # 0 and 3.0e-3
kx=500./(X1-X0); ky=3.0e-3/(Y0-YREF)
R,G,B=a[:,:,0],a[:,:,1],a[:,:,2]
solid =(R>200)&(G<80)&(B<80)
anyred=(R>200)&(R-G>8)&(R-B>8)
region=anyred|solid
def groups(m,gap=3):
    i=np.where(m)[0]
    if not len(i): return []
    o=[];c=[i[0]]
    for j in i[1:]:
        if j-c[-1]<=gap: c.append(j)
        else: o.append(c); c=[j]
    o.append(c); return o
xs=range(113,993)
cen={};lt={};lb={};prev=None
for x in xs:
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
    if not(-1<=E<=580): continue
    t=min(bt.get(x,lt[x]),lt[x]); b=max(bb.get(x,lb[x]),lb[x])
    res=(lt[x]-t)>1 or (b-lb[x])>1
    rec.append((E,(Y0-cen[x])*ky,(Y0-t)*ky,(Y0-b)*ky,int(res)))
A=np.array(rec); out=[]
for e in np.arange(0,581,1.):
    m=(A[:,0]>=e-0.75)&(A[:,0]<=e+0.75)
    if m.sum(): out.append((e,A[m,1].mean(),A[m,2].mean(),A[m,3].mean(),int(round(A[m,4].mean()))))
with open(DST,'w',newline='') as f:
    w=csv.writer(f); w.writerow(["E_keV","intensity","intensity_upper","intensity_lower","band_resolved"])
    for e,c,u,d,r in out: w.writerow([f"{e:.0f}",f"{c:.6e}",f"{u:.6e}",f"{d:.6e}",r])
print("rows %d (E up to %.0f keV) ; band resolved on %d"%(len(out),out[-1][0],sum(r[4] for r in out)))
print("\n   E      central      upper       lower    res   halfwidth")
for e,c,u,d,r in out:
    if int(e)%50==0 or e in (1,5,15,30):
        print("  %4.0f  %.4e  %.4e  %.4e    %d   %5.2f%%"%(e,c,u,d,r,100*(u-d)/2/c if c>0 else 0))
