"""Digitise the 'Exchange - 1st Forb.' dashed red curve from the 85Kr figure.

Two dashed reds are present: the target (255,0,0) and 'No Exchange' (249,127,114).
Only pure red is selected. The figure has two panels; only the upper one
(y 68-409) is used. The curve is dashed, so gaps between dashes are filled by
interpolation and flagged. There is no shaded band in this figure.

    python trace85.py [source.png] [output.csv]

The source image is a crop of the published figure and is NOT redistributed
here; see PROVENANCE.md. An overlay <output>_check.png is written beside the
CSV so the trace can be checked against the figure.
"""
import sys, os
_H = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_H, 'Selection_153.png')
DST = sys.argv[2] if len(sys.argv) > 2 else os.path.join(_H, 'kr85_exchange_1stforb.csv')
OVL = os.path.splitext(DST)[0] + "_check.png"
from PIL import Image
import numpy as np, csv
a=np.array(Image.open(SRC).convert('RGB')).astype(int)
X0,X1=131.5,898.0           # 0 and 600 keV
Y0,YREF=410.0,67.0          # 0 and 3.0e-3
kx=600./(X1-X0); ky=3.0e-3/(Y0-YREF)
R,G,B=a[:,:,0],a[:,:,1],a[:,:,2]
tgt=(R>200)&(G<80)&(B<80)
tgt[:68,:]=False; tgt[410:,:]=False              # upper panel only
def groups(m,gap=3):
    i=np.where(m)[0]
    if not len(i): return []
    o=[];c=[i[0]]
    for j in i[1:]:
        if j-c[-1]<=gap: c.append(j)
        else: o.append(c); c=[j]
    o.append(c); return o
xs=list(range(133,1008))
raw={}; prev=None
for x in xs:
    g=groups(tgt[:,x])
    if not g: continue
    k=max(g,key=len) if prev is None else min(g,key=lambda z:abs(np.mean(z)-prev))
    if prev is not None and abs(np.mean(k)-prev)>18: continue   # legend / stray
    raw[x]=np.mean(k); prev=raw[x]
have=sorted(raw); yv=np.interp(xs,have,[raw[h] for h in have])
rec=[]
for x,y in zip(xs,yv):
    E=(x-X0)*kx
    if not(-1<=E<=700): continue
    rec.append((E,(Y0-y)*ky,int(x in raw)))
A=np.array(rec); out=[]
for e in np.arange(0,691,1.):
    m=(A[:,0]>=e-0.6)&(A[:,0]<=e+0.6)
    if m.sum(): out.append((e,A[m,1].mean(),int(round(A[m,2].mean()))))
with open(DST,'w',newline='') as f:
    w=csv.writer(f); w.writerow(["E_keV","intensity","on_dash"])
    for e,c,d in out: w.writerow([f"{e:.0f}",f"{c:.6e}",d])
print("rows %d (E up to %.0f keV); %d on a dash, %d interpolated"
      %(len(out),out[-1][0],sum(r[2] for r in out),len(out)-sum(r[2] for r in out)))
print("\n   E      intensity    on_dash")
for e,c,d in out:
    if int(e)%50==0 or e in (1,3,6,10,20):
        print("  %4.0f   %.4e     %d"%(e,c,d))
