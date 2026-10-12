"""Digitise the ER band (navy) from Fig. 2 of the Letter.
Median is the solid navy curve; the 10% and 90% percentiles are the dashed
navy curves above and below it. Output: band_published.csv.

    python trace_band.py [path-to-fig2.png]

The source image is not redistributed here; supply a crop of Fig. 2 of the
Letter. The calibration constants below are the pixel coordinates of 0 and
800 phd on the x axis and of log10 S2c = 5.0 and 3.0 on the y axis, read from
the tick labels; they must be re-measured for a different crop.

band_published.csv in this archive was produced independently of this script,
by the authors digitising the three curves by hand. The two agree on the
median to ~0.005 dex, which is the only check available on either.
"""
from PIL import Image
import numpy as np, csv
import sys
FIG = sys.argv[1] if len(sys.argv) > 1 else 'fig2.png'   # Fig. 2 of the Letter
a=np.array(Image.open(FIG).convert('RGB')).astype(int)
X0,X1,Y50,Y30 = 163.0, 1153.0, 48.5, 777.5      # 0 and 800 phd ; 5.0 and 3.0 dex
kx = 800.0/(X1-X0); ky = 2.0/(Y30-Y50)
R,G,B = a[:,:,0],a[:,:,1],a[:,:,2]
navy = (B>90)&(B<190)&(R<90)&(G<90)             # (0,0,127) and its antialiasing
navy[:40,:]=False
def groups(col,gap=4):
    i=np.where(col)[0]
    if not len(i): return []
    o=[];c=[i[0]]
    for j in i[1:]:
        if j-c[-1]<=gap: c.append(j)
        else: o.append(c); c=[j]
    o.append(c); return o
xs=range(int(X0)+2,int(X1)-1)
raw={}
for x in xs:
    gs=groups(navy[:,x])
    if len(gs)>=1: raw[x]=[float(np.mean(gg)) for gg in gs]
# the median is the middle of three; where only 1-2 groups are visible (dash
# gaps) we defer and interpolate afterwards
tri={x:v for x,v in raw.items() if len(v)==3}
print("columns with all three curves resolved: %d of %d"%(len(tri),len(list(xs))))
xt=np.array(sorted(tri))
hi=np.array([tri[x][0] for x in xt])   # smallest y = highest S2c = 90th pct
md=np.array([tri[x][1] for x in xt])
lo=np.array([tri[x][2] for x in xt])
grid=np.arange(60,801,20.0)
px=lambda E: X0+E/kx
out=[]
for S in grid:
    xp=px(S)
    if xp<xt[0] or xp>xt[-1]: continue
    q90=(np.interp(xp,xt,hi)-Y50)*ky*-1.0+5.0
    q50=(np.interp(xp,xt,md)-Y50)*ky*-1.0+5.0
    q10=(np.interp(xp,xt,lo)-Y50)*ky*-1.0+5.0
    out.append((S,q10,q50,q90))
with open('band_published_trace.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['S1c_phd','log10S2c_10','log10S2c_50','log10S2c_90'])
    for r in out: w.writerow(['%.0f'%r[0]]+['%.4f'%v for v in r[1:]])
print("rows written:",len(out))
print("\n   S1c    10%     50%     90%")
for S,a_,b_,c_ in out:
    if int(S)%100==0: print("  %5.0f  %6.3f  %6.3f  %6.3f"%(S,a_,b_,c_))
