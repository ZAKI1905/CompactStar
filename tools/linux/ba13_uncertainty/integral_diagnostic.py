#!/usr/bin/env python3
"""Independent onset-split fixed-background Cv, event and benchmark L integrals.
No fitted error percentages: quadrature spreads are diagnostics, never global bounds.
"""
import json,sys,time,math
from pathlib import Path
import numpy as np
S=Path(__file__).resolve().parents[3];R=S.parents[1]/'ba13-uncertainty'
sys.path.insert(0,str(S/'tests/analysis'))
from chemical_trackr_budget import Model
KB=8.61733326214518e-11;CONV=1.602176634e33

def integrate(p,order,bisect=False):
 d=np.loadtxt(p/'profile.tsv',skiprows=1);model=Model(p);r,m,nu,nb=d[:,:4].T
 ons=np.interp(model.onsets,nb[::-1],r[::-1]);cuts=np.unique(np.r_[0,r,ons])
 if bisect:cuts=np.sort(np.r_[cuts,(cuts[:-1]+cuts[1:])/2])
 x,w=np.polynomial.legendre.leggauss(order);t=(x+1)/2;w=w/2
 px,pw=np.polynomial.legendre.leggauss(12);pt=(px+1)/2;pw=pw/2
 sums=[[] for _ in range(8)]
 for k in range(0,len(cuts)-1,256):
  a=cuts[k:min(k+256,len(cuts)-1),None];b=cuts[k+1:min(k+256,len(cuts)-1)+1,None]
  rr=a+(b-a)*t;ww=np.broadcast_to((b-a)*w,rr.shape).copy()
  for onset in ons:
   left=b[:,0]==onset;right=a[:,0]==onset
   rr[left]=b[left]-(b[left]-a[left])*t*t;ww[left]=(b[left]-a[left])*2*t*w
   rr[right]=a[right]+(b[right]-a[right])*t*t;ww[right]=(b[right]-a[right])*2*t*w
  mm=np.interp(rr,r,m);mm=np.where(rr<r[0],m[0]*(rr/r[0])**3,mm)
  nv=np.interp(rr,r,nu);bn=np.interp(rr,r,nb);ns=model.state(bn)
  pf=model.hc*np.cbrt(3*np.pi**2*ns);cv=np.sum(pf*np.hypot(pf,model.mass),axis=-1)/(3*model.hc**3)
  nodal=np.stack([np.interp(rr,r,d[:,j]) for j in range(6,10)],axis=-1)
  pn=model.hc*np.cbrt(3*np.pi**2*nodal);cvn=np.sum(pn*np.hypot(pn,model.mass),axis=-1)/(3*model.hc**3)
  volume=ww*4*np.pi*rr**2/np.sqrt(1-2*mm/rr);lapse=np.exp(nv)
  measure=volume*lapse*ns[...,0]*1e54;mu=np.hypot(model.mass[0],pf[...,0])
  # Historical P1 uses its explicitly hard-coded momentum constant. Independent momentum integration.
  p1pf=197.3269804*np.cbrt(3*np.pi**2*ns[...,0])
  average=sum(3*pp*pp*weight*np.hypot(939.56542052,p1pf*pp) for pp,weight in zip(pt,pw))
  vals=[volume/lapse*cv*1e15*KB**2*CONV,volume/lapse*cvn*1e15*KB**2*CONV,
   measure,measure*lapse*mu,measure*lapse*average,
   volume*np.exp(-6*nv)*(rr<ons[0])*1e-36,
   volume*np.exp(-6*nv)*(rr<ons[1])*2e-36,
   volume/lapse*cv*1e15*KB**2*CONV*(np.abs(bn-model.onsets[0])<model.onsets[0]*.01)]
  for out,v in zip(sums,vals):out.append(float(np.sum(v,dtype=np.longdouble)))
 v=list(map(math.fsum,sums))
 return dict(Cstar_slope=v[0],nodal_density_Cstar_slope=v[1],neutron_measure=v[2],mu_B=v[3]/v[2],P0=v[3]/v[2],P1=v[4]/v[2],P2=0.,Ltilde_Me=v[5],Ltilde_Mmu=v[6],neutron_onset_Cstar_slope=v[7],order=order,bisect=bisect)
records=[]
for i in (0,10,20):
 for case in ('nominal','r160000','r320000'):
  orders=[(16,False),(32,False),(64,False),(96,False),(32,True)] if case=='nominal' else [(16,False),(32,False)]
  for order,bisect in orders:
   start=time.monotonic();v=integrate(R/'runs'/f'{case}-{i}'/'profile',order,bisect)
   records.append(dict(index=i,case=case,seconds=time.monotonic()-start,**v));print(i,case,order,bisect,records[-1]['seconds'],flush=True)
  # An interrupted experiment retains each completed target/case without publishing COMPLETE.
  (R/'integral-diagnostic.partial.json').write_text(json.dumps(records,indent=2)+'\n')
(R/'integral-diagnostic.partial.json').rename(R/'integral-diagnostic.json')
