#!/usr/bin/env python3
"""Direct event integrals and support-coordinate mapping on retained profiles."""
import csv,hashlib,json,os,pathlib,resource,sys,time
import numpy as np
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1];O=R/'ba13-headroom';raw=R/'ba13-recovery/run1/raw'
sys.path.insert(0,str(S/'tools/linux/ba13_onset'))
from equilibrium import momenta,HC
if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
start=time.process_time();wall=time.time();results=[]
metrics=list(csv.DictReader((raw/'profile_metrics.tsv').open(),delimiter='\t'))
for i in range(21):
 p=raw/'chemical-profiles'/f'target-{i}'/'profile.tsv';d=np.loadtxt(p,skiprows=1);r,m,nu,nb=d[:,:4].T
 ref=json.loads((O/f'profile-v2-{i}.json').read_text());edge=ref['r_nu'];mui=ref['mu_inf']
 original=[v for v in metrics if int(v['index'])==i];q=np.array([float(v['q']) for v in original]);rr=q*edge
 mm=np.interp(rr,r,m);mm=np.where(rr<r[0],m[0]*(rr/r[0])**3,mm);nv=np.interp(rr,r,nu)
 ratio=np.divide(2*mm,rr,out=np.zeros_like(rr),where=rr>0);lam=-.5*np.log1p(-ratio)
 mapped=np.r_[rr,mm,nv,lam]
 integrals={}
 for order in ([16,32] if i in [0,10,20] else [16]):
  x,w=np.polynomial.legendre.leggauss(order);t=(x+1)/2;w=w/2
  z,zw=np.polynomial.legendre.leggauss(24);z=(z+1)/2;zw=zw/2
  cuts=np.unique(np.r_[0,r[r<edge],edge]);total=np.zeros(3,dtype=np.longdouble)
  for j in range(0,len(cuts)-1,512):
   a=cuts[j:min(j+512,len(cuts)-1),None];b=cuts[j+1:min(j+512,len(cuts)-1)+1,None];rr=a+(b-a)*t;ww=np.broadcast_to((b-a)*w,rr.shape).copy()
   mask=b[:,0]==edge;rr[mask]=b[mask]-(b[mask]-a[mask])*t*t;ww[mask]=(b[mask]-a[mask])*2*t*w
   mm=np.interp(rr,r,m);mm=np.where(rr<r[0],m[0]*(rr/r[0])**3,mm);nv=np.interp(rr,r,nu)
   mom=momenta(mui*np.exp(-nv));nn=(mom[...,0]/HC)**3/(3*np.pi**2)
   vol=ww*4*np.pi*rr*rr/np.sqrt(1-2*mm/rr);measure=vol*np.exp(nv)*nn
   pf=197.3269804*np.cbrt(3*np.pi**2*nn);average=np.zeros_like(rr)
   for zz,wz in zip(z,zw):average+=3*wz*zz*zz*np.hypot(939.56542052,pf*zz)
   total[0]+=np.sum(measure,dtype=np.longdouble)
   total[1]+=np.sum(measure*np.exp(nv)*np.hypot(939.56542052,mom[...,0]),dtype=np.longdouble)
   total[2]+=np.sum(measure*np.exp(nv)*average,dtype=np.longdouble)
  integrals[str(order)]={'P0':float(total[1]/total[0]),'P1':float(total[2]/total[0])}
 results.append(dict(index=i,integrals=integrals,mapped_metric=mapped.tolist(),profile_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
result=dict(records=results,cpu_seconds=time.process_time()-start,wall_seconds=time.time()-wall,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,job=os.environ['SLURM_JOB_ID'],source_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest())
text=json.dumps(result,allow_nan=False)
with (O/'supplement.json').open('x') as f:f.write(text+'\n')
print(result['cpu_seconds'],result['wall_seconds'])
