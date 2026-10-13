"""70-digit EOS oracle using exact binary input constants; no production EOS."""
import json,sys
from pathlib import Path
import mpmath as mp
mp.mp.dps=70
MN,MP,ME,MM=map(lambda x:mp.mpf(float(x)),[939.56542051999998,938.27208815999995,.51099894999999995,105.65837550000001]);HC=mp.mpf(float(197.32698045930249))
def pf(mu,m):return mp.sqrt(max((mu-m)*(mu+m),mp.mpf(0)))
def gas(p,m):
 mu=mp.sqrt(m*m+p*p);a=mp.asinh(p/m)
 return ((mu*p*(2*p*p-3*m*m)+3*m**4*a)/(24*mp.pi**2*HC**3),(mu*p*(2*p*p+m*m)-m**4*a)/(8*mp.pi**2*HC**3))
def momenta(mu):
 p=mp.sqrt(max((mu-MP-ME)*(mu+MP+ME)*(mu-MP+ME)*(mu+MP-ME)/(4*mu*mu),0))
 return [pf(mu,MN),p,p,mp.mpf(0)] # all audited onset/base states below muon onset

def atmu(mu):
 ps=momenta(mu);g=[gas(p,m) for p,m in zip(ps,[MN,MP,ME,MM])];n=(ps[0]**3+ps[1]**3)/(3*mp.pi**2*HC**3)
 return n,sum(v[0] for v in g),sum(v[1] for v in g)
def bisect(f,lo,hi,iters=250):
 for _ in range(iters):
  mid=(lo+hi)/2
  if f(mid)>0:hi=mid
  else:lo=mid
 return (lo+hi)/2

def atn(n):
 mu=bisect(lambda mu:atmu_density(mu)-n,MP+ME,MN+1)
 return atmu(mu),mu

def atmu_density(mu):
 p=momenta(mu);return (p[0]**3+p[1]**3)/(3*mp.pi**2*HC**3)

def constants():
 n,p,e=atmu(MN);eb=mp.mpf('1e10')*mp.mpf('2.99792458e10')**2/mp.mpf('1.602176634e33')
 mub=bisect(lambda u:atmu(u)[2]-eb,MN,MN+1);nb,pb,_=atmu(mub)
 points=[]
 for u in [MP+ME+mp.mpf('.1'),MN,MN+mp.mpf('.0001'),mub]:
  n1,p1,e1=atmu(u);der=mp.diff(lambda v:atmu(v)[1],u)
  points.append(dict(mu=str(u),enthalpy_relative=str((e1+p1)/(n1*u)-1),dP_n_relative=str(der/n1-1)))
 return dict(onset_n=float(n),P_on=float(p),eps_on=float(e),mu_b=float(mub),P_b=float(pb),n_b=float(nb),eps_b=float(eb),precision=70,thermodynamics=points)

def audit(root):
 import numpy as np
 results=[];old=[(7.3567289037320691e-09,7.3567289037320765e-09),(7.3567289037320832e-09,7.3567289037320889e-09),(7.3567289037336159e-09,7.3567289037336225e-09),(7.3567289037337119e-09,7.3567289037337201e-09)]
 pairs=[('historical32768',a,b) for a,b in old]
 for file in sorted((root/'tables').glob('*/capped.tsv')):
  d=np.loadtxt(file,skiprows=1);k=int(np.searchsorted(d[:,2],7.3567289037328318e-9))
  for j in [k-2,k-1,k,k+1]:pairs.append((file.parent.name,float(d[j,2]),float(d[j+1,2])))
 for name,a,b in pairs:
  va,ua=atn(mp.mpf(a));vb,ub=atn(mp.mpf(b));pa,ea=va[1:];pb,eb=vb[1:]
  results.append(dict(case=name,a=a,b=b,high_precision_order=bool(pb>pa and eb>ea),rounded_order=bool(float(pb)>float(pa) and float(eb)>float(ea)),pressure_ulp_gap=float((pb-pa)/mp.mpf(float(np.spacing(float(pa))))),energy_ulp_gap=float((eb-ea)/mp.mpf(float(np.spacing(float(ea))))),provider_pf_equal=bool(float(HC)*np.cbrt(3*np.pi**2*a)==float(HC)*np.cbrt(3*np.pi**2*b))))
 return results
if __name__=='__main__':
 root=Path(sys.argv[1]);(root/'high-precision.json').write_text(json.dumps(dict(constants=constants(),pairs=audit(root)),indent=2)+'\n')
