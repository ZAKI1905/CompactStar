"""Validation-only cold neutral equilibrium reconstructions. Not an owner/bound."""
import csv,math,sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
S=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(S/'tests/analysis'))
from chemical_trackr_budget import Model,matrix,zsolve,integrate
MASS=np.array([939.56542051999998,938.27208815999995,.51099894999999995,105.65837550000001])
MN,MP,ME,MM=MASS;HC=197.32698045930249;KB=8.61733326214518e-11;CONV=1.602176634e33
TON=7.3567289037328318e-9;MUM=np.hypot(MP,np.sqrt(MM*MM-ME*ME))+MM

def mu_pe(n):
 p=HC*np.cbrt(3*np.pi**2*n);return np.hypot(MP,p)+np.hypot(ME,p)

def momenta(mu):
 mu=np.asarray(mu);out=np.zeros(mu.shape+(4,));nomu=mu<MUM
 # Algebraic two-body inversion, independent of review's 200-step bisection.
 v=mu[nomu];p2=((v-MP-ME)*(v+MP+ME)*(v-MP+ME)*(v+MP-ME))/(4*v*v)
 out[nomu,1]=out[nomu,2]=np.sqrt(np.maximum(p2,0))
 if np.any(~nomu):
  v=mu[~nomu];lo=np.full(v.shape,MM);hi=v-MP
  for _ in range(52):
   mid=(lo+hi)/2;pe=np.sqrt(np.maximum((mid-ME)*(mid+ME),0));pm=np.sqrt(np.maximum((mid-MM)*(mid+MM),0));pp=np.cbrt(pe**3+pm**3)
   pos=np.hypot(MP,pp)+mid>v;hi=np.where(pos,mid,hi);lo=np.where(pos,lo,mid)
  lam=(lo+hi)/2;pe=np.sqrt((lam-ME)*(lam+ME));pm=np.sqrt(np.maximum((lam-MM)*(lam+MM),0))
  out[~nomu,1]=np.cbrt(pe**3+pm**3);out[~nomu,2]=pe;out[~nomu,3]=pm
 # Match review's neutron formula; stable factorized inverse checked separately.
 out[...,0]=np.sqrt(np.maximum(mu*mu-MN*MN,0));return out

def gas_energy(p):
 # Independent momentum quadrature; stable even at small p/m.
 x,w=leggauss(24);t=(x+1)/2;w=w/2
 return np.sum(p[...,None]**3/(np.pi**2*HC**3)*w*t*t*np.hypot(MASS[:,None],p[...,None]*t),axis=(-2,-1))

def edges(d):
 r,m,nu,nb=d[:,:4].T;k=int(np.searchsorted(-nb,-TON));a=k-1
 assert 0<a<k<len(r)-1 and np.all(np.diff(nu)>0)
 mi=np.exp(nu[k])*mu_pe(nb[k]);alt=np.exp(nu[k+1])*mu_pe(nb[k+1]);no=np.log(mi/MN)
 assert nu[a]<=no<=nu[k], 'physical onset outside bracket'
 rn=float(np.interp(no,nu,r));rp=r[a]+(nb[a]-TON)/(nb[a]-nb[k])*(r[k]-r[a])
 return dict(k=k,r_nu=rn,r_nb=float(rp),mu_inf=float(mi),mu_inf_altref=float(alt),r_nu_altref=float(np.interp(np.log(alt/MN),nu,r)),R=float(r[-1]))

def cache_factor(T):
 grid=np.exp((1-np.arange(160)/159)*math.log(1e-5)+np.arange(160)/159*math.log(1.))
 v=np.asarray(T)*KB;j=np.searchsorted(grid,v)-1;w=np.log(v/grid[j])/np.log(grid[j+1]/grid[j]);return ((1-w)*grid[j]+w*grid[j+1])/v

def integrate_fields(d,e,order=24,radial_composition=False):
 r,m,nu,nb=d[:,:4].T;ons=[e['r_nu'],float(np.interp(np.log(e['mu_inf']/MUM),nu,r))]
 cuts=np.unique(np.r_[0.,r,ons,e['r_nb']]);x,w=leggauss(order);t=(x+1)/2;w=w/2
 total=np.zeros(7,dtype=np.longdouble)
 for j in range(0,len(cuts)-1,1024):
  a=cuts[j:min(j+1024,len(cuts)-1),None];b=cuts[j+1:min(j+1024,len(cuts)-1)+1,None]
  rr=a+(b-a)*t;ww=np.broadcast_to((b-a)*w,rr.shape).copy()
  for o in ons:
   left=b[:,0]==o;right=a[:,0]==o
   rr[left]=b[left]-(b[left]-a[left])*t*t;ww[left]=(b[left]-a[left])*2*t*w
   rr[right]=a[right]+(b[right]-a[right])*t*t;ww[right]=(b[right]-a[right])*2*t*w
  mm=np.interp(rr,r,m);mm=np.where(rr<r[0],m[0]*(rr/r[0])**3,mm);nv=np.interp(rr,r,nu)
  vol=ww*4*np.pi*rr**2/np.sqrt(1-2*mm/rr)
  l=vol*np.exp(-6*nv)*1e-36
  total[0]+=np.sum(np.where(rr<e['r_nu'],l,0),dtype=np.longdouble)
  total[1]+=np.sum(np.where(rr<e['r_nb'],l,0),dtype=np.longdouble)
  p=momenta(e['mu_inf']*np.exp(-nv));s=np.hypot(MASS,p)*p/(np.pi**2*HC**3)
  total[2]+=np.sum(vol*np.exp(-nv)*np.sum(s,axis=-1)*np.pi**2/3*1e15*KB**2*CONV,dtype=np.longdouble)
  den=s[...,1]+s[...,2]+s[...,3]
  response=np.stack([s[...,0],s[...,2]*(s[...,1]+s[...,3])/den,-s[...,2]*s[...,3]/den,s[...,3]*(s[...,1]+s[...,2])/den],axis=-1)
  total[3:]+=np.sum((vol*np.exp(-nv))[...,None]*response,axis=(0,1),dtype=np.longdouble)*1e54
 return {'L_nu':float(total[0]),'L_nb':float(total[1]),'C_nu':float(total[2]),'G_nu':matrix(np.array(total[3:],float)).tolist(),'Z_nu':zsolve(matrix(np.array(total[3:],float))).tolist()}

def analyze(path,hp,order=24,extra=False):
 d=np.loadtxt(path/'profile/profile.tsv',skiprows=1);r,m,nu,nb=d[:,:4].T;e=edges(d);e.update(integrate_fields(d,e,order))
 rows=list(csv.DictReader((path/'thermal.tsv').open(),delimiter='\t'));T=np.array([float(v['Tinf_K']) for v in rows]);C=np.array([float(v['Cstar']) for v in rows]);F=cache_factor(T)
 nodalp=HC*np.cbrt(3*np.pi**2*d[:,6:10]);nodals=np.sum(nodalp*np.hypot(MASS,nodalp),axis=1)/(3*HC**3)
 f=4*np.pi*r*r/np.sqrt(1-2*m/r)*np.exp(-nu)*nodals*1e15*KB**2*CONV
 ct=float(np.trapezoid(f,r))
 mui=np.where(nb<=TON,mu_pe(nb),np.hypot(MN,nodalp[:,0]))*np.exp(nu)
 ent=(d[:,10]+d[:,11])/nb*np.exp(nu)
 rhs=(m+4*np.pi*r**3*d[:,4])/(r*(r-2*m))
 # Node secants vs midpoint geometric RHS; a diagnostic, not independent accuracy.
 sec=np.diff(nu)/np.diff(r);midrhs=(rhs[:-1]+rhs[1:])/2
 # Actual production selection uses profile energy column21, not provider column10.
 geom=1.3235754613110692e-6;epsb=hp['eps_b'];ib=np.where(d[:,21]>=epsb*geom)[0][-1]
 nub=np.log(e['mu_inf']/hp['mu_b']);rb=float(np.interp(nub,nu,r));nulog=float(np.interp(np.log(epsb*geom),np.log(d[::-1,21]),nu[::-1]));rlog=float(np.interp(nulog,nu,r))
 # independent physical P threshold applied to TOV P, no mu_inf identity
 rpressure=float(np.interp(np.log(hp['P_on']*geom),np.log(d[::-1,4]),r[::-1]));rbpressure=float(np.interp(np.log(hp['P_b']*geom),np.log(d[::-1,4]),r[::-1]))
 tbnode=T*np.exp(-nu[ib]);tb=T*np.exp(-nub);g=np.array([float(v['g14']) for v in rows]);ts=(1e24*g*(1.81*tb/1e8)**2.42)**.25;tsp=(1e24*g*(1.81*tbnode/1e8)**2.42)**.25
 actualts=np.array([float(v['Tsurface']) for v in rows])
 re=e['r_nu'];me=np.interp(re,r,m);ne=np.interp(re,r,nu)
 e.update(C_trap=ct,C_production_slope=float(C[0]/T[0]),cache_factor=F.tolist(),thermal=[dict(T=float(t),C_prod=float(c),C_nu=float(t*e['C_nu']),C_decache=float(c/f),Tb_node=float(tn),Tb_nu=float(tb_),Ts_prod=float(ts_),Ts_nu=float(tsn),Ts_inf_nu=float(tsn*np.exp(nu[-1]))) for t,c,f,tn,tb_,ts_,tsn in zip(T,C,F,tbnode,tb,actualts,ts)],cache_law_residual=float(np.max(np.abs(C/(T*F)/(C[0]/T[0]/F[0])-1))),nodal_vs_decache=float(ct/(C[0]/T[0]/F[0])-1),mu_constancy=float(np.max(np.abs(mui/e['mu_inf']-1))),enthalpy_mu_disagreement=float(np.max(np.abs(ent/mui-1))),metric_secant_max=float(np.max(np.abs(sec[3:]/midrhs[3:]-1))),ib=int(ib),nu_node=float(nu[ib]),nu_b=float(nub),Tb_factor=float(np.exp(-nub)),Tb_node_factor=float(np.exp(-nu[ib])),rb=rb,rb_logeps=rlog,rb_pressure=rbpressure,r_on_pressure=rpressure,Ts_node_reproduction=float(np.max(np.abs(actualts/tsp-1))),edge_sensitivity=float(1e-36*np.exp(-6*ne)*4*np.pi*re**2/np.sqrt(1-2*me/re)/e['L_nu']))
 star=next(csv.DictReader((path/'star.tsv').open(),delimiter='\t'));e['star']={k:float(v) for k,v in star.items()}
 if extra:
  model=Model(path/'profile');gg=matrix(integrate(d,model,32)[0]);e['G_prod']=gg.tolist();e['Z_prod']=zsolve(gg).tolist()
  e['order48']=integrate_fields(d,e,48)
  # composition sampled on nodes, then interpolated separately from spatial rule
  x,w=leggauss(24);t=(x+1)/2;w=w/2;val=np.longdouble(0);val_nb=np.longdouble(0)
  for j in range(0,len(r)-1,1024):
   a=r[j:min(j+1024,len(r)-1),None];b=r[j+1:min(j+1024,len(r)-1)+1,None];rr=a+(b-a)*t
   mm=np.interp(rr,r,m);nv=np.interp(rr,r,nu);nn=np.interp(rr,r,nb)
   vol=(b-a)*w*4*np.pi*rr**2/np.sqrt(1-2*mm/rr)*np.exp(-nv)*1e15*KB**2*CONV
   ns=np.stack([np.interp(rr,r,d[:,k]) for k in range(6,10)],axis=-1);p=HC*np.cbrt(3*np.pi**2*ns)
   val+=np.sum(vol*np.sum(p*np.hypot(MASS,p),axis=-1)/(3*HC**3),dtype=np.longdouble)
   ns=model.state(nn);p=HC*np.cbrt(3*np.pi**2*ns)
   val_nb+=np.sum(vol*np.sum(p*np.hypot(MASS,p),axis=-1)/(3*HC**3),dtype=np.longdouble)
  e['C_species_linear']=float(val);e['C_nb_equilibrium']=float(val_nb)
 return e
