#!/usr/bin/env python3
"""Read-only validation of actual 21 nominal profiles; never regenerates a star."""
import csv,hashlib,json,os,pathlib,resource,sys,time
import numpy as np
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1];O=R/'ba13-headroom';raw=R/'ba13-recovery/run1/raw'
sys.path.insert(0,str(S/'tools/linux/ba13_onset'))
from equilibrium import edges,integrate_fields,momenta,cache_factor,MASS,HC,KB,CONV,MN,mu_pe,TON
sys.path.insert(0,str(S/'tests/analysis'))
from chemical_trackr_budget import Model,matrix,integrate,zsolve
from enclosure import inverse_majorant,paired_inverse,L

def table(name):return [{k:float(v) for k,v in r.items()} for r in csv.DictReader((raw/name).open(),delimiter='\t')]
def compute(i):
 start=time.process_time();wall=time.time();p=raw/'chemical-profiles'/f'target-{i}';d=np.loadtxt(p/'profile.tsv',skiprows=1);model=Model(p)
 r,m,nu,nb=d[:,:4].T;e=edges(d);e.update(integrate_fields(d,e,24))
 c=json.loads((p/'certificate.json').read_text());coef=table('coefficients.tsv')[i];star=table('target_stars.tsv')[i];thermal=[v for v in table('thermal_surface.tsv') if v['index']==i]
 # Full independently assembled fixed-background chemical response, same original cut.
 gp=matrix(integrate(d,model,16)[0]);z=zsolve(gp)
 center=abs(matrix(model.response(np.array([nb[0]]))[0]))*4*np.pi/3*r[0]**3*np.exp(-nu[0])*1e54
 errors={k:np.array(c[k]) for k in ['background_error','center_error','tail_error']}
 # Retain inherited background as an input hypothesis, not a newly validated bound.
 _,maj,rem,rho=inverse_majorant(gp,sum(errors.values()))
 e.update(G_fixed_profile=gp.tolist(),Z_fixed_profile=z.tolist(),G_background_discrepancy=(np.array(e['G_nu'])-gp).tolist(),
   center_current=center.tolist(),center_inherited=c['center_error'],tail=c['tail'],tail_error=c['tail_error'],background_inherited=c['background_error'],
   inverse_conditional_majorant=(abs(L).T@maj@abs(L)).tolist(),inverse_nonlinear_remainder=(abs(L).T@rem@abs(L)).tolist(),inverse_rho=rho,
   coefficient_Z=[[coef[f'Z{j}{k}'] for k in range(2)] for j in range(2)])
 if i in (0,10,20):
  e['order48']=integrate_fields(d,e,48);g32=matrix(integrate(d,model,32)[0]);e['G_order32']=g32.tolist();e['Z_order32']=zsolve(g32).tolist()
 # Integral comparison for ALL numbers: production nodal values vs neutral mu reconstruction.
 x,w=np.polynomial.legendre.leggauss(16);t=(x+1)/2;w=w/2
 nums=np.zeros(4,dtype=np.longdouble);nodal=np.zeros(4,dtype=np.longdouble);baryon=np.longdouble(0)
 # Every cell split at physical onsets; high-side square map removes sqrt singularities.
 from equilibrium import MUM
 ons=[e['r_nu'],float(np.interp(np.log(e['mu_inf']/MUM),nu,r))]
 cuts=np.unique(np.r_[0,r,ons]);muon_L=np.longdouble(0)
 for j in range(0,len(cuts)-1,512):
  a=cuts[j:min(j+512,len(cuts)-1),None];b=cuts[j+1:min(j+512,len(cuts)-1)+1,None];rr=a+(b-a)*t;ww=np.broadcast_to((b-a)*w,rr.shape).copy()
  for edge in ons:
   left=b[:,0]==edge;right=a[:,0]==edge
   rr[left]=b[left]-(b[left]-a[left])*t*t;ww[left]=(b[left]-a[left])*2*t*w
   rr[right]=a[right]+(b[right]-a[right])*t*t;ww[right]=(b[right]-a[right])*2*t*w
  mm=np.interp(rr,r,m);mm=np.where(rr<r[0],m[0]*(rr/r[0])**3,mm);nv=np.interp(rr,r,nu)
  vol=ww*4*np.pi*rr*rr/np.sqrt(1-2*mm/rr)
  ns=(momenta(e['mu_inf']*np.exp(-nv))/HC)**3/(3*np.pi**2)
  ns_node=np.stack([np.interp(rr,r,d[:,k]) for k in range(6,10)],axis=-1)
  nums+=np.sum(vol[...,None]*ns,axis=(0,1),dtype=np.longdouble)*1e54
  nodal+=np.sum(vol[...,None]*ns_node,axis=(0,1),dtype=np.longdouble)*1e54
  baryon+=np.sum(vol*np.interp(rr,r,nb),dtype=np.longdouble)*1e54
  muon_L+=np.sum(np.where(rr<ons[1],vol*np.exp(-6*nv)*2e-36,0),dtype=np.longdouble)
 T=np.array([v['Tinf_K'] for v in thermal]);C=np.array([v['Cstar_erg_K'] for v in thermal]);F=cache_factor(T)
 pf=HC*np.cbrt(3*np.pi**2*d[:,6:10]);slope=np.sum(pf*np.hypot(MASS,pf),axis=1)/(3*HC**3)
 f=4*np.pi*r*r/np.sqrt(1-2*m/r)*np.exp(-nu)*slope*1e15*KB**2*CONV
 trap=float(np.trapezoid(f,r));mu=np.where(nb<=TON,mu_pe(nb),np.hypot(MN,pf[:,0]))*np.exp(nu)
 hp=json.loads((R/'ba13-onset-pilot/constants.json').read_text());ib=np.where(d[:,21]>=hp['eps_b']*model.geom)[0][-1];nub=np.log(e['mu_inf']/hp['mu_b'])
 # Fixed reconstructed piecewise-linear geometry gives a conservative support-shell bound.
 k=e['k'];a,b=r[k-1:k+1];vmax=4*np.pi*b*b/np.sqrt(1-2*max(m[k-1:k+1])/a)*np.exp(-6*min(nu[k-1:k+1]))*1e-36
 e.update(index=i,N_nu=list(map(float,nums)),N_nodal=list(map(float,nodal)),B_nodal=float(baryon),B_functional_relative=float(baryon/star['B_solved']-1),Lmu_nu=float(muon_L),
  production=coef,star=star,C_trapezoid=trap,C_decache=float(C[0]/T[0]/F[0]),cache_factors=F.tolist(),cache_reproduction=float(max(abs(C/(T*F)/trap-1))),
  mu_profile_constancy=float(max(abs(mu/e['mu_inf']-1))),
  neutron_cell=[float(a),float(b)],support_cell_L_upper=float(vmax*(b-a)),support_cell_relative=float(vmax*(b-a)/e['L_nu']),
  envelope_node=int(ib),Tb_node_factor=float(np.exp(-nu[ib])),Tb_physical_factor=float(np.exp(-nub)),
  Tb_offset=float(np.exp(nub-nu[ib])-1),Ts_offset=float(np.exp(.605*(nub-nu[ib]))-1),
  envelope_cell_logT_jump=float(nu[ib+1]-nu[ib]),
  hash_inputs={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [p/'profile.tsv',p/'model.txt',p/'certificate.json',raw/'coefficients.tsv',raw/'thermal_surface.tsv',raw/'target_stars.tsv']},
  cpu_seconds=time.process_time()-start,wall_seconds=time.time()-wall,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 return e
if __name__=='__main__':
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});i=int(sys.argv[1]);data=compute(i);data['job']=os.environ['SLURM_JOB_ID']
 text=json.dumps(data,indent=2,allow_nan=False)
 with (O/f'profile-v2-{i}.json').open('x') as f:f.write(text+'\n')
 print(i,data['cpu_seconds'],data['wall_seconds'])
