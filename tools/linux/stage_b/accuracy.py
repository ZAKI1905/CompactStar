#!/usr/bin/env python3
"""Direct full-error enclosures and explicitly partial original budget accounting.

An omitted arithmetic-only term is reported missing, never silently set to zero.
"""
import bisect,concurrent.futures,json,math,pathlib,sys
from decimal import Decimal as D,localcontext
from piecewise import MODEL,Coupled,domain,SCALE,END
from local_reference import advance
from high_precision import rows,lift
from run import save,digest
POWERS=('P_dir_eq','P_dir_actual','LH','DeltaLnu','DeltaPbeta','Lnu_eq','Lnu_full','Lout','Lgamma','Lother','Pnet')
NATIVE={'P_dir_eq':'P_dir_eq_erg_s','P_dir_actual':'P_dir_actual_erg_s','LH':'LH_erg_s','DeltaLnu':'DeltaLnu_erg_s','DeltaPbeta':'DeltaPbeta_erg_s','Lnu_eq':'Lnu_eq_erg_s','Lnu_full':'Lnu_full_erg_s','Lout':'L_out_fluid_inf_erg_s','Lgamma':'Lgamma_erg_s','Lother':'Lother_erg_s','Pnet':'Pnet_erg_s'}

def powers(m,y):
 from verify_chemical_reference import PI
 T,C,direct,LH,nu,photon,rates=m.powers(y)
 extra=sum(m.luminosity[j]*T**8*sum(D(v)*(y[j+1]/(m.kb*T*PI))**(2*k+2) for k,v in enumerate((22020,5670,420,9)))/11513 for j in range(2))
 return dict(P_dir_eq=m.pdir,P_dir_actual=direct,LH=LH,DeltaLnu=extra,DeltaPbeta=LH-extra,Lnu_eq=sum(m.luminosity)*T**8,Lnu_full=nu,Lout=D(0),Lgamma=photon,Lother=D(0),Pnet=direct+LH-nu-photon)

def assess_one(task):
 out,role,path=task
 with localcontext() as ctx:
  ctx.prec=110;m=Coupled(MODEL,role);box=domain(m);tier=int(path.name.split('-')[1]);rt=D(['1e-7','1e-9','1e-11'][tier]);at=[D(['1e-12','1e-14','1e-16'][tier]),D(['1e-18','1e-20','1e-22'][tier]),D(['1e-18','1e-20','1e-22'][tier])]
  hp=json.loads((out/f'hp-{role}-1000000.json').read_text())['records'];ht=[D(r['t_s']) for r in hp]
  forcing=D(json.loads((out/'preflight.json').read_text())['roles'][role]['normalization_rhs_bound'])
  norm=json.loads((MODEL/'normalization.json').read_text())['roles'][role]
  R=D('.1');cp=powers(m,(box['hi']+R,*[box['Y'][j]+R*SCALE[j+1] for j in range(2)]));Pmax=sum(abs(cp[k]) for k in ('P_dir_actual','LH','Lnu_full','Lgamma'))
  power_gradient=6*Pmax/R
  def reference(t):
   i=bisect.bisect_right(ht,t)-1;r=hp[i];y=tuple(D(r[k]) for k in ('x','e','mu'));u=D(r['U'])
   if ht[i]!=t:_,y,u=advance(m,box,ht[i],y,u,t)
   u+=forcing*t*(box['L']*t).exp()
   return y,u
  def component_u(u,j):return D(0) if role=='control' and j>0 else u*SCALE[j]
  data=rows(path/'diagnostics-1920.tsv');unc=rows(path/'uncertainties-1920.tsv')
  if len(data)!=1921 or len(unc)!=1921:return dict(role=role,configuration=path.name,status='INCOMPLETE')
  refs=[rows(out/'native'/role/f'reference-{k}-observations.tsv') for k in range(5)]
  assert all(len(r)==1921 for r in refs)
  maxima=dict(direct_state=D(0),known_terms_state=D(0),known_terms_reference_allocation=D(0),power=D(0),power_reference_allocation=D(0),accepted_direct_state=D(0));worst={};reference_temperature=[];reference_contraction_exceptions=0
  for index,(d,urow) in enumerate(zip(data,unc)):
   t=lift(d['t_s']);y,u=reference(t);pv=powers(m,y);reference_errors=[]
   for j,key in enumerate(('x','e','mu')):
    v=lift(d[('x_state','eta_e_MeV','eta_mu_MeV')[j]]);rv=[lift(r[index][key]) for r in refs]
    differences=[abs(rv[a]-rv[b]) for a,b in ((0,1),(1,2),(3,4),(2,4))]
    total=max(abs(z-y[j])+component_u(u,j) for z in rv)
    known=2*max(differences)+64*lift(math.ulp(float(max(map(abs,rv)))))+total
    if role=='control' and j>0:assert v==y[j]==0 and all(z==0 for z in rv);known=D(0)
    budget=at[j]+rt*max(abs(v),abs(rv[-1]));ultra=D('1e-16' if j==0 else '1e-22')+D('1e-11')*max(map(abs,rv))
    checks=dict(direct_state=(abs(v-y[j])+component_u(u,j))/(at[j]+rt*max(abs(v),abs(y[j]))),known_terms_state=(abs(v-rv[-1])+known+lift(urow[('Ux','Ue','Umu')[j]]))/budget,known_terms_reference_allocation=known/(D('.05')*ultra))
    for k,z in checks.items():
     if z>maxima[k]:maxima[k]=z;worst[k]=dict(time_s=t,component=key,ratio=z)
    if differences[1]>differences[0]/2:reference_contraction_exceptions+=1
   T=D('1e8')*y[0].exp();reference_temperature.append(abs(lift(refs[-1][index]['Tinf'])-T)+T*(u.exp()-1))
   G=max(D(1),sum(abs(lift(d[NATIVE[k]])) for k in ('P_dir_actual','LH','DeltaLnu','Lnu_eq','Lgamma','Lother','Lout')))
   output_scaled=max(lift(urow[k])/SCALE[j] for j,k in enumerate(('Ux','Ue','Umu')))
   for k in POWERS:
    normalization=D(0)
    if k in ('P_dir_eq','P_dir_actual','Pnet'):normalization+=D(norm['direct_normalization_error'])
    if k in ('Lgamma','Pnet'):normalization+=D(norm['photon_normalization_error'])*(D.from_float(2.42)*box['hi']).exp()
    pu=power_gradient*u+normalization;ou=power_gradient*output_scaled
    if k in ('Lout','Lother') or (role=='control' and k in ('P_dir_eq','P_dir_actual','LH','DeltaLnu','DeltaPbeta')):pu=ou=D(0)
    rv=[lift(r[index][k]) for r in refs];native=lift(d[NATIVE[k]])
    ru=max(abs(v-pv[k])+pu for v in rv)+2*max(abs(rv[a]-rv[b]) for a,b in ((0,1),(1,2),(3,4),(2,4)))+64*lift(math.ulp(float(max(map(abs,rv)))))
    checks=dict(power=(abs(native-rv[-1])+ru+ou)/(D('1e-7' if tier==0 else '1e-9')*G),power_reference_allocation=ru/(D('5e-11')*G))
    for key,z in checks.items():
     if z>maxima[key]:maxima[key]=z;worst[key]=dict(time_s=t,channel=k,ratio=z)
  history=rows(path/'main.accepted_states.tsv')
  for r in history:
   t=lift(r['t_after_s']);y,u=reference(t)
   for j,k in enumerate(('x_state','eta_e_MeV','eta_mu_MeV')):
    v=lift(r[k]);z=(abs(v-y[j])+component_u(u,j))/(at[j]+rt*max(abs(v),abs(y[j])))
    if z>maxima['accepted_direct_state']:maxima['accepted_direct_state']=z;worst['accepted_direct_state']=dict(time_s=t,component=k,ratio=z)
  restart=None
  if (path/'restart-comparison.tsv').exists():
   tail=rows(path/'restart-tail/main.accepted_states.tsv');a=lift(tail[0]['t_before_s']);prefix=next(r for r in history if lift(r['t_after_s'])==a)
   ya,ua=reference(a);ye,ue=reference(END)
   inherited=max((abs(lift(prefix[k])-ya[j])+component_u(ua,j))/SCALE[j] for j,k in enumerate(('x_state','eta_e_MeV','eta_mu_MeV')))*(box['L']*(END-a)).exp()
   final=rows(path/'restart-comparison.tsv')[0];ratios=[]
   for j,k in enumerate(('x','e','mu')):
    v=lift(final[k]);extra=D(0) if role=='control' and j>0 else inherited*SCALE[j]
    ratios.append((abs(v-ye[j])+component_u(ue,j)+extra)/(at[j]+rt*max(abs(v),abs(ye[j]))))
   restart=dict(time_s=a,prefix_inherited_scaled_bound=inherited,ratios=ratios)
  result=dict(role=role,configuration=path.name,status='PARTIAL_ORIGINAL_GATES_WITH_DIRECT_FULL_ERROR_ENCLOSURES',maxima=maxima,worst=worst,restart=restart,accepted_steps=len(history),maximum_reference_temperature_bound_K=max(reference_temperature),endpoint_reference_temperature_bound_K=reference_temperature[-1],reference_temperature_bounds_K=reference_temperature,noncontracting_state_comparisons=reference_contraction_exceptions,missing=['An actual-stage accumulated-arithmetic allowance is still required in the prescribed reference U_ref; known_terms_* excludes that unknown positive term.','Contraction exceptions need that independent arithmetic assessment.','Five native reference observations at every accepted main/restart endpoint are not supplied by this direct HP comparison.'],history_sha256=digest(path/'main.accepted_states.tsv'))
  target=out/'accuracy';target.mkdir(exist_ok=True);save(target/f'{role}-{path.name}.json',result)
  return {k:v for k,v in result.items() if k!='reference_temperature_bounds_K'}

if __name__=='__main__':
 out=pathlib.Path(sys.argv[1]);workers=int(sys.argv[2]);tasks=[]
 for role in ('source','control'):
  for p in sorted((out/'native'/role).glob('tier-*')):
   if (p/'diagnostics-1920.tsv').exists() and len(rows(p/'diagnostics-1920.tsv'))==1921:tasks.append((out,role,p))
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:results=list(pool.map(assess_one,tasks))
 save(out/'accuracy.json',dict(status='PRELIMINARY_ORIGINAL_REFERENCE_GATES_INCOMPLETE',workers=workers,configurations=results))
