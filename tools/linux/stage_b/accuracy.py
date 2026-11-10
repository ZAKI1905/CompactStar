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
  arithmetic=[]
  for k in range(5):
   f=out/'interpolation-roundoff'/f'{role}-reference-{k}.json'
   if f.exists():arithmetic.append({D(r['time_s']):[D(v) for v in r['accumulated_arithmetic_enclosure']] for r in json.loads(f.read_text())['records']})
  complete_arithmetic=len(arithmetic)==5
  accepted_refs={}
  for f in (out/f'{role}-supplement').glob('accepted-*.tsv'):
   rr=rows(f);assert len(rr)==5;accepted_refs[lift(rr[0]['t'])]=rr
  def au(k,t,j):
   if not complete_arithmetic:return D(0)
   if role=='control' and j>0:return D(0)
   return arithmetic[k][t][j]
  maxima=dict(direct_state=D(0),known_terms_state=D(0),known_terms_reference_allocation=D(0),power=D(0),power_reference_allocation=D(0),accepted_direct_state=D(0));worst={};reference_temperature=[];reference_contraction_exceptions=0;unresolved_contractions=[];unresolved_power_contractions=[]
  for index,(d,urow) in enumerate(zip(data,unc)):
   t=lift(d['t_s']);y,u=reference(t);pv=powers(m,y);reference_errors=[]
   for j,key in enumerate(('x','e','mu')):
    v=lift(d[('x_state','eta_e_MeV','eta_mu_MeV')[j]]);rv=[lift(r[index][key]) for r in refs]
    differences=[abs(rv[a]-rv[b]) for a,b in ((0,1),(1,2),(3,4),(2,4))]
    total=max(abs(z-y[j])+component_u(u,j) for z in rv)
    known=2*max(differences)+64*lift(math.ulp(float(max(map(abs,rv)))))+max(au(k,t,j) for k in range(5))+total
    if role=='control' and j>0:assert v==y[j]==0 and all(z==0 for z in rv);known=D(0)
    budget=at[j]+rt*max(abs(v),abs(rv[-1]));ultra=D('1e-16' if j==0 else '1e-22')+D('1e-11')*max(map(abs,rv))
    checks=dict(direct_state=(abs(v-y[j])+component_u(u,j))/(at[j]+rt*max(abs(v),abs(y[j]))),known_terms_state=(abs(v-rv[-1])+known+lift(urow[('Ux','Ue','Umu')[j]]))/budget,known_terms_reference_allocation=known/(D('.05')*ultra))
    for k,z in checks.items():
     if z>maxima[k]:maxima[k]=z;worst[k]=dict(time_s=t,component=key,ratio=z)
    if differences[1]>differences[0]/2:
     reference_contraction_exceptions+=1
     if not complete_arithmetic or not (differences[0]<=au(0,t,j)+au(1,t,j) and differences[1]<=au(1,t,j)+au(2,t,j)):unresolved_contractions.append(dict(time_s=t,component=key,differences=differences[:2]))
   T=D('1e8')*y[0].exp();reference_temperature.append(abs(lift(refs[-1][index]['Tinf'])-T)+T*(u.exp()-1))
   G=max(D(1),sum(abs(lift(d[NATIVE[k]])) for k in ('P_dir_actual','LH','DeltaLnu','Lnu_eq','Lgamma','Lother','Lout')))
   output_scaled=max(lift(urow[k])/SCALE[j] for j,k in enumerate(('Ux','Ue','Umu')))
   reference_power_values=[powers(m,tuple(lift(r[index][q]) for q in ('x','e','mu'))) for r in refs]
   for k in POWERS:
    normalization=D(0)
    if k in ('P_dir_eq','P_dir_actual','Pnet'):normalization+=D(norm['direct_normalization_error'])
    if k in ('Lgamma','Pnet'):normalization+=D(norm['photon_normalization_error'])*(D.from_float(2.42)*box['hi']).exp()
    pu=power_gradient*u+normalization;ou=power_gradient*output_scaled
    if k in ('Lout','Lother') or (role=='control' and k in ('P_dir_eq','P_dir_actual','LH','DeltaLnu','DeltaPbeta')):pu=ou=D(0)
    rv=[lift(r[index][k]) for r in refs];native=lift(d[NATIVE[k]])
    individual_arithmetic=[power_gradient*max(au(i,t,j)/SCALE[j] for j in range(3))+abs(lift(refs[i][index][k])-reference_power_values[i][k])+normalization+D('1e-60')*max(D(1),abs(reference_power_values[i][k])) for i in range(5)]
    arithmetic_power=max(individual_arithmetic)
    if k in ('Lout','Lother') or (role=='control' and k in ('P_dir_eq','P_dir_actual','LH','DeltaLnu','DeltaPbeta')):arithmetic_power=D(0)
    ru=arithmetic_power+max(abs(v-pv[k])+pu for v in rv)+2*max(abs(rv[a]-rv[b]) for a,b in ((0,1),(1,2),(3,4),(2,4)))+64*lift(math.ulp(float(max(map(abs,rv)))))
    d0=abs(rv[0]-rv[1]);d1=abs(rv[1]-rv[2])
    if d1>d0/2 and (not complete_arithmetic or not (d0<=individual_arithmetic[0]+individual_arithmetic[1] and d1<=individual_arithmetic[1]+individual_arithmetic[2])):unresolved_power_contractions.append(dict(time_s=t,channel=k))
    checks=dict(power=(abs(native-rv[-1])+ru+ou)/(D('1e-7' if tier==0 else '1e-9')*G),power_reference_allocation=ru/(D('5e-11')*G))
    for key,z in checks.items():
     if z>maxima[key]:maxima[key]=z;worst[key]=dict(time_s=t,channel=k,ratio=z)
  history=rows(path/'main.accepted_states.tsv');accepted=[];restart=None
  def accepted_check(r,prefix_extra=None,label='main'):
   t=lift(r['t_after_s']);y,u=reference(t);ratios=[];allocations=[];direct=[];unresolved=[]
   for j,k in enumerate(('x_state','eta_e_MeV','eta_mu_MeV')):
    v=lift(r[k]);extra=(prefix_extra or (D(0),)*3)[j]
    direct.append((abs(v-y[j])+component_u(u,j)+extra)/(at[j]+rt*max(abs(v),abs(y[j]))))
    if t in accepted_refs:
     rv=[lift(q[('x','e','mu')[j]]) for q in accepted_refs[t]]
     differences=[abs(rv[a]-rv[b]) for a,b in ((0,1),(1,2),(3,4),(2,4))]
     ru=max(abs(q-y[j])+component_u(u,j) for q in rv)+2*max(differences)+max(au(i,t,j) for i in range(5))+64*lift(math.ulp(float(max(map(abs,rv)))))
     if role=='control' and j>0:assert all(q==0 for q in rv) and v==0;ru=D(0)
     ratios.append((abs(v-rv[-1])+ru+extra)/(at[j]+rt*max(abs(v),abs(rv[-1]))))
     allocations.append(ru/(D('.05')*(D('1e-16' if j==0 else '1e-22')+D('1e-11')*max(map(abs,rv)))))
     if differences[1]>differences[0]/2 and (not complete_arithmetic or not (differences[0]<=au(0,t,j)+au(1,t,j) and differences[1]<=au(1,t,j)+au(2,t,j))):unresolved.append(j)
   return dict(time_s=t,kind=label,prefix_extra=prefix_extra or [D(0)]*3,direct_ratios=direct,prescribed_state_ratios=ratios,reference_allocations=allocations,unresolved_contractions=unresolved)
  for r in history:accepted.append(accepted_check(r))
  maxima['accepted_direct_state']=max(max(r['direct_ratios']) for r in accepted)
  if (path/'restart-comparison.tsv').exists():
   tail=rows(path/'restart-tail/main.accepted_states.tsv');a=lift(tail[0]['t_before_s']);prefix=next(r for r in history if lift(r['t_after_s'])==a)
   ya,ua=reference(a)
   initial=max((abs(lift(prefix[k])-ya[j])+component_u(ua,j))/SCALE[j] for j,k in enumerate(('x_state','eta_e_MeV','eta_mu_MeV')))
   records=[]
   for r in tail:
    inherited=initial*(box['L']*(lift(r['t_after_s'])-a)).exp()
    records.append(accepted_check(r,[D(0) if role=='control' and j>0 else inherited*SCALE[j] for j in range(3)],'restart'))
   restart=dict(time_s=a,prefix_initial_scaled_bound=initial,records=records)
  accepted_power=[];missing_native_power=False
  for label,records in [('main',accepted),('restart',restart['records'] if restart else [])]:
   if not records:continue
   native_path=out/f'{role}-native-points'/(path.name+('-restart' if label=='restart' else '')+'.tsv')
   if not native_path.exists():missing_native_power=True;continue
   native_rows=rows(native_path);assert len(native_rows)==len(records)
   for q,r in zip(native_rows,records):
    t=lift(q['t']);assert t==r['time_s']
    if t not in accepted_refs:continue
    y,u=reference(t);pv=powers(m,y);ref=accepted_refs[t]
    reference_power_values=[powers(m,tuple(lift(q[k]) for k in ('x','e','mu'))) for q in ref]
    G=max(D(1),sum(abs(lift(q[k])) for k in ('P_dir_actual','LH','DeltaLnu','Lnu_eq','Lgamma','Lother','Lout')))
    inherited=power_gradient*max(r['prefix_extra'][j]/SCALE[j] for j in range(3));largest=D(0);allocation=D(0)
    for k in POWERS:
     normalization=D(0)
     if k in ('P_dir_eq','P_dir_actual','Pnet'):normalization+=D(norm['direct_normalization_error'])
     if k in ('Lgamma','Pnet'):normalization+=D(norm['photon_normalization_error'])*(D.from_float(2.42)*box['hi']).exp()
     individual_arithmetic=[power_gradient*max(au(i,t,j)/SCALE[j] for j in range(3))+abs(lift(ref[i][k])-reference_power_values[i][k])+normalization+D('1e-60')*max(D(1),abs(reference_power_values[i][k])) for i in range(5)]
     pu=power_gradient*u+normalization;ar=max(individual_arithmetic);extra=inherited
     if k in ('Lout','Lother') or (role=='control' and k in ('P_dir_eq','P_dir_actual','LH','DeltaLnu','DeltaPbeta')):pu=ar=extra=D(0)
     rv=[lift(v[k]) for v in ref];native=lift(q[k])
     ru=max(abs(v-pv[k])+pu for v in rv)+ar+2*max(abs(rv[a]-rv[b]) for a,b in ((0,1),(1,2),(3,4),(2,4)))+64*lift(math.ulp(float(max(map(abs,rv)))))
     d0=abs(rv[0]-rv[1]);d1=abs(rv[1]-rv[2])
     if d1>d0/2 and (not complete_arithmetic or not (d0<=individual_arithmetic[0]+individual_arithmetic[1] and d1<=individual_arithmetic[1]+individual_arithmetic[2])):unresolved_power_contractions.append(dict(time_s=t,channel=k,kind=label))
     largest=max(largest,(abs(native-rv[-1])+ru+extra)/(D('1e-7' if tier==0 else '1e-9')*G));allocation=max(allocation,ru/(D('5e-11')*G))
    accepted_power.append(dict(time_s=t,kind=label,budget_ratio=largest,reference_allocation_ratio=allocation))
  missing=[]
  if unresolved_power_contractions:missing.append('Power reference contraction exceptions not bounded by actual arithmetic.')
  if missing_native_power:missing.append('Native power evaluations absent at accepted main/restart endpoints.')
  if not complete_arithmetic:missing.append('Actual observation accumulated arithmetic is missing; known_terms_* excludes that unknown positive term.')
  if any(len(r['prescribed_state_ratios'])!=3 for r in accepted+(restart['records'] if restart else [])):missing.append('Native reference observations missing at accepted main/restart endpoints.')
  if unresolved_contractions or any(r['unresolved_contractions'] for r in accepted+(restart['records'] if restart else [])):missing.append('State reference contraction exceptions not bounded by actual arithmetic.')
  state_records=accepted+(restart['records'] if restart else [])
  gates=dict(power_refinement_contraction=not unresolved_power_contractions,accepted_powers=not missing_native_power and all(r['budget_ratio']<=1 for r in accepted_power),accepted_reference_power_allocation=not missing_native_power and all(r['reference_allocation_ratio']<=1 for r in accepted_power),scientific_states=maxima['known_terms_state']<=1,reference_state_allocation=maxima['known_terms_reference_allocation']<=1,powers=maxima['power']<=1,reference_power_allocation=maxima['power_reference_allocation']<=1,accepted_restart_states=all(max(r['prescribed_state_ratios'],default=D('Infinity'))<=1 for r in state_records),accepted_reference_allocation=all(max(r['reference_allocations'],default=D('Infinity'))<=1 for r in state_records),state_refinement_contraction=not unresolved_contractions and all(not r['unresolved_contractions'] for r in state_records))
  result=dict(role=role,configuration=path.name,status='CONDITIONAL_STATE_POWER_BUDGET_CANDIDATE' if not missing and all(gates.values()) else 'PRELIMINARY_INCOMPLETE_OR_FAILED_BUDGET',gates=gates,complete_arithmetic=complete_arithmetic,accepted=accepted,accepted_power=accepted_power,unresolved_contractions=unresolved_contractions,unresolved_power_contractions=unresolved_power_contractions,maxima=maxima,worst=worst,restart=restart,accepted_steps=len(history),maximum_reference_temperature_bound_K=max(reference_temperature),endpoint_reference_temperature_bound_K=reference_temperature[-1],reference_temperature_bounds_K=reference_temperature,noncontracting_state_comparisons=reference_contraction_exceptions,missing=missing,history_sha256=digest(path/'main.accepted_states.tsv'))
  target=out/'accuracy';target.mkdir(exist_ok=True);save(target/f'{role}-{path.name}.json',result)
  return {k:v for k,v in result.items() if k not in ('reference_temperature_bounds_K','accepted','accepted_power')}

if __name__=='__main__':
 out=pathlib.Path(sys.argv[1]);workers=int(sys.argv[2]);tasks=[]
 for role in ('source','control'):
  for p in sorted((out/'native'/role).glob('tier-*')):
   if (p/'diagnostics-1920.tsv').exists() and len(rows(p/'diagnostics-1920.tsv'))==1921:tasks.append((out,role,p))
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:results=list(pool.map(assess_one,tasks))
 save(out/'accuracy.json',dict(status='CONDITIONAL_BUDGET_CANDIDATE_REQUIRING_ALL_OTHER_GATES',workers=workers,configurations=results))
