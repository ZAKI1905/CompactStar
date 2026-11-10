#!/usr/bin/env python3
"""Accepted endpoints, convergence, matched norms and normalization closure.

This supplements (and does not weaken) the original production acceptance tests.
All analytic enclosures retain the assumptions stated in ENCLOSURE.md.
"""
import bisect,json,math,pathlib,sys
from decimal import Decimal as D,localcontext
from assess import enclosed_reference,POWERS
from high_precision import SOURCE,rows,lift,save,digest,Coupled
from thermal_enclosure import constants
from verify_chemical_reference import corner_uncertainty

def reference_u(values,t,arithmetic,component,observation_bounds=None):
 bounds=[]
 for a in arithmetic:
  ts=[lift(r['time_s']) for r in a['records']]
  bounds.append(D(a['records'][min(bisect.bisect_left(ts,t),len(ts)-1)]['accumulated_arithmetic_enclosure'][component]) if t else D(0))
 if observation_bounds is not None:bounds=[max(a,b) for a,b in zip(bounds,observation_bounds)]
 diffs=[abs(values[a]-values[b]) for a,b in ((0,1),(1,2),(3,4),(2,4))]
 floors=[64*lift(math.ulp(float(abs(v)))) for v in values]
 rounding=max(bounds+floors)
 # Refinement contraction is required only above the independently enclosed
 # arithmetic level. A sum encloses a difference between independent paths.
 under=[diffs[j]<=bounds[a]+bounds[b] for j,(a,b) in enumerate(((0,1),(1,2),(3,4),(2,4)))]
 contract=diffs[1]<=diffs[0]/2 or (under[0] and under[1])
 return 2*max(diffs)+rounding,dict(time_s=t,component=component,differences=diffs,roundoff=rounding,below_roundoff=under,rk8pd_contraction_or_roundoff=contract,independent_method_roundoff=under[2] and under[3])

def check(out,partial=False):
 result=dict(status='CONDITIONAL; NOT FORMALLY QUALIFIED',roles={})
 norm=json.loads((out/'normalization.json').read_text())
 paired={};H=D(31557600)
 for role,base in [('source',SOURCE),('control',out/'control')]:
  model,c,hp=enclosed_reference(out,role)
  refs=[rows(base/f'reference-{i}-observations.tsv') for i in range(5)]
  main=rows(base/'tier-0-cap-0/diagnostics-1920.tsv');native_u=rows(base/'tier-0-cap-0/uncertainties-1920.tsv')
  arithmetic=[json.loads((out/f'roundoff/{role}-reference-{i}.json').read_text()) for i in range(5)]
  observation_arithmetic=[json.loads((out/f'interpolation-roundoff/{role}-reference-{i}.json').read_text()) for i in range(5)]
  checks=[];accepted=[];accepted_power=[];state_data=[];energies={};energy_u={};reference_temperature_u=[]
  bounds={k:[] for k in POWERS};power_values={k:[] for k in POWERS}
  for i,d in enumerate(main):
   t=lift(d['t_s']);h=hp[t];U=[];state=[lift(d[k]) for k in ('x_state','eta_e_MeV','eta_mu_MeV')]
   reference_temperature_u.append(abs(lift(refs[4][i]['Tinf'])-h['powers']['Tinf'])+h['powers']['Tinf']*(h['u'][0].exp()-1))
   ref=[lift(refs[4][i][k]) for k in ('x','e','mu')]
   for j,k in enumerate(('x','e','mu')):
    vals=[lift(r[i][k]) for r in refs];u,a=reference_u(vals,t,arithmetic,j,[D(r['records'][i]['accumulated_arithmetic_enclosure'][j]) for r in observation_arithmetic])
    u+=max(abs(v-h['y'][j])+h['u'][j] for v in vals)
    checks.append(a);U.append(u+lift(native_u[i][('Ux','Ue','Umu')[j]]))
   state_data.append((state,ref,U))
   for k in POWERS:
    val=lift(refs[4][i][k]);pu=6*c['complex_channel_bounds'][k]/c['R']*h['scaled_u']
    if k in ('P_dir_eq','P_dir_actual','Pnet'):pu+=D(norm['roles'][role]['direct_normalization_error'])
    if k in ('Lgamma','Pnet'):pu+=D(norm['roles'][role]['photon_normalization_error'])*(lift(2.42)*c['real_box']['x_max']).exp()
    if role=='control' and k in ('P_dir_eq','P_dir_actual','LH','DeltaLnu','DeltaPbeta','Lother','Lout'):
     assert val==0 and h['powers'][k]==0
     pu=D(0)
    bounds[k].append(abs(val-h['powers'][k])+pu);power_values[k].append(val)
  for path in sorted((out/f'accepted-observations/{role}').glob('point-*.tsv')):
   rr=rows(path);assert len(rr)==6
   t=lift(rr[0]['t']);h=hp[t]
   for j,k in enumerate(('x','e','mu')):
    vals=[lift(r[k]) for r in rr[1:]];u,a=reference_u(vals,t,arithmetic,j)
    u+=max(abs(v-h['y'][j])+h['u'][j] for v in vals);checks.append(a)
    y=lift(rr[0][k]);budget=D('1e-16' if j==0 else '1e-22')+D('1e-11')*max(abs(y),abs(vals[-1]))
    alloc=D('1e-16' if j==0 else '1e-22')+D('1e-11')*max(map(abs,vals))
    accepted.append(dict(time_s=t,component=k,ultra_ratio=(abs(y-vals[-1])+u)/budget,reference_allocation_ratio=u/(D('.05')*alloc)))
   G=max(D(1),sum(abs(lift(rr[0][k])) for k in ('P_dir_actual','LH','DeltaLnu','Lnu_eq','Lgamma','Lother','Lout')))
   for k in POWERS:
    pu=6*c['complex_channel_bounds'][k]/c['R']*h['scaled_u']
    if k in ('P_dir_eq','P_dir_actual','Pnet'):pu+=D(norm['roles'][role]['direct_normalization_error'])
    if k in ('Lgamma','Pnet'):pu+=D(norm['roles'][role]['photon_normalization_error'])*(lift(2.42)*c['real_box']['x_max']).exp()
    vals=[lift(r[k]) for r in rr[1:]]
    u=max(abs(v-h['powers'][k])+pu for v in vals)+2*max(abs(vals[a]-vals[b]) for a,b in ((0,1),(1,2),(3,4),(2,4)))+64*lift(math.ulp(float(max(map(abs,vals)))))
    accepted_power.append(dict(time_s=t,channel=k,ultra_ratio=(abs(lift(rr[0][k])-vals[-1])+u)/(D('1e-9')*G),reference_allocation_ratio=u/(D('5e-11')*G)))
  assert len(accepted)==3*len(rows(base/'tier-0-cap-0/main.accepted_states.tsv'))
  dt=H/1920
  for k in POWERS:
   energies[k]=sum(dt*(a+b)/2 for a,b in zip(power_values[k],power_values[k][1:]))
   propagated=sum(dt*(a+b)/2 for a,b in zip(bounds[k],bounds[k][1:]))
   # Cauchy second derivative bound for a channel composed with the exact flow.
   quad=H*dt**2/12*(2*c['complex_channel_bounds'][k]/c['analytic_time_radius_s']**2)
   if role=='control' and k in ('P_dir_eq','P_dir_actual','LH','DeltaLnu','DeltaPbeta','Lother','Lout'):quad=D(0)
   energy_u[k]=dict(propagated_state_and_evaluation_erg=propagated,analytic_trapezoid_erg=quad,total_erg=propagated+quad)
  ledger=rows(base/'tier-0-cap-0/R20-grids.tsv');diffs=[abs(lift(a['residual'])-lift(b['residual'])) for a,b in zip(ledger,ledger[1:])]
  assessment=json.loads((out/('assessment-partial.json' if partial else 'assessment-final.json')).read_text())['roles'][role]
  rounding=D(assessment['energy']['U_quad_arithmetic_erg'])+D(assessment['energy']['equilibrium_count_component']['supported_component_bound_erg'])
  quad_ok=diffs[-1]<=diffs[-2]/2 or max(diffs)<=2*rounding
  floor=max(lift(r['finite_T_omitted_floor_erg_s']) for r in main)
  temp_floor=floor*H/c['real_box']['C_min']
  gates=dict(accepted_state=max(r['ultra_ratio'] for r in accepted)<=1,accepted_reference=max(r['reference_allocation_ratio'] for r in accepted)<=1,accepted_power=max(r['ultra_ratio'] for r in accepted_power)<=1,accepted_power_reference=max(r['reference_allocation_ratio'] for r in accepted_power)<=1,reference_contraction=all(r['rk8pd_contraction_or_roundoff'] and r['independent_method_roundoff'] for r in checks),R20_contraction=quad_ok,reconstruction_event_allocation=D(assessment['energy']['U_reconstruction_erg'])/D(assessment['energy']['N_erg'])<=D('5e-6'))
  result['roles'][role]=dict(gates=gates,accepted_state_worst=max(accepted,key=lambda r:r['ultra_ratio']),accepted_reference_worst=max(accepted,key=lambda r:r['reference_allocation_ratio']),accepted_power_worst=max(accepted_power,key=lambda r:r['ultra_ratio']),accepted_power_reference_worst=max(accepted_power,key=lambda r:r['reference_allocation_ratio']),reference_convergence_checks=len(checks),reference_convergence_exceptions=[r for r in checks if not(r['rk8pd_contraction_or_roundoff'] and r['independent_method_roundoff'])],R20_nested_differences=diffs,R20_independent_rounding_erg=rounding,energies_erg=energies,energy_uncertainties=energy_u,maximum_finite_T_omitted_floor_erg_s=floor,finite_T_floor_one_year_temperature_scale_K=temp_floor,finite_T_scope='Conservative scale of a declared omitted physical-model term, not numerical roundoff and not a requalified stellar uncertainty bound.',maximum_frozen_utilization=max(lift(r['max_frozen_utilization']) for r in main),maximum_depletion=max(abs(lift(r['DeltaB_over_B0'])) for r in main),accepted_steps=len(accepted)//3)
  if role=='control':
   local_control=[]
   for i in range(1,1920):
    directory=base/f'tier-0-cap-0/checkpoints-1920/point-{i}'
    cp=json.loads((directory/'checkpoint.json').read_text());t=lift(cp['time_s']);left=lift(cp['main_left_time_s']);h=hp[t];hl=hp[left]
    # Exact-zero chemistry is the independently owner-checked invariant of
    # this zero-source context. Only thermal prefix error can propagate here.
    assert cp['O2_state'][1:]==[0,0] and cp['main_left_state'][1:]==[0,0]
    prefix=abs(lift(cp['main_left_state'][0])-hl['y'][0])+hl['u'][0]
    local_x=abs(lift(cp['O2_state'][0])-h['y'][0])+h['u'][0]+prefix*(c['L']*(t-left)).exp()
    sr={r['channel']:r for r in rows(directory/'state.tsv')};dr={r['channel']:r for r in rows(directory/'diagnostics.tsv')}
    assert all(lift(dr[k]['O2'])==0 for k in ('Echem','LH','DeltaLnu','DeltaPbeta','P_dir_eq','P_dir_actual'))
    local_control.append(dict(point=i,thermal_error=local_x,ratio=local_x/lift(sr['x']['propagated_state_radius'])))
   gates['local_reconstruction']=all(r['ratio']<=1 for r in local_control)
   result['roles'][role]['local_reconstruction_worst']=max(local_control,key=lambda r:r['ratio'])
  paired[role]=dict(states=state_data,powers=power_values,power_bounds=bounds,temperature_u=reference_temperature_u,temperature=[lift(r['Tinf']) for r in refs[4]])
 pair_ratios=[];state_signs=[0,0,0];power_signs={k:0 for k in POWERS}
 for i in range(1921):
  sy,sr,su=paired['source']['states'][i];cy,cr,cu=paired['control']['states'][i]
  for j in range(3):
   budget=2*D('1e-16' if j==0 else '1e-22')+D('1e-11')*(max(abs(sy[j]),abs(sr[j]))+max(abs(cy[j]),abs(cr[j])))
   delta=sr[j]-cr[j];u=su[j]+cu[j]
   pair_ratios.append((abs((sy[j]-cy[j])-delta)+u)/budget)
   if i and u<=D('.1')*abs(delta):state_signs[j]+=1
  for k in POWERS:
   delta=paired['source']['powers'][k][i]-paired['control']['powers'][k][i]
   u=paired['source']['power_bounds'][k][i]+paired['control']['power_bounds'][k][i]
   if i and delta and u<=D('.1')*abs(delta):power_signs[k]+=1
 result['paired']=dict(worst_ultra_state_ratio=max(pair_ratios),state_sign_resolved_nonzero_samples=state_signs,power_sign_resolved_nonzero_samples=power_signs,noninitial_samples=1920,uncertainties_added_without_cancellation=True)
 temperature_u=[a+b for a,b in zip(paired['source']['temperature_u'],paired['control']['temperature_u'])]
 result['paired'].update(plotted_reference_endpoint_temperature_uncertainty_K=temperature_u[-1],plotted_reference_maximum_temperature_uncertainty_K=max(temperature_u),temperature_sign_resolved_nonzero_samples=sum(temperature_u[i]<=D('.1')*abs(paired['source']['temperature'][i]-paired['control']['temperature'][i]) for i in range(1,1921)))
 result['paired_energies']={}
 for k in POWERS:
  delta=result['roles']['source']['energies_erg'][k]-result['roles']['control']['energies_erg'][k]
  u=sum(result['roles'][role]['energy_uncertainties'][k]['total_erg'] for role in ('source','control'))
  result['paired_energies'][k]=dict(source_minus_control_erg=delta,conditional_uncertainty_erg=u,sign_status='EXACT_ZERO' if delta==u==0 else 'CONDITIONALLY_RESOLVED' if u<=D('.1')*abs(delta) else 'SIGN_UNRESOLVED')
 # Close the local continuous enclosure and coupled endpoint tests to the
 # primary normalization formula. Retained left states remain exact inputs.
 local=[];model=Coupled(out,'source');c=constants(model);rhs=D(norm['roles']['source']['uniform_scaled_rhs_perturbation'])
 for path in sorted((out/'thermal-enclosure').glob('point-*.json')):
  a=json.loads(path.read_text());i=a['point'];directory=next(p for p in (out/'sample'/f'point-{i}',out/'full'/f'point-{i}') if p.exists())
  cp=json.loads((directory/'packet/checkpoint.json').read_text());span=lift(cp['time_s'])-lift(cp['main_left_time_s']);u=D(a['reference_scaled_error_bound'])+rhs*span*(c['L']*span).exp()
  y=tuple(D(v) for v in a['reference_final']);sr={r['channel']:r for r in rows(directory/'packet/state.tsv')};dr={r['channel']:r for r in rows(directory/'packet/diagnostics.tsv')}
  ratios=[(abs(lift(cp['O2_state'][j])-y[j])+u*c['scale'][j])/lift(sr[k]['propagated_state_radius']) for j,k in enumerate(('x','eta_e','eta_mu'))]
  ue,_=corner_uncertainty(model.energy,y[1:],tuple(u*c['scale'][j] for j in (1,2)))
  er=(abs(lift(dr['Echem']['O2'])-model.energy(y[1:]))+ue)/lift(dr['Echem']['U'])
  local.append(dict(point=i,thermal=D(a['uniform_thermal_enclosure'])+rhs*span*(c['L']*span).exp(),state_ratio=max(ratios),energy_ratio=er))
 assert len(local)==1919
 result['local_primary_formula_closure']=dict(points=len(local),worst_thermal=max(local,key=lambda r:r['thermal']),worst_state=max(local,key=lambda r:r['state_ratio']),worst_energy=max(local,key=lambda r:r['energy_ratio']),conditional_pass=all(r['thermal']<=D('1e-12') and max(r['state_ratio'],r['energy_ratio'])<=1 for r in local))
 result['all_additional_conditional_gates_pass']=all(all(r['gates'].values()) for r in result['roles'].values()) and max(pair_ratios)<=1 and result['local_primary_formula_closure']['conditional_pass']
 return result

if __name__=='__main__':
 out=pathlib.Path(sys.argv[1])
 partial='--partial' in sys.argv
 with localcontext() as ctx:ctx.prec=110;result=check(out,partial)
 result['partial_matrix']=partial
 save(out/('supplement-partial.json' if partial else 'supplement.json'),result)
 print(result['all_additional_conditional_gates_pass'],{k:v['gates'] for k,v in result['roles'].items()})
 if not partial:assert result['all_additional_conditional_gates_pass'],'an additional conditional numerical gate failed'
