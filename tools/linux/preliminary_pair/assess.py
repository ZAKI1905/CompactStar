#!/usr/bin/env python3
"""Full numerical assessment against an analytically enclosed frozen-input model.

Every scientific PASS remains conditional until the model-normalization and
analytic-implementation assumptions are reviewed. Production criteria unchanged.
"""
import argparse,bisect,json,math,pathlib,sys
from decimal import Decimal as D,localcontext
from high_precision import Coupled,SOURCE,rows,lift,save,digest,ROOT
from thermal_enclosure import constants
from verify_chemical_reference import corner_uncertainty
from assess_source_retry import equilibrium_count_roundoff

KEYS=('x','e','mu')
POWERS=('P_dir_eq','P_dir_actual','LH','DeltaLnu','DeltaPbeta','Lnu_eq','Lnu_full','Lgamma','Lother','Pnet','Lout')

def power(model,y):
 T,C,direct,LH,nu,photon,rates=model.powers(y)
 eq=sum(model.luminosity)*T**8;delta=nu-eq
 # Exact polynomial identity in the owner-checked equilibrium subspace.
 # Avoid creating a Decimal summation residue by subtracting two equivalent
 # evaluations of the equilibrium luminosity. Production output is untouched.
 if y[1]==0 and y[2]==0:delta=D(0);nu=eq
 return dict(P_dir_eq=model.pdir,P_dir_actual=direct,LH=LH,DeltaLnu=delta,DeltaPbeta=LH-delta,Lnu_eq=eq,Lnu_full=nu,Lgamma=photon,Lother=D(0),Pnet=direct+LH-nu-photon,Lout=D(0),Tinf=T,Cstar=C,Echem=model.energy(y[1:]))

def enclosed_reference(out,role):
 model=Coupled(out,role);c=constants(model)
 hp=json.loads((out/f'high-precision/{role}-n8-p110.json').read_text())
 normalization=json.loads((out/'normalization.json').read_text())['roles'][role]
 result={};u=D(0);last=D(0)
 for row in hp['records']:
  t=D(row['t_s']);y=tuple(D(row[k]) for k in KEYS)
  if t>last:
   h=(t-last)/8;q=h/c['analytic_time_radius_s'];assert q<D('.01')
   local=c['R']*q**5/(1-q)+D('1e-80')
   for _ in range(8):u=(c['L']*h).exp()*u+local
  reach=u+c['M']*(t-last)
  assert c['real_box']['x_min']<y[0]-reach and y[0]+reach<c['real_box']['x_max']
  assert all(abs(y[j+1])+reach*c['scale'][j+1]<c['real_box']['eta_abs'][j] for j in range(2))
  # Variation of constants uses the anchor-model Lipschitz constant. The
  # primary-formula normalization difference is a uniformly bounded forcing.
  normalization_u=D(normalization['uniform_scaled_rhs_perturbation'])*t*(c['L']*t).exp()
  total=u+normalization_u
  result[t]=dict(y=y,u=tuple(total*s for s in c['scale']),scaled_u=total,analytic_RK4_u=u,normalization_u=normalization_u,powers=power(model,y))
  last=t
 return model,c,result

def assess(out,partial=False):
 result=dict(status='CONDITIONAL_NUMERICAL_ASSESSMENT; NOT FORMALLY QUALIFIED',source_qualified=False,control_qualified=False,Stage_B_launched=False,roles={},conditional_assumptions=['Frozen binary64 Cstar cells, geometry cache operands and chemical coefficients remain the accepted numerical inputs; background/parameter uncertainty is not requalified.','Photon/direct-power normalization is closed to its frozen-operand primary formula with an explicit propagated perturbation.','Analytic complex-domain/RK4/Hermite derivation and Decimal error padding are new unreviewed verification candidates.','No event occurred; actual crossing qualification is not established.'])
 paired=[]
 for role,base in [('source',SOURCE),('control',out/'control')]:
  model,c,hp=enclosed_reference(out,role)
  normalization=json.loads((out/'normalization.json').read_text())['roles'][role]
  normalization_power=D(normalization['direct_normalization_error'])+D(normalization['photon_normalization_error'])*(lift(2.42)*c['real_box']['x_max']).exp()
  refs=[rows(base/f'reference-{i}-observations.tsv') for i in range(5)]
  assert all(len(r)==1921 for r in refs)
  arithmetic=[json.loads((out/f'roundoff/{role}-reference-{i}.json').read_text()) for i in range(5)]
  observation_arithmetic=[json.loads((out/f'interpolation-roundoff/{role}-reference-{i}.json').read_text()) for i in range(5)]
  atimes=[[lift(r['time_s']) for r in a['records']] for a in arithmetic]
  diagnostics=rows(base/'tier-0-cap-0/diagnostics-1920.tsv')
  native_unc=rows(base/'tier-0-cap-0/uncertainties-1920.tsv')
  assert len(diagnostics)==len(native_unc)==1921
  states=[];power_ratios=[];state_ratios=[];temperature_unc=[];state_error_bound=[];power_error_bound=[];r18=[];u_reference=[]
  for i,d in enumerate(diagnostics):
   t=lift(d['t_s']);h=hp[t];ref=[lift(refs[4][i][k]) for k in KEYS]
   y=[lift(d[k]) for k in ('x_state','eta_e_MeV','eta_mu_MeV')]
   state_error_bound.append([abs(y[j]-h['y'][j])+h['u'][j] for j in range(3)])
   U=[]
   for j,k in enumerate(KEYS):
    vals=[lift(r[i][k]) for r in refs]
    differences=[abs(vals[a]-vals[b]) for a,b in ((0,1),(1,2),(3,4),(2,4))]
    accumulated=max(D(a['records'][min(bisect.bisect_left(times,t),len(times)-1)]['accumulated_arithmetic_enclosure'][j]) for a,times in zip(arithmetic,atimes)) if t else D(0)
    accumulated=max(accumulated,max(D(a['records'][i]['accumulated_arithmetic_enclosure'][j]) for a in observation_arithmetic))
    floor=64*lift(math.ulp(float(max(map(abs,vals)))))
    # Total independently enclosed reference error also covers interpolation.
    # Adding it to prefix arithmetic is conservative double counting, intentional.
    interpolation=max(abs(v-h['y'][j])+h['u'][j] for v in vals)
    u=2*max(differences)+max(floor,accumulated)+interpolation
    U.append(u)
    ultra=D('1e-16' if j==0 else '1e-22')+D('1e-11')*max(map(abs,vals))
    states.append(dict(time_s=t,component=k,reference_allocation_ratio=u/(D('.05')*ultra),reference_U=u,roundoff_component=accumulated,reference_total_error_enclosure=interpolation))
   u_reference.append(U)
   for tier in range(3):
    atol=[('1e-12','1e-18'),('1e-14','1e-20'),('1e-16','1e-22')][tier];rt=D(['1e-7','1e-9','1e-11'][tier])
    ratios=[]
    for j,key in enumerate(('Ux','Ue','Umu')):
     budget=D(atol[0 if j==0 else 1])+rt*max(abs(y[j]),abs(ref[j]))
     ratios.append((abs(y[j]-ref[j])+U[j]+lift(native_unc[i][key]))/budget)
    state_ratios.append(dict(time_s=t,tier=tier,ratios=ratios))
   G=max(D(1),sum(abs(lift(d[k])) for k in ('P_dir_actual_erg_s','LH_erg_s','DeltaLnu_erg_s','Lnu_eq_erg_s','Lgamma_erg_s','Lother_erg_s','L_out_fluid_inf_erg_s')))
   packet_unc={r['channel']:lift(r['U']) for r in rows(base/f'tier-0-cap-0/checkpoints-1920/point-{i}/diagnostics.tsv')}
   errors={}
   for k in POWERS:
    pu=6*c['complex_channel_bounds'][k]/c['R']*h['scaled_u']
    if k in ('P_dir_eq','P_dir_actual','Pnet'):pu+=D(normalization['direct_normalization_error'])
    if k in ('Lgamma','Pnet'):pu+=D(normalization['photon_normalization_error'])*(lift(2.42)*c['real_box']['x_max']).exp()
    native='L_out_fluid_inf_erg_s' if k=='Lout' else k+'_erg_s'
    errors[k]=abs(lift(d[native])-h['powers'][k])+pu
    rv=[lift(refs[z][i][k]) for z in range(5)]
    referr=max(abs(v-h['powers'][k])+pu for v in rv)+2*max(abs(rv[a]-rv[b]) for a,b in ((0,1),(1,2),(3,4),(2,4)))+64*lift(math.ulp(float(max(map(abs,rv)))))
    power_ratios.append(dict(time_s=t,channel=k,reference_allocation_ratio=referr/(D('5e-11')*G),main_refined_ultra_ratio=(errors[k]+referr+packet_unc.get(k,D(0)))/(D('1e-9')*G)))
   power_error_bound.append(errors)
   temperature_unc.append(abs(lift(d['Tinf_K'])-h['powers']['Tinf'])+h['powers']['Tinf']*(h['u'][0].exp()-1))
   r18.append(abs(lift(d['R18_residual_erg_s']))/G)
  ledger=rows(base/'tier-0-cap-0/R20-grids.tsv');assert len(ledger)==4
  last=ledger[-1];N=lift(last['N']);H=D(31557600);dt=H/1920
  count=equilibrium_count_roundoff(base/'tier-0-cap-0/diagnostics-1920.tsv',N) if role=='source' else {'supported_component_bound_erg':D(0),'scope':'Exact-zero count history.'}
  quad=2*max(abs(lift(a['residual'])-lift(b['residual'])) for a,b in zip(ledger,ledger[1:]))
  analytic_quad=H*dt**2/12*(2*c['complex_outgoing_power_bound']/c['analytic_time_radius_s']**2)
  # Thermal integration is a trapezoid in T, not a time quadrature.
  cell=c['cell'];Tmin=min(lift(d['Tinf_K']) for d in diagnostics);Tmax=max(lift(d['Tinf_K']) for d in diagnostics)
  Cmax=max(lift(d['Cstar_erg_K']) for d in diagnostics)*(1+D('1e-12'))
  thermal_quad=sum(abs(lift(b['Tinf_K'])-lift(a['Tinf_K']))**3 for a,b in zip(diagnostics,diagnostics[1:]))*abs(model.slopes[cell])/Tmin**2/12
  Umain=Cmax*(temperature_unc[0]+temperature_unc[-1])
  chemical_hp_unc,_=corner_uncertainty(model.energy,hp[H]['y'][1:],hp[H]['u'][1:])
  Umain+=(abs(lift(diagnostics[-1]['Echem_MeV'])-hp[H]['powers']['Echem'])+chemical_hp_unc)*model.erg
  Umain+=sum(dt/2*sum(a[k]+b[k] for k in ('Lnu_full','Lgamma','Lother','Lout')) for a,b in zip(power_error_bound,power_error_bound[1:]))
  # Exact Decimal reevaluation of the original published-grid operations.
  exact_du=sum((lift(a['Cstar_erg_K'])+lift(b['Cstar_erg_K']))*(lift(b['Tinf_K'])-lift(a['Tinf_K']))/2 for a,b in zip(diagnostics,diagnostics[1:]))
  exact_out=sum(dt/2*sum(lift(a[k])+lift(b[k]) for k in ('Lnu_full_erg_s','Lgamma_erg_s','Lother_erg_s','L_out_fluid_inf_erg_s')) for a,b in zip(diagnostics,diagnostics[1:]))
  quadrature_rounding=abs(exact_du-lift(last['delta_U']))+abs(exact_out-lift(last['outgoing_integral']))
  exact_eq=model.erg*lift(diagnostics[0]['mu_B_inf_MeV'])*(lift(diagnostics[-1]['B_count'])-lift(diagnostics[0]['B_count']))
  exact_chem=model.erg*(lift(diagnostics[-1]['Echem_MeV'])-lift(diagnostics[0]['Echem_MeV']))
  quadrature_rounding+=abs(lift(last['residual'])-(exact_eq+exact_chem+exact_du+exact_out))
  # Also account for Cstar evaluation at the published rounded temperatures.
  capacity_eval=D(0)
  for a,b in zip(diagnostics,diagnostics[1:]):
   def dc(d):
    x=(lift(d['Tinf_K'])/D('1e8')).ln();ideal=model.capacity[cell]+model.slopes[cell]*(x-model.knots[cell]);return abs(lift(d['Cstar_erg_K'])-ideal)
   capacity_eval+=abs(lift(b['Tinf_K'])-lift(a['Tinf_K']))*(dc(a)+dc(b))/2
  Umain+=capacity_eval
  Uquad=quad+analytic_quad+thermal_quad+quadrature_rounding
  subtotal=(abs(lift(last['residual']))+D(count['supported_component_bound_erg'])+lift(last['reconstruction_uncertainty'])+Umain+Uquad)/N
  restarts=[]
  for path in sorted(base.glob('tier-*/restart-comparison.tsv')):
   values=rows(path);tier=int(path.parent.name.split('-')[1]);atol=[('1e-12','1e-18'),('1e-14','1e-20'),('1e-16','1e-22')][tier];rt=D(['1e-7','1e-9','1e-11'][tier])
   tail=rows(path.parent/'restart-tail/main.accepted_states.tsv');left=lift(tail[0]['t_before_s'])
   prefix=next(r for r in rows(path.parent/'main.accepted_states.tsv') if lift(r['t_after_s'])==left)
   prefix_unc=[abs(lift(prefix[k])-hp[left]['y'][j])+hp[left]['u'][j] for j,k in enumerate(('x_state','eta_e_MeV','eta_mu_MeV'))]
   inherited=max(prefix_unc[j]/c['scale'][j] for j in range(3))*(c['L']*(H-left)).exp()
   ratios=[]
   for j,k in enumerate(KEYS):
    y=lift(values[0][k]);budget=D(atol[0 if j==0 else 1])+rt*max(abs(y),abs(hp[H]['y'][j]))
    ratios.append((abs(y-lift(refs[4][-1][k]))+u_reference[-1][j]+inherited*c['scale'][j])/budget)
   restarts.append(dict(path=str(path),prefix_uncertainty=prefix_unc,ratios=ratios))
  hashes={}
  for name in ('main.accepted_states.tsv','diagnostics-1920.tsv','uncertainties-1920.tsv','R20-grids.tsv'):
   paths=sorted(base.glob('tier-*/'+name))
   if not partial:assert len(paths)==12
   complete=[p for p in paths if name!='R20-grids.tsv' or len(rows(p))==4]
   hs={str(p):digest(p) for p in complete}
   if not partial:assert len(hs)==12 and len(set(hs.values()))==1,(role,name,'matrix not identical')
   hashes[name]=hs
  r=dict(reference_allocation_worst=max(states,key=lambda r:r['reference_allocation_ratio']),state_comparison_worst=max(state_ratios,key=lambda r:max(r['ratios'])),power_reference_worst=max(power_ratios,key=lambda r:r['reference_allocation_ratio']),power_main_worst=max(power_ratios,key=lambda r:r['main_refined_ultra_ratio']),maximum_R18_over_G=max(r18),maximum_temperature_error_enclosure_K=max(temperature_unc),endpoint_temperature_error_enclosure_K=temperature_unc[-1],global_HP_final_scaled_enclosure=hp[H]['scaled_u'],constants=c,restarts=restarts,matrix_hashes=hashes,energy=dict(original_R20_residual_erg=lift(last['residual']),N_erg=N,equilibrium_count_component=count,U_main_erg=Umain,U_reconstruction_erg=lift(last['reconstruction_uncertainty']),U_event_erg=D(0),U_quad_predeclared_erg=quad,U_quad_analytic_outgoing_erg=analytic_quad,U_quad_thermal_erg=thermal_quad,U_quad_arithmetic_erg=quadrature_rounding,U_quad_total_erg=Uquad,conditional_total_over_N=subtotal,threshold=D('2e-4'),actual_event_scope='No events; zero event term applies only to this trajectory.'),state_error_bound_max=[max(v[j] for v in state_error_bound) for j in range(3)],reference_U_max=[max(v[j] for v in u_reference) for j in range(3)],power_error_bound_max={k:max(v[k] for v in power_error_bound) for k in POWERS})
  result['roles'][role]=r;paired.append(temperature_unc)
  r['gates']=dict(reference_allocation=r['reference_allocation_worst']['reference_allocation_ratio']<=1,state_accuracy=max(r['state_comparison_worst']['ratios'])<=1,power_reference=r['power_reference_worst']['reference_allocation_ratio']<=1,power_accuracy=r['power_main_worst']['main_refined_ultra_ratio']<=1,R18=max(r18)<=D('1e-10'),R20=subtotal<=D('2e-4'),restart=len(restarts)==3 and all(max(x['ratios'])<=1 for x in restarts))
 result['paired_temperature_error_enclosure_K']=max(a+b for a,b in zip(*paired))
 result['endpoint_paired_temperature_error_enclosure_K']=paired[0][-1]+paired[1][-1]
 result['all_computed_conditional_gates_pass']=all(all(r['gates'].values()) for r in result['roles'].values())
 result['partial_matrix']=partial
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('out',type=pathlib.Path);p.add_argument('--partial',action='store_true');a=p.parse_args()
 with localcontext() as ctx:ctx.prec=110;result=assess(a.out,a.partial)
 save(a.out/('assessment-partial.json' if a.partial else 'assessment-final.json'),result)
 print(result['all_computed_conditional_gates_pass'],{k:v['gates'] for k,v in result['roles'].items()})
 if not a.partial:assert result['all_computed_conditional_gates_pass'],'a required conditional numerical gate failed'
