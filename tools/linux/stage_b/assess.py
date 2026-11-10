#!/usr/bin/env python3
"""Partial, explicit Stage B assessment. Never treats missing gates as zero/PASS."""
import concurrent.futures,csv,json,math,pathlib,sys
from decimal import Decimal as D,localcontext
from piecewise import ROOT,MODEL,Coupled,END,YEAR,SCALE,domain,rows,lift
from run import save,digest

KEYS=('x','e','mu')
def physics(task):
 role,record=task
 with localcontext() as ctx:
  ctx.prec=110;m=Coupled(MODEL,role);y=tuple(D(record[k]) for k in KEYS)
  T,C,direct,LH,nu,photon,rates=m.powers(y)
  # Evaluate excess neutrinos directly, not by subtracting nearly equal powers.
  from verify_chemical_reference import PI
  extra=sum(m.luminosity[j]*T**8*sum(D(v)*(y[j+1]/(m.kb*T*PI))**(2*k+2) for k,v in enumerate((22020,5670,420,9)))/11513 for j in range(2))
  return dict(t_s=D(record['t_s']),T=T,C=C,x=y[0],e=y[1],mu=y[2],U=D(record['U']),Pdir=direct,LH=LH,DeltaLnu=extra,Lnu=nu,Lgamma=photon,Echem=m.energy(y[1:])*m.erg,Pnet=direct+LH-nu-photon)

def assess(out,workers):
 pre=json.loads((out/'preflight.json').read_text());result=dict(status='PRELIMINARY_NOT_STAGE_B_QUALIFIED',roles={},missing_obligations=['Actual-stage accumulated arithmetic including event trials and reference observation interpolation.','Full per-configuration main/accepted endpoint/restart budgets and all eleven power allocations.','Every local reconstruction enclosure, crossing-neighbor and event uncertainty assessment.','Complete original R20 quadrature/main/reconstruction/event uncertainties.'],resources=dict(assessment_workers=workers))
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
  for role in ('source','control'):
   hp=json.loads((out/f'hp-{role}-1000000.json').read_text())
   values=list(pool.map(physics,[(role,r) for r in hp['records']],chunksize=16))
   with localcontext() as ctx:
    ctx.prec=110;m=Coupled(MODEL,role);box=domain(m);forcing=D(pre['roles'][role]['normalization_rhs_bound'])
    for v in values:
     v['U_primary']=v['U']+forcing*v['t_s']*(box['L']*v['t_s']).exp()
     v['U_T']=v['T']*(v['U_primary'].exp()-1)
    energies={k:sum((b['t_s']-a['t_s'])*(a[k]+b[k])/2 for a,b in zip(values,values[1:])) for k in ('Pdir','LH','DeltaLnu','Lnu','Lgamma','Pnet')}
    checks=[];reference_errors=[]
    for i in range(5):
     path=out/'native'/role/f'reference-{i}-observations.tsv'
     if not path.exists():continue
     data=rows(path)
     if len(data)!=1921:continue
     maximum=[D(0)]*3
     for r,v in zip(data,values):
      assert lift(r['t'])==v['t_s']
      for j,k in enumerate(KEYS):maximum[j]=max(maximum[j],abs(lift(r[k])-v[k])+v['U_primary']*SCALE[j])
     reference_errors.append(dict(index=i,maximum_total_error=maximum))
    native=out/'native'/role
    configurations=[]
    for path in sorted(native.glob('tier-*-cap-*')):
     histories=path/'main.accepted_states.tsv';eventfile=path/'events.tsv';d=dict(name=path.name)
     if histories.exists():
      hist=rows(histories);d.update(accepted_steps=len(hist),history_sha256=digest(histories),final=hist[-1] if hist else None)
     if eventfile.exists():d['events']=[r for r in rows(eventfile) if int(r['knot'])>=0]
     if (path/'R20-grids.tsv').exists():d['R20']=rows(path/'R20-grids.tsv')
     for g in (240,480,960,1920):
      u=path/f'uncertainties-{g}.tsv'
      if u.exists():d[f'grid_{g}_rows']=len(rows(u))
     configurations.append(d)
    result['roles'][role]=dict(endpoint=values[-1],integrated_powers_trapezoid_erg=energies,HP_crossings=hp['crossings'],HP_steps=hp['steps'],HP_bridges=hp['bridges'],native_reference_errors=reference_errors,configurations=configurations,energy_scope='Original uniform-grid trapezoids for description, not fully qualified quadrature.',history_distinct_count=len({r.get('history_sha256') for r in configurations if 'history_sha256' in r}))
    with (out/f'curve-{role}.csv').open('w') as f:
     w=csv.DictWriter(f,fieldnames=values[0]);w.writeheader();w.writerows(values)
    if role=='source':source=values
    else:control=values
 with localcontext() as ctx:
  ctx.prec=110
  paired=[dict(t_s=a['t_s'],year=a['t_s']/YEAR,T_source=a['T'],T_control=b['T'],delta_T=a['T']-b['T'],U_delta=a['U_T']+b['U_T']) for a,b in zip(source,control)]
  result['paired_endpoint']=paired[-1]
  result['maximum_paired_HP_uncertainty_K']=max(v['U_delta'] for v in paired)
  with (out/'temperature.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=paired[0]);w.writeheader();w.writerows(paired)
  # Energy/C scale is a descriptive sensitivity diagnostic, not a new no-heating trajectory.
  heat=result['roles']['source']['integrated_powers_trapezoid_erg']['LH']
  beta=heat-result['roles']['source']['integrated_powers_trapezoid_erg']['DeltaLnu']
  result['chemical_temperature_scale_K']=dict(heating_energy_over_C_min=heat/min(v['C'] for v in source),heating_energy_over_C_max=heat/max(v['C'] for v in source),net_beta_energy_over_C_initial=beta/source[0]['C'],scope='Descriptive energy/heat-capacity scale; no independently qualified separated temperature response or counterfactual evolution is claimed.')
 save(out/'assessment.json',result)
if __name__=='__main__':assess(pathlib.Path(sys.argv[1]),int(sys.argv[2]))
