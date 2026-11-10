#!/usr/bin/env python3
"""Compact receipts from full external evidence, with no qualification promotion."""
import concurrent.futures,json,pathlib,sys
from decimal import Decimal as D,localcontext
from run import ROOT,digest,save
out=pathlib.Path(sys.argv[1])
def read(p):return json.loads(p.read_text())
with localcontext() as ctx:
 ctx.prec=110
 assessment=read(out/'assessment-final.json');supplement=read(out/'supplement.json')
 assert assessment['all_computed_conditional_gates_pass'] and supplement['all_additional_conditional_gates_pass']
 with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:local=list(pool.map(read,sorted((out/'local-hp').glob('point-*.json'))))
 assert len(local)==1919 and all(r['status']=='CONDITIONAL_PASS' for r in local)
 chemical=dict(points=len(local),rhs_components_checked=sum(r['rhs_fidelity']['components_checked'] for r in local),maximum_rhs_difference_ulp=[max(D(r['rhs_fidelity']['maximum_difference_ulp'][j]) for r in local) for j in range(2)],rhs_limit_ulp=2,maximum_state_ratio=max(D(v) for r in local for v in r['chemical_state_ratios']),maximum_energy_ratio=max(D(r['chemical_energy_ratio']) for r in local),maximum_measured_thermal_rho=max(D(r['thermal_rho']) for r in local),maximum_refinement_ratio=max(D(c['ratio']) for r in local for x in r['convergence'] for c in x['checks'] if c['ratio'] is not None))
 hp={}
 for role in ('source','control'):
  vals=[read(out/f'high-precision/{role}-n{n}-p110.json') for n in (1,2,4,8)];p80=read(out/f'high-precision/{role}-n8-p80.json')
  detail={}
  for k in ('x','e','mu'):
   diffs=[max(abs(D(a[k])-D(b[k])) for a,b in zip(vals[j]['records'],vals[j+1]['records'])) for j in range(3)]
   detail[k]=dict(maximum_adjacent_differences=diffs,maximum_80_110_precision_difference=max(abs(D(a[k])-D(b[k])) for a,b in zip(p80['records'],vals[-1]['records'])),interpretation='Measured cross-check; analytic remainder, not agreement, supplies the enclosure.')
  hp[role]=detail
 resources={}
 for kind in ('reference','tier'):
  resources[kind]=[dict(name=p.parent.name,**read(p)) for p in sorted(out.glob(kind+'-*/completion.json'))]
 for name in ('sample','full','source-history','control-history','source-interpolation','control-interpolation','interpolation-roundoff','accepted'):
  p=out/f'{name}-completion.json';r=read(p);resources[name]={k:r[k] for k in ('wall_seconds','elapsed_s','exit_code','job','workers','observations') if k in r}
 for name in ('local-hp-sample','local-hp-full','thermal-enclosure-all'):
  r=read(out/(name+'.json'));resources[name]={k:r[k] for k in ('workers','elapsed_s','status')}
 roundoff={}
 for role in ('source','control'):
  roundoff[role]={}
  for name in ['main']+[f'reference-{i}' for i in range(5)]:
   r=read(out/f'roundoff/{role}-{name}.json');old=read(out/f'roundoff-before-box-checks/{role}-{name}.json')
   assert r['final_enclosure']==old['final_enclosure']
   roundoff[role][name]={k:r[k] for k in ('final_enclosure','maximum_rhs_defect','maximum_replay_difference','elapsed_s')}
 result=dict(status='PRELIMINARY — NOT NUMERICALLY QUALIFIED',all_computed_conditional_gates_pass=True,formal_Stage_A_qualified=False,Stage_B_launched=False,source_evidence_recovered=True,source_packets=43248,source_unique_reconstructions=1921,source_strict_interiors=1919,chemical=chemical,high_precision_cross_checks=hp,roundoff=roundoff,resources=resources,retained_native_source_trajectories_rerun=False,independent_HP_original_t0_references_computed=True,notes=['All new analytic and arithmetic verification is an unratified scientific candidate under GOVERNANCE section 5.','No event crosses a Cstar knot in either one-year trajectory.','Three control tiers ran concurrently with 24 physical workers each; local high-precision chemistry used 48 workers.','Slurm accounting is disabled; wrapper elapsed times are measured, aggregate peak RSS and exact total CPU usage are unavailable.'],code_sha256={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'tools/linux/preliminary_pair').iterdir()) if p.suffix in ('.py','.cpp','.sbatch')})
 save(out/'summary.json',result)
 print('PASS compact evidence assembled; formal qualification remains false')
