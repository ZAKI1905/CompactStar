"""Score the committed rules; no fitted acceptance bands or production writes."""
import json,collections,resource,time
from pathlib import Path
import numpy as np
S=Path(__file__).resolve().parents[3];R=S.parents[1];OUT=R/'ba13-onset-pilot';DOC=S/'docs/validation/ba13-onset-pilot/evidence'
d=json.loads((OUT/'discovery.json').read_text());pred=json.loads((DOC/'predictions.json').read_text());ep=json.loads((DOC/'eos-predictions.json').read_text())
for p in (OUT/'analysis').glob('*.json'):
 case,i=p.stem.rsplit('-',1);d[f'{i}:{case}']=json.loads(p.read_text())
scores=[]
def add(rule,key,observed=None,band=None,eligible=True,reason=None):
 status='INCONCLUSIVE' if observed is None or band is None or not eligible else ('PASS' if observed<=band else 'FAIL')
 scores.append(dict(rule=rule,key=key,status=status,observed=observed,band=band,reason=reason))
def get(i,c,k):
 o=d.get(f'{i}:{c}')
 if o is None:return None
 if k.startswith('Z') and len(k)==3:return o['Z_nu'][int(k[1])][int(k[2])]
 if k in ['M','N_mu']:return o['star'][k]
 return o[k]
def paired(i,c,k):
 a,b=get(i,c,k),get(0,c,k);return None if a is None or b is None else a/b-1
for key,p in pred['absolute'].items():
 i,k=key.split(':');x=get(int(i),'r640000',k)
 add('P1',key,abs(x-p['discovery'][-1]) if x is not None else None,p['sharp'],p['eligible'])
for key,p in pred['paired'].items():
 i,k=key.split(':');x=paired(int(i),'r640000',k)
 add('P2-next',key,abs(x-p['discovery'][-1]) if x is not None else None,max(p['sharp'],1e-7),p['eligible'],'paired order eligibility is separate from stability')
 add('P2-tail',key,p['conditional_tail'],1e-7,p['eligible'])
 add('paired-stability-only',key,abs(x-p['discovery'][-1]) if x is not None else None,1e-7)
for c in ['r640000','nominal']+[f'r{n}' for n in range(80001,80006)]:
 for i in [10,20]:
  a,b=d.get(f'0:{c}'),d.get(f'{i}:{c}')
  if a is None or b is None:add('P3-L',f'{i}:{c}');add('P3-Tb',f'{i}:{c}');continue
  delta=b['L_nb']/a['L_nb']-b['L_nu']/a['L_nu']
  motion=(b['r_nb']-a['r_nb'])-(b['r_nu']-a['r_nu']);model=a['edge_sensitivity']*motion
  add('P3-L',f'{i}:{c}',abs(delta-model),.02*abs(delta))
  delta=b['Tb_node_factor']/a['Tb_node_factor']-b['Tb_factor']/a['Tb_factor'];model=(b['nu_b']-a['nu_b'])-(b['nu_node']-a['nu_node'])
  add('P3-Tb',f'{i}:{c}',abs(delta-model),.02*abs(delta))
for key,p in pred['secondary'].items():
 i,c,k=key.split(':');v=paired(int(i),c,k);add('P4-secondary',key,abs(v-p['prediction']) if v is not None else None,p['band'])
for key,p in ep.items():
 i,k=key.split(':');i=int(i);a=get(i,'nominal',k);b=get(i,'e16384c',k);c=get(i,'e32768c',k);f=get(i,'e65536c',k)
 add('P5-16384',key,abs(b-a) if b is not None else None,p['sharp'],p['eligible'])
 add('P5-32768-holdout',key,abs(c-b) if b is not None and c is not None else None,2*p['q']*abs(b-a) if b is not None else None,p['eligible'])
 add('P5-65536-consistency',key,abs(f-c) if c is not None and f is not None else None,2*p['q']*abs(c-b) if c is not None and b is not None else None,p['eligible'])
for i in [0,10,20]:
 for k in ['L_nu','C_nu','r_nu','R','M','N_mu']:
  a,b=get(i,'e32768c',k),get(i,'e32768d12',k)
  add('P5-depth',f'{i}:{k}',abs(b/a-1) if a is not None and b is not None else None,1e-12)
for k in ['L_nu','C_nu','r_nu','R','M','N_mu']:
 a,b,c,f=[get(20,case,k) for case in ['nominal','r640000','e16384c','cross']]
 add('P6',k,abs(f-b-c+a) if all(v is not None for v in [a,b,c,f]) else None,.1*max(abs(b-a),abs(c-a)) if all(v is not None for v in [a,b,c]) else None)
phases=['nominal']+[f'r{n}' for n in range(80001,80006)]
for k in ['L_nu','C_nu','r_nu']:
 for i in [0,10,20]:
  v=[get(i,c,k) for c in phases];add('P7-absolute',f'{i}:{k}',(max(v)-min(v))/abs(v[0]) if all(x is not None for x in v) else None,1e-9)
 for i in [10,20]:
  v=[paired(i,c,k) for c in phases];add('P7-paired',f'{i}:{k}',max(v)-min(v) if all(x is not None for x in v) else None,1e-10)
for i in [10,20]:
 # Check same target-reference index at each phase; distinct N indices do not
 # denote common physical nodes. The review's intended statement is testable
 # as equality of paired node drift across phases where target/reference share index.
 v=[]
 for c in phases:
  a,b=d.get(f'0:{c}'),d.get(f'{i}:{c}')
  if a and b and a['ib']==b['ib']:v.append(b['Tb_node_factor']/a['Tb_node_factor']-1)
 add('P7-Tb-shared-target-node',str(i),max(v)-min(v) if len(v)>1 else None,1e-10,reason=f'{len(v)} phases with reference and target on same index')
# Preregistered EOS availability and independent physical cross-check completion.
for n,cap in [(2048,11),(4096,11),(8192,11),(16384,11),(32768,11),(65536,11),(32768,12)]:
 p=OUT/f'table-{n}-{cap}.log'
 add('P5-order-and-identity',f'{n}:{cap}',0. if p.exists() and p.read_text().startswith('PASS') else None,0.)
hp=json.loads((OUT/'high-precision.json').read_text())
for j,row in enumerate(hp['pairs']):
 add('HP-order',f'{row["case"]}:{j}',0. if row['high_precision_order'] and row['rounded_order'] else 1.,0.)
# Full all-knot accounting and reconstruction separations.
metrics={'withheld_ratios':{},'discovery_reproduction':{},'thermal':{},'onset':{},'oracle':{},'tangents':{},'conditional_estimates':{}}
rev=json.loads((R/'ba13-uncertainty-independent-review-20261012/edge_coordinates.json').read_text());cv=json.loads((R/'ba13-uncertainty-independent-review-20261012/cstar_nu.json').read_text())
metrics['discovery_reproduction']={'L_relative_max':max(abs(d[k]['L_nu']/v['L_nu']-1) for k,v in rev.items()),'C_relative_max':max(abs(d[k]['C_nu']/v['nu_slope']-1) for k,v in cv.items())}
for key,o in d.items():
 i,case=key.split(':');base=d.get(f'0:{case}')
 therm=[]
 for j,t in enumerate(o['thermal']):
  row=dict(t,absolute_cache_bias=o['cache_factor'][j]-1,absolute_C_spatial=t['C_decache']/t['C_nu']-1,absolute_Tb=t['Tb_node']/t['Tb_nu']-1,absolute_Ts=t['Ts_prod']/t['Ts_nu']-1)
  if base:
   b=base['thermal'][j];row.update(paired_C_prod=t['C_prod']/b['C_prod']-1,paired_C_nu=t['C_nu']/b['C_nu']-1,paired_C_error=t['C_prod']/b['C_prod']-t['C_nu']/b['C_nu'],paired_Tb_error=t['Tb_node']/b['Tb_node']-t['Tb_nu']/b['Tb_nu'],paired_Ts_error=t['Ts_prod']/b['Ts_prod']-t['Ts_nu']/b['Ts_nu'])
  therm.append(row)
 metrics['thermal'][key]=therm
 metrics['onset'][key]={k:o[k] for k in ['r_nu','r_nb','mu_constancy','enthalpy_mu_disagreement','metric_secant_max','rb','rb_logeps','rb_pressure','r_on_pressure','cache_law_residual','nodal_vs_decache','Ts_node_reproduction']}
 metrics['onset'][key].update(absolute_L=o['L_nb']/o['L_nu']-1,B_residual=o['star']['residual']/o['star']['B_target'],B_bracket=o['star']['bracket']/o['star']['B_target'])
 if base:metrics['onset'][key].update(paired_L_prod=o['L_nb']/base['L_nb']-1,paired_L_nu=o['L_nu']/base['L_nu']-1,edge_motion_cm=(o['r_nu']-base['r_nu'])*1e5)
 if 'order48' in o:
  metrics['onset'][key]['quadrature_relative']={k:o['order48'][k]/o[k]-1 for k in ['L_nu','C_nu']}
  metrics['onset'][key]['composition']={k:o[k]/o['C_nu']-1 for k in ['C_trap','C_species_linear','C_nb_equilibrium']}
for i in [0,10,20]:
 for tol in ['2e-10','2e-12']:
  p=OUT/f'oracle-{i}-{tol}.log'
  if not p.exists():continue
  try:v=list(map(float,p.read_text().split()));assert len(v)==9
  except (ValueError,AssertionError):metrics['oracle'][f'{i}:{tol}']={'status':'REFUSED','message':p.read_text()};continue
  keys=['B','R','M','r_nu','L_nu','C_nu','rb','mu_inf','Hc'];o=dict(zip(keys,v));o['status']='COMPLETE';metrics['oracle'][f'{i}:{tol}']=o
  if f'{i}:r640000' in d:o['vs640_relative']={k:o[k]/get(i,'r640000',k)-1 for k in ['R','M','r_nu','L_nu','C_nu']}
for i in [0,10,20]:
 fine=metrics['oracle'].get(f'{i}:2e-12');coarse=metrics['oracle'].get(f'{i}:2e-10')
 if fine and coarse and fine['status']==coarse['status']=='COMPLETE':fine['tolerance_relative']={k:fine[k]/coarse[k]-1 for k in ['R','M','r_nu','L_nu','C_nu','mu_inf']}
 if i and fine and fine['status']=='COMPLETE':
  base=metrics['oracle'].get('0:2e-12')
  if base and base['status']=='COMPLETE':fine['paired']={k:fine[k]/base[k]-1 for k in ['R','M','r_nu','L_nu','C_nu','mu_inf']}
for i in [10,20]:
 if f'{i}:r640000' in d:
  r640=paired(i,'r640000','N_mu');r320=paired(i,'r320000','N_mu');dx=abs(d[f'{i}:nominal']['star']['B_target']/d['0:nominal']['star']['B_target']-1)
  metrics['tangents'][str(i)]={'paired_Nmu_radial_increment':r640-r320,'increment_per_B_fraction':(r640-r320)/dx,'status':'secant structure sensitivity only; no full tangent bound'}
for key,p in pred['absolute'].items():
 i,k=key.split(':');v=get(int(i),'r640000',k);last=p['increments'][-1]
 if v is not None and last:
  q=(v-p['discovery'][-1])/last
  metrics['withheld_ratios'][key]={'ratio':q,'effective_order':float(-np.log2(abs(q))) if q else None}
for k,p in pred['absolute'].items():metrics['conditional_estimates'][k]={'at320_relative':p['conditional_tail']/abs(p['discovery'][-1]) if p['conditional_tail'] is not None else None,'at640_relative':p['q']*p['conditional_tail']/abs(p['discovery'][-1]) if p['conditional_tail'] is not None else None,'eligible':p['eligible'],'condition':'future increments contract by frozen q; excludes common EOS/solver errors'}
for i in [0,10,20]:
 for tol in ['2e-10','2e-12']:
  row=metrics['oracle'].get(f'{i}:{tol}',{})
  add('independent-oracle-completion',f'{i}:{tol}',0. if row.get('status')=='COMPLETE' else None,0.)
# Independent validation cannot be inferred from same-star agreement. Preserve
# all oracle differences and inherited background diagnostics for adjudication.
res=json.loads((OUT/'resources.json').read_text());counts=collections.Counter(s['status'] for s in scores)
summary=dict(status='FAIL' if counts['FAIL'] else 'INCONCLUSIVE' if counts['INCONCLUSIVE'] else 'CONDITIONAL_PILOT_SCORES_PASS_NOT_GLOBAL_BOUNDS',counts=dict(counts),by_rule={rule:dict(collections.Counter(s['status'] for s in scores if s['rule']==rule)) for rule in sorted(set(s['rule'] for s in scores))},scores=scores,resources=res,missing_independent_bounds=['absolute common EOS/solver enclosure','uniform remainder beyond last refinement','inter-sample and inner-half response remainder','full tangent/background/tail paired bounds'])
(OUT/'scores.json').write_text(json.dumps(summary,indent=2)+'\n');(OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k in ['status','counts','by_rule']},indent=2))
u=resource.getrusage(resource.RUSAGE_SELF);(OUT/'score-resource.json').write_text(json.dumps(dict(cpu_seconds=u.ru_utime+u.ru_stime,max_rss_kib=u.ru_maxrss),indent=2)+'\n')
