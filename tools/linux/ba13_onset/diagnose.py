"""Post-outcome diagnosis from the fixed cases. No new levels, bands or solves."""
import json,resource
from pathlib import Path
import numpy as np
S=Path(__file__).resolve().parents[3];R=S.parents[1];out=R/'ba13-onset-pilot';metrics=json.loads((out/'metrics.json').read_text());scores=json.loads((out/'scores.json').read_text());d=json.loads((out/'discovery.json').read_text())
for p in (out/'analysis').glob('*.json'):
 c,i=p.stem.rsplit('-',1);d[f'{i}:{c}']=json.loads(p.read_text())
print('FAILURES')
for row in scores['scores']:
 if row['status']=='FAIL':print(row)
print('RATIOS',json.dumps(metrics['withheld_ratios']))
print('ORACLES',json.dumps(metrics['oracle'],indent=1))
print('CONDITIONAL ESTIMATES',json.dumps(metrics['conditional_estimates'],indent=1))
ph=['nominal']+[f'r{n}' for n in range(80001,80006)]
diag={'status':'POST_OUTCOME_DIAGNOSTIC_NOT_NEW_VALIDATION','phase':{},'eos':{},'thermal':{},'nominal':{},'absolute_oracle':{},'conditioning':{},'Z':{}}
for i in [0,10,20]:
 a=d[f'{i}:nominal'];arr=[]
 for c in ph:
  o=d[f'{i}:{c}'];row={'case':c,'B_residual':o['star']['residual']/o['star']['B_target'],'B_bracket':o['star']['bracket']/o['star']['B_target']}
  row.update({k:o[k]/a[k]-1 for k in ['L_nu','C_nu','r_nu','mu_inf','R']})
  row.update({k:o['star'][k]/a['star'][k]-1 for k in ['M','N_mu','rho_c']})
  row.update(mu_constancy=o['mu_constancy'],anchor_change=o['mu_inf_altref']/o['mu_inf']-1)
  arr.append(row)
 diag['phase'][str(i)]=arr
 diag['eos'][str(i)]=[{**{'case':c},**{k:d[f'{i}:{c}'][k]/a[k]-1 for k in ['L_nu','C_nu','r_nu','mu_inf']},'B_residual':d[f'{i}:{c}']['star']['residual']/d[f'{i}:{c}']['star']['B_target']} for c in ['e16384c','e32768c','e65536c','e32768d12']]
 diag['nominal'][str(i)]=metrics['onset'][f'{i}:nominal'];diag['thermal'][str(i)]=dict(first=metrics['thermal'][f'{i}:nominal'][0],cache_bias_max=max(v['absolute_cache_bias'] for v in metrics['thermal'][f'{i}:nominal']),max_paired_C_error=max(abs(v['paired_C_error']) for v in metrics['thermal'][f'{i}:nominal']))
 fine=metrics['oracle'][f'{i}:2e-12'];coarse=metrics['oracle'][f'{i}:2e-10'];diag['absolute_oracle'][str(i)]=dict(nominal_relative={k:a[k]/fine[k]-1 for k in ['L_nu','C_nu','r_nu','R','mu_inf','rb']},fine_relative=fine['vs640_relative'],tolerance_relative=fine['tolerance_relative'])
 if i:
  b=metrics['oracle']['0:2e-12'];base=d['0:nominal'];diag['absolute_oracle'][str(i)]['paired_nominal_discrepancy']={k:a[k]/base[k]-fine[k]/b[k] for k in ['L_nu','C_nu','r_nu','R','mu_inf','rb']}
 diag['Z'][str(i)]={k:d[f'{i}:r640000'][k] for k in ['Z_nu','Z_prod','G_nu','G_prod']}
 diag['Z'][str(i)]['condition_G']=float(np.linalg.cond(np.array(d[f'{i}:r640000']['G_nu'])))
# Analyze actual fixed EOS conditioning from authenticated original table, not an
# independent interpolant or additional refinement experiment.
table=np.loadtxt(R/'ba13-uncertainty/runs/nominal-0/profile/freegas.tsv',skiprows=1)
n=table[:,2];p=table[:,1];good=(n[1:]!=n[:-1])&(p[1:]>p[:-1]);gamma=np.diff(np.log(p))[good]/np.diff(np.log(n))[good];mid=np.sqrt(n[:-1]*n[1:])[good];t=7.3567289037328318e-9
mask=(mid>1.001*t)&(mid<7*t);diag['conditioning']=dict(soft_layer_gamma_min=float(np.min(gamma[mask])),soft_layer_gamma_max=float(np.max(gamma[mask])),note='finite log secants, not a proved lower bound')
print('DIAGNOSTIC',json.dumps(diag,indent=1))
(out/'diagnosis.json').write_text(json.dumps(diag,indent=2)+'\n');u=resource.getrusage(resource.RUSAGE_SELF);(out/'diagnosis-resource.json').write_text(json.dumps(dict(cpu_seconds=u.ru_utime+u.ru_stime,max_rss_kib=u.ru_maxrss),indent=2)+'\n')
