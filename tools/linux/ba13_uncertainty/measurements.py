#!/usr/bin/env python3
import csv,json
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')))
def table(p):return list(csv.DictReader(p.open(),delimiter='\t'))
rows=[]
for p in sorted((R/'runs').glob('*')):
 if not (p/'COMPLETE.json').exists():continue
 case,i=p.name.rsplit('-',1);i=int(i)
 for tab,keys in [('star.tsv',['R','M','N_n','N_e','N_mu','rho_c']),('thermal.tsv',['Cstar','Tsurface','Tsurface_inf','g14'])]:
  nom=table(R/'runs'/f'nominal-{i}'/tab);now=table(p/tab)
  b0=table(R/'runs'/'nominal-0'/tab);pair=R/'runs'/f'{case}-0'/tab
  pairrows=table(pair) if pair.exists() else None
  for k,row in enumerate(now):
   for key in keys:
    v=float(row[key]);v0=float(nom[k][key]);paired=None
    if pairrows:paired=v/float(pairrows[k][key])-v0/float(b0[k][key])
    rows.append(dict(case=case,index=i,table=tab,knot=k,quantity=key,value=v,nominal_value=v0,relative_change=(v-v0)/abs(v0),paired_depletion_ratio_change=paired,bound_status='MEASURED_SPREAD_NOT_BOUND'))
with (R/'measurements.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
ints=json.loads((R/'integral-diagnostic.json').read_text());out=[]
for i in (10,20):
 for key in ('Cstar_slope','mu_B','P0','P1','Ltilde_Me','Ltilde_Mmu'):
  values={}
  for case in ('nominal','r160000','r320000'):
   v=next(v[key] for v in ints if v['case']==case and v['index']==i and v['order']==32 and not v['bisect'])
   b=next(v[key] for v in ints if v['case']==case and v['index']==0 and v['order']==32 and not v['bisect'])
   values[case]=v/b-1
  out.append(dict(index=i,quantity=key,depletion_ratios=values,change_finest_minus_nominal=values['r320000']-values['nominal'],status='PAIRED_DIAGNOSTIC_NOT_JOINT_ERROR_BOUND'))
(R/'paired-integrals.json').write_text(json.dumps(out,indent=2)+'\n')
print('measurement rows',len(rows))
