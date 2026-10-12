#!/usr/bin/env python3
"""Freeze pilot predictions before examining withheld resolution levels."""
import csv,json,math
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')));records=[]
for i in (0,10,20):
 for axis,cases in [('radial',['r20000','r40000','nominal','r160000']),('eos',['e2048','e4096','nominal','e16384'])]:
  for table,keys in [('star.tsv',['R','M','N_n','N_e','N_mu']),('thermal.tsv',['Cstar','Tsurface','Tsurface_inf','g14'])]:
   tables=[list(csv.DictReader((R/'runs'/f'{c}-{i}'/table).open(),delimiter='\t')) for c in cases]
   for row in range(len(tables[0])):
    for key in keys:
     v=[float(t[row][key]) for t in tables];d=[abs(v[j+1]-v[j]) for j in range(3)]
     q=[d[j+1]/d[j] if d[j] else None for j in range(2)]
     orders=[-math.log2(z) if z and z>0 else None for z in q]
     eligible=all(z is not None and 0<z<=.75 for z in q) and min(orders)>0 and max(orders)<=2*min(orders)
     # Require increments resolved above conservative arithmetic floor.
     eligible=eligible and min(d)>128*2.220446049250313e-16*max(map(abs,v))
     tail=2*d[-1]/(1-max(q)) if eligible else None
     records.append(dict(index=i,axis=axis,table=table,knot=row,quantity=key,values=v,differences=d,ratios=q,orders=orders,
       estimator_eligible=eligible,finest_remainder_prediction=tail,nominal_error_prediction=abs(v[3]-v[2])+tail if eligible else None,
       status='CONDITIONAL_EMPIRICAL_HOLDOUT_REQUIRED' if eligible else 'INCONCLUSIVE'))
with (R/'holdout-predictions.json').open('x') as f:json.dump(records,f,indent=2);f.write('\n')
print('eligible',sum(r['estimator_eligible'] for r in records),'of',len(records))
