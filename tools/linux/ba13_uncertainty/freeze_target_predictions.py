#!/usr/bin/env python3
import csv,json,math
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=S.parents[1]/'ba13-uncertainty'
records=[]
for i in (5,15):
 for key in ('R','M','N_n','N_e','N_mu'):
  values=[float(next(csv.DictReader((R/'runs'/f'{c}-{i}'/'star.tsv').open(),delimiter='\t'))[key]) for c in ('r20000','r40000','nominal','r160000')]
  d=[abs(values[j+1]-values[j]) for j in range(3)];q=[d[j+1]/d[j] for j in range(2)];orders=[-math.log2(z) for z in q]
  eligible=all(0<z<=.75 for z in q) and min(orders)>0 and max(orders)<=2*min(orders) and min(d)>128*2.220446049250313e-16*max(map(abs,values))
  records.append(dict(index=i,axis='radial',table='star.tsv',quantity=key,knot=0,values=values,differences=d,ratios=q,orders=orders,estimator_eligible=eligible,finest_remainder_prediction=2*d[-1]/(1-max(q)) if eligible else None,nominal_error_prediction=abs(values[-1]-values[-2])+2*d[-1]/(1-max(q)) if eligible else None))
with (R/'target-holdout-predictions.json').open('x') as f:json.dump(records,f,indent=2);f.write('\n')
print('eligible',sum(v['estimator_eligible'] for v in records),'of',len(records))
