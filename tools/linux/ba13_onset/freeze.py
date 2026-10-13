"""Freeze deterministic predictions from authenticated DISCOVERY JSON only."""
import json,sys
from pathlib import Path
S=Path(__file__).resolve().parents[3]; R=S.parents[1]
rev=R/'ba13-uncertainty-independent-review-20261012'
e=json.loads((rev/'edge_coordinates.json').read_text());c=json.loads((rev/'cstar_nu.json').read_text());z=json.loads((rev/'z_radial.json').read_text())
cases=['r20000','r40000','nominal','r160000','r320000'];eps=2.220446049250313e-16
def predict(v):
 scale=max(map(abs,v));d=[b-a for a,b in zip(v,v[1:])];floor=128*eps*max(scale,1e-300)
 ratios=[abs(b/a) if a else None for a,b in zip(d,d[1:])]
 q=max([x for x in ratios[-2:] if x is not None],default=1.)
 eligible=all(abs(x)>floor for x in d[-3:]) and all(d[j]*d[j-1]>0 for j in range(len(d)-2,len(d))) and all(x is not None and 0<x<.5 for x in ratios[-2:])
 return dict(discovery=v,increments=d,ratios=ratios,arithmetic_floor=floor,eligible=eligible,q=q,prediction=v[-1]+d[-1]/3,sharp=2*q*abs(d[-1]),loose=2*abs(d[-1])/(1-q) if q<1 else None,conditional_tail=2*q*abs(d[-1])/(1-q) if q<1 else None)
def val(i,k,field):
 if field in ['L_nu','r_nu']:return e[f'{i}:{k}'][field]
 if field=='C_nu':return c[f'{i}:{k}']['nu_slope']
 a,b=map(int,field[1:]);return z[f'{i}:{k}']['Z_nu'][a][b]
if __name__=='__main__':
 out={'provenance':'review discovery; not pilot evidence','absolute':{},'paired':{},'secondary':{}}
 for i in [0,10,20]:
  for field in ['L_nu','r_nu','C_nu','Z00','Z01','Z11']:out['absolute'][f'{i}:{field}']=predict([val(i,k,field) for k in cases])
 for i in [10,20]:
  for field in ['L_nu','C_nu']:out['paired'][f'{i}:{field}']=predict([val(i,k,field)/val(0,k,field)-1 for k in cases])
 # Physical targets are equally spaced to represented precision; read exact bytes.
 import csv
 rows=list(csv.DictReader((R/'ba13-recovery/run1/raw/target_stars.tsv').open(),delimiter='\t'))
 print('target columns',list(rows[0]))
 # Historical fraction column is recorded abscissa; physical B ratio is frozen too.
 Bkey=next(k for k in rows[0] if k in ['B_target','target_B'])
 b=[float(r[Bkey]) for r in rows];x=[(v-b[0])/b[0] for v in b]
 for i in [5,15]:
  for field in ['L_nu','C_nu']:
   for k in cases:
    p=x[i]/x[20]*(val(20,k,field)/val(0,k,field)-1)
    out['secondary'][f'{i}:{k}:{field}']={'prediction':p,'band':1e-9+1e-3*abs(p),'x_ratio':x[i]/x[20]}
 (S/'docs/validation/ba13-onset-pilot/evidence/predictions.json').write_text(json.dumps(out,indent=2)+'\n')
