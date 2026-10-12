#!/usr/bin/env python3
"""Report captured structural drift without converting an ODE rtol into a bound."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--reference',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
def read(p):
 with p.open() as f:return list(csv.DictReader((line for line in f if line.strip() and not line.startswith('#')),delimiter='\t'))
x,y=read(a.reference),read(a.candidate);assert len(x)==len(y) and x and set(x[0])==set(y[0])
if any(u.get('target_M_Msun')!=v.get('target_M_Msun') for u,v in zip(x,y)): raise ValueError('target identities changed')
selections={k:all(float(u[k])==float(v[k]) for u,v in zip(x,y)) for k in ('epsilon_c_gcm3','ec_km^-2') if k in x[0]}
report={}
for key in x[0]:
 rows=[]
 for i,(u,v) in enumerate(zip(x,y)):
  try:old,new=float(u[key]),float(v[key])
  except ValueError:
   assert u[key]==v[key],(i,key,u[key],v[key]);continue
  if not (math.isfinite(old) and math.isfinite(new)): raise ValueError('nonfinite structural capture')
  rows.append(dict(row=i,reference=old,candidate=new,absolute=abs(old-new),relative=abs(old-new)/max(abs(old),abs(new)) if old or new else 0))
 if rows:report[key]=dict(max_absolute=max(r['absolute'] for r in rows),max_relative=max(r['relative'] for r in rows),changed=sum(r['reference']!=r['candidate'] for r in rows))
a.output.write_text(json.dumps(dict(status='DIAGNOSTIC_NO_GLOBAL_ERROR_BOUND',qualification='UNQUALIFIED until quantity-specific propagated/refinement uncertainty is supplied; local ODE tolerance alone is insufficient',central_selection_payload_equal=selections,inputs={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (a.reference,a.candidate)},columns=report),indent=2)+'\n')
