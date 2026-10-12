#!/usr/bin/env python3
"""Per-target center/tail and independent fixed-background response diagnostics."""
import json,sys,time
from pathlib import Path
import numpy as np
S=Path(__file__).resolve().parents[3];R=S.parents[1]/'ba13-uncertainty';OLD=S.parents[1]/'ba13-recovery'
sys.path.insert(0,str(S/'tests/analysis'))
from chemical_trackr_budget import Model,matrix,integrate,zsolve,zerror,tail
records=[]
for i in range(21):
 p=OLD/'run1/raw/chemical-profiles'/f'target-{i}';d=np.loadtxt(p/'profile.tsv',skiprows=1);model=Model(p)
 c=json.loads((p/'certificate.json').read_text());r,m,nu,nb=d[0,:4]
 center=np.abs(matrix(model.response(np.array([nb]))[0]))*4*np.pi/3*r**3*np.exp(-nu)*1e54
 inherited=np.array(c['center_error']);positive=inherited>0
 record=dict(index=i,center_current_profile_estimate=center.tolist(),center_inherited=inherited.tolist(),
  center_estimate_within_inherited=bool(np.all(center<=inherited)),
  max_center_ratio=float(np.max(center[positive]/inherited[positive])),
  background_inherited=c['background_error'],tail_current=c['tail_error'],
  center_status='ESTIMATE_ONLY_NO_INDEPENDENT_CENTER_REMAINDER',background_status='UNBOUNDED_PER_TARGET_REFINEMENT',
  tail_status='SAME_CUT_ANALYTIC_COMPARISON_PRESERVED',integrals=None)
 if i in (0,10,20):
  ladder={};start=time.monotonic()
  for order in (16,32,64,96):
   g=matrix(integrate(d,model,order)[0]);ladder[str(order)]=dict(G=g.tolist(),Z=zsolve(g).tolist())
  g=matrix(integrate(d,model,32,partition_level=1)[0]);ladder['32bisect']=dict(G=g.tolist(),Z=zsolve(g).tolist())
  record['integrals']=ladder;record['seconds']=time.monotonic()-start
 records.append(record)
 print('chemical target',i,'center ratio',record['max_center_ratio'],flush=True)
(R/'chemical-diagnostic.json').write_text(json.dumps(dict(status='PILOT_PARTIAL_BOUNDS_ONLY',records=records),indent=2)+'\n')
