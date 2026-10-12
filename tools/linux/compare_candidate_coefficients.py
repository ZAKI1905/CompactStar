#!/usr/bin/env python3
"""ADR-0020 componentwise agreement conditional on fresh refinement validation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--mac',type=Path,required=True);p.add_argument('--linux',type=Path,required=True);p.add_argument('--validation',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
mac,linux=[json.loads(f.read_text()) for f in (a.mac,a.linux)]
v=json.loads(a.validation.read_text());assert v['PASS'] and v['repeat_requested'] and len(set(v['candidate_sha256']))==1
assert hashlib.sha256(a.linux.read_bytes()).hexdigest()==v['candidate_sha256'][0], 'candidate differs from independently repeated artifact'
assert linux['candidate_only'] is True and not linux.get('governed_baseline',False)
assert mac['units']==linux['units'] and mac['fixture']==linux['fixture']
def flat(x):return [v for y in x for v in flat(y)] if isinstance(x,list) else [float(x)]
components={}
for key in ('G','Q','Z','I','W'):
 arrays=[flat(d[k][key]) for k in ('values','numerical_error') for d in (mac,linux)]
 assert len(set(map(len,arrays)))==1
 rows=[]
 for x,y,ex,ey in zip(*arrays):
  assert all(math.isfinite(z) for z in (x,y,ex,ey)), 'nonfinite coefficient or uncertainty'
  assert ex>=0 and ey>=0
  error=abs(x-y);bound=ex+ey
  rows.append(dict(mac=x,linux=y,absolute=error,relative=error/max(abs(x),abs(y)) if x or y else 0.,combined_uncertainty=bound,utilization=error/bound if bound else None,pass_candidate=error<=bound))
 components[key]=rows
passed=all(r['pass_candidate'] for rs in components.values() for r in rs)
result=dict(status='PASS_CANDIDATE_COMPONENT_AGREEMENT' if passed else 'FAIL',not_owner_ratified=True,
    scope='Structure-1 static coefficients; empirical refinement uncertainty, not rigorous enclosure or native Phase-6 frozen certificate',
    input_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (a.mac,a.linux,a.validation)},components=components,
    provenance_note='Actual regenerated onset/provider/EOS provenance is preserved in native artifact; no string substitution')
a.output.write_text(json.dumps(result,indent=2)+'\n')
assert passed
