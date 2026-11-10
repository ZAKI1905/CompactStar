#!/usr/bin/env python3
"""Compact provenance from retained launch/build/allocation receipts."""
import json,pathlib,re,sys
from run import save,digest
out=pathlib.Path(sys.argv[1]);roots=[pathlib.Path(p) for p in sys.argv[2:]];result=dict(allocations=[],builds=[],launches=[],completions=[])
for root in roots:
 for p in sorted(root.glob('allocation-*.txt')):
  s=p.read_text();fields={}
  for k in ('JobId','ArrayJobId','ArrayTaskId','JobName','JobState','RunTime','TimeLimit','Partition','NodeList','NumCPUs','CPUs/Task','MinMemoryNode','SubmitTime','StartTime'):
   match=re.search(r'(?<!\S)'+re.escape(k)+r'=([^\s]+)',s)
   if match:fields[k]=match.group(1)
  result['allocations'].append(dict(path=str(p),sha256=digest(p),fields=fields))
 for p in sorted(root.glob('*build*/build.json')):
  result['builds'].append(dict(path=str(p),sha256=digest(p),bindings=json.loads(p.read_text())))
 for p in sorted(root.glob('**/*invocation.json')):
  result['launches'].append(dict(path=str(p),sha256=digest(p),record=json.loads(p.read_text())))
 for p in sorted(root.glob('**/*completion.json')):
  result['completions'].append(dict(path=str(p),sha256=digest(p),record=json.loads(p.read_text())))
result['scope']='Preserved Slurm allocation and executable/input receipts. Slurm accounting was unavailable; these are allocations and measured wall times, not invented CPU accounting. Physical workers: 24 per reconstruction task, requiring 48 allocated SMT threads on these 48-core/96-thread nodes. Reference trajectories used one worker each; independent HP and assessment work is recorded separately.'
save(out/'resources.json',result)
