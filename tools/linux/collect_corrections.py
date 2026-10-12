#!/usr/bin/env python3
"""Collect bounded run receipts, compact text evidence, exact identities; no promotion."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
p=argparse.ArgumentParser();p.add_argument('--toolchain',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--jobs',default='14463,14464,14465,14466,14467,14468,14469,14470,14471');a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=False)
root=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
summary=dict(status='CANDIDATE_ONLY',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),tests={},artifacts={},jobs={})
for path in sorted((a.toolchain/'corrections').rglob('*')):
 if not path.is_file():continue
 rel=path.relative_to(a.toolchain/'corrections')
 # Record every run-local file identity; copy only compact evidence, no build
 # products or scientific input tree. Producer binaries are identified, not imported.
 if path.suffix in ('.json','.xml','.log','.txt','.tsv','.exit'):
  summary['artifacts'][str(rel)]={'sha256':sha(path),'bytes':path.stat().st_size}
 if path.suffix=='.xml':
  try:tree=ET.parse(path)
  except ET.ParseError:continue
  summary['tests'][str(rel)]=[dict(name=t.get('name'),seconds=t.get('time'),status='FAIL' if t.find('failure') is not None or t.get('status')=='fail' else 'SKIP' if t.find('skipped') is not None else 'PASS') for t in tree.findall('.//testcase')]
 # Do not archive bulky repetition tables or scientific profile/input trees.
 if path.suffix not in ('.json','.xml','.log','.exit'):continue
 if any(x in rel.parts for x in ('thermal','central','tangent-sequence')):continue
 if path.stat().st_size>4_000_000:continue
 dest=a.output/'runs'/rel
 dest=dest.with_name(dest.name+'.gz');dest.parent.mkdir(parents=True,exist_ok=True)
 with dest.open('wb') as f:
  with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as out:out.write(path.read_bytes())
ids=a.jobs
res=subprocess.run(['sacct','-j',ids,'--units=K','--parsable2','--format=JobID,JobName,State,ExitCode,Elapsed,TotalCPU,MaxRSS,AllocCPUS,NodeList'],capture_output=True,text=True)
summary['sacct_returncode']=res.returncode
(a.output/'sacct.txt').write_text(res.stdout+res.stderr)
(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(a.output/'SHA256SUMS').write_text(''.join(sha(f)+'  '+str(f.relative_to(a.output))+'\n' for f in sorted(a.output.rglob('*')) if f.is_file()))
print(a.output)
