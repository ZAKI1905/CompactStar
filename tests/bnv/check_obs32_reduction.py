#!/usr/bin/env python3
"""Bounded read-only negative checks; run from the repository root. No ODEs."""
import importlib.util,tempfile,json,gzip,shutil
from pathlib import Path
sp=importlib.util.spec_from_file_location('reduce_obs32','tests/bnv/analyze_obs32_forensics.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
base=Path('docs/validation/evidence/obs32-forensics');root=base/'runs/primary/output';supp=base/'runs/supplement/output';results={}
with tempfile.TemporaryDirectory(prefix='obs32-readonly-negative-') as tmp:
 tmp=Path(tmp);ev=tmp/'evidence';shutil.copytree(base,ev)
 w=json.loads((ev/'saved-witness.json').read_text());w['D_O1_x']='1e-12';(ev/'saved-witness.json').write_text(json.dumps(w))
 try:m.analyze(root,ev,supp)
 except AssertionError:results['incorrect_archived_budget_rejected']=True
 else:raise RuntimeError('mutant budget accepted')
 tr=tmp/'traces';shutil.copytree(root,tr);p=tr/'unsplit-O1.rhs.tsv.gz';s=gzip.decompress(p.read_bytes()).decode().splitlines();fields=s[1].split('\t');fields[-1]='0';s[1]='\t'.join(fields);p.write_bytes(gzip.compress(('\n'.join(s)+'\n').encode(),mtime=0))
 try:m.analyze(tr,base,supp)
 except AssertionError:results['noncontaining_cell_rejected']=True
 else:raise RuntimeError('mutant cell accepted')
 results['compressed_reduction_equals_saved']=m.analyze(root,base,supp)==json.loads((base/'analysis.json').read_text())
 results['ODE_integrations']=0
 print(json.dumps(results,indent=2))
