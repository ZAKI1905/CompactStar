#!/usr/bin/env python3
"""Exact candidate CPU/runtime/arithmetic class comparison, not model-label inference."""
import argparse
import json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('receipts',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
rows=[]
for path in a.receipts:
 d=json.loads(path.read_text())
 features=d['cpu']['flags'];assert len(features)==1,'heterogeneous CPUs inside node'
 relevant=sorted(f for f in features[0].split() if f.startswith(('avx','sse','fma','xsave')) or f in ('f16c','osxsave','fxsr'))
 libs={Path(k).name:v for k,v in d['loaded_files'].items() if Path(k).name.startswith(('libm.so','libc.so','ld-linux','libstdc++','libgcc_s'))}
 numeric='\n'.join(s for s in d['raw_probe'].splitlines() if not s.startswith('loaded\t'))
 rows.append(dict(path=str(path),host=d['host'],types=d['types'],features=relevant,loaded=libs,dispatch=d['dispatch_environment'],numeric=numeric,compiler=d['compiler']))
assert len({r['host'] for r in rows})>=3,'three representative nodes required'
for key in ('types','features','loaded','dispatch','numeric','compiler'):
 assert all(r[key]==rows[0][key] for r in rows), 'runtime/arithmetic disagreement: '+key
a.output.write_text(json.dumps(dict(status='PASS_EXACT_CANDIDATE_RUNTIME_CLASS',nodes=sorted({r['host'] for r in rows}),reference=rows[0],receipts=[r['path'] for r in rows]),indent=2)+'\n')
