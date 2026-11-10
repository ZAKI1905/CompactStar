#!/usr/bin/env python3
"""Verify every original packet maps to the recovered unique checkpoint, no inference by sampling."""
import concurrent.futures,json,pathlib,sys,time
from run import BASE,digest,save
out=pathlib.Path(sys.argv[1]);source=BASE/'stage-a-parallel-20261110-r1/source-stage/run/trajectory/source'
start=time.monotonic();canonical={}
ignored={'checkpoint_index','main_rtol','main_atol','main_maximum_step_s'}
for i in range(1921):
 p=next(p for p in (out/'sample'/f'point-{i}',out/'full'/f'point-{i}') if p.exists())
 assert (p/'consistency.txt').read_text().startswith('PASS')
 metadata=json.loads((p/'packet/checkpoint.json').read_text())
 canonical[i]=({name:digest(p/'packet'/name) for name in ('state.tsv','diagnostics.tsv')},{k:v for k,v in metadata.items() if k not in ignored})

def check(tier):
 failures=[];count=0
 for n in (240,480,960,1920):
  for i in range(n+1):
   index=i*(1920//n);expected,meta=canonical[index];p=tier/f'checkpoints-{n}/point-{i}'
   for name,h in expected.items():
    if digest(p/name)!=h:failures.append(str(p/name))
   actual=json.loads((p/'checkpoint.json').read_text())
   if {k:v for k,v in actual.items() if k not in ignored}!=meta:failures.append(str(p/'checkpoint.json'))
   count+=1
 return dict(tier=tier.name,packets=count,failures=failures,main_sha256=digest(tier/'main.accepted_states.tsv'))

with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:results=list(pool.map(check,sorted(source.glob('tier-*'))))
assert sum(r['packets'] for r in results)==43248
assert len({r['main_sha256'] for r in results})==1
result=dict(status='PASS' if all(not r['failures'] for r in results) else 'FAIL',original_packets=43248,unique_reconstructed_packets=1921,records=results,ignored_metadata_fields=sorted(ignored),elapsed_s=time.monotonic()-start,meaning='Every retained state/diagnostic packet byte-matches the unique reconstruction at its exact time; remaining metadata matches after excluding only tier/cap/index fields.')
save(out/'coverage.json',result);assert result['status']=='PASS'
