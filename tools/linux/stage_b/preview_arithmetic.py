#!/usr/bin/env python3
"""Early read-only role audit; final assessment independently repeats all roles."""
import concurrent.futures,functools,json,os,pathlib,sys,time
from decimal import Decimal as D
from piecewise import ROOT,lift
from apply_residuals import audit
import observation_residuals as observations
from run import save

out=pathlib.Path(sys.argv[1]);role=sys.argv[2];workers=int(sys.argv[3])
assert role in ('source','control')
target=out/(role+'-arithmetic-preview');target.mkdir(exist_ok=False)
start=time.monotonic()
def reference(k):
 r=audit(role,f'reference-{k}',out/f'{role}-supplement/reference-{k}',method='rk8pd' if k<3 else 'rkf45')
 save(target/f'reference-{k}.json',r)
 return dict(reference=k,status=r['status'],final_enclosure=r['final_enclosure'])
@functools.lru_cache(maxsize=5)
def prefix(directory,selected_role,k):
 assert directory==out and selected_role==role
 r=json.loads((target/f'reference-{k}.json').read_text())
 return {D(0):(D(0),)*3,**{lift(v['time_s']):tuple(D(x) for x in v['accumulated_arithmetic_enclosure']) for v in r['records']}}
observations.prefix=prefix  # Identical recurrence; preview has its own input/output directory.
with concurrent.futures.ProcessPoolExecutor(max_workers=min(workers,5)) as pool:references=list(pool.map(reference,range(5)))
tasks=[]
for k in range(5):
 paths=sorted((out/f'{role}-interpolation/reference-{k}').glob('point-*'));assert len(paths)>=1921
 tasks.extend((out,role,k,p) for p in paths)
with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:records=list(pool.map(observations.point,tasks,chunksize=8))
for k in range(5):save(target/f'observations-{k}.json',[r for r in records if r['reference']==k])
save(target/'completion.json',dict(status='CONDITIONAL_ARITHMETIC_PREVIEW',role=role,references=references,observations=len(records),workers=workers,elapsed_s=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID'),scope='Read-only early check using the same actual-stage and observation recurrence as the final full-matrix assessment. It does not qualify state/power budgets or replace the final audit.'))
