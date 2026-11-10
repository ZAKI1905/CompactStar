#!/usr/bin/env python3
"""Piecewise reference arithmetic through every observation, carrying its prefix."""
import concurrent.futures,functools,json,pathlib,sys,time
from decimal import Decimal as D,localcontext
from residuals import old,jacobian,MODEL
from high_precision import rows,lift
from run import save
@functools.lru_cache(maxsize=10)
def prefix(out,role,k):
 r=json.loads((out/f'roundoff/{role}-reference-{k}.json').read_text())
 return {D(0):(D(0),)*3,**{lift(v['time_s']):tuple(D(x) for x in v['accumulated_arithmetic_enclosure']) for v in r['records']}}
def point(task):
 out,role,k,p=task;old.jacobian_bound=jacobian
 with localcontext() as ctx:
  ctx.prec=110;before=prefix(out,role,k)
  if (p/'endpoint.txt').exists():
   t=lift((p/'endpoint.txt').read_text().split()[0]);u=before[t];rhs=(D(0),)*3
  else:
   row=rows(p/'maps.tsv')[0];t=lift(row['right']);initial=before[lift(row['left'])]
   r=old.audit(MODEL,role,f'reference-{k}-'+p.name,p,initial,'rk8pd',False);u=r['final_enclosure'];rhs=r['maximum_rhs_defect']
  return dict(role=role,reference=k,point=p.name,time_s=str(t),accumulated_arithmetic_enclosure=u,maximum_rhs_defect=rhs)
if __name__=='__main__':
 out=pathlib.Path(sys.argv[1]);workers=int(sys.argv[2]);start=time.monotonic();tasks=[]
 for role in ('source','control'):
  for k in range(5):
   paths=sorted((out/f'{role}-interpolation/reference-{k}').glob('point-*'));assert len(paths)>=1921
   tasks.extend((out,role,k,p) for p in paths)
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:results=list(pool.map(point,tasks,chunksize=8))
 target=out/'interpolation-roundoff';target.mkdir(exist_ok=False)
 for role in ('source','control'):
  for k in range(5):
   records=[r for r in results if r['role']==role and r['reference']==k]
   save(target/f'{role}-reference-{k}.json',dict(status='CONDITIONAL_PIECEWISE_ARITHMETIC_ENCLOSURE',records=records))
 save(out/'interpolation-roundoff-completion.json',dict(status='CONDITIONAL_PASS',observations=len(results),workers=workers,elapsed_s=time.monotonic()-start))
