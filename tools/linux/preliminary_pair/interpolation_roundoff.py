#!/usr/bin/env python3
"""Actual reference observation arithmetic including the accumulated prefix."""
import concurrent.futures,functools,json,pathlib,sys,time
from decimal import Decimal as D,localcontext
from high_precision import rows,lift,save,digest
from roundoff import audit

@functools.lru_cache(maxsize=10)
def prefix(out,role,k):
 a=json.loads((out/f'roundoff/{role}-reference-{k}.json').read_text())
 return {D(0):(D(0),)*3,**{lift(r['time_s']):tuple(D(v) for v in r['accumulated_arithmetic_enclosure']) for r in a['records']}}

def point(args):
 out,role,k,i=args
 with localcontext() as ctx:
  ctx.prec=110;p=out/f'{role}-interpolation/reference-{k}/point-{i}';before=prefix(out,role,k)
  if (p/'endpoint.txt').exists():
   t=lift((p/'endpoint.txt').read_text().split()[0]);u=before[t];rhs=(D(0),)*3
  else:
   row=rows(p/'maps.tsv')[0];t=lift(row['right']);initial=before[lift(row['left'])]
   a=audit(out,role,f'reference-{k}-point-{i}',p,initial,'rk8pd',False)
   u=a['final_enclosure'];rhs=a['maximum_rhs_defect']
  return dict(role=role,reference=k,point=i,time_s=t,accumulated_arithmetic_enclosure=u,maximum_rhs_defect=rhs)

if __name__=='__main__':
 out=pathlib.Path(sys.argv[1]);workers=int(sys.argv[2]);start=time.monotonic()
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:result=list(pool.map(point,[(out,role,k,i) for role in ('source','control') for k in range(5) for i in range(1921)],chunksize=4))
 target=out/'interpolation-roundoff';target.mkdir(exist_ok=False)
 for role in ('source','control'):
  for k in range(5):
   records=[r for r in result if r['role']==role and r['reference']==k];assert len(records)==1921
   save(target/f'{role}-reference-{k}.json',dict(status='CONDITIONAL_ARITHMETIC_ENCLOSURE',records=records,scope='Accumulated retained reference prefix plus every actual observation interpolation stage; exact endpoints bypass interpolation.'))
 save(out/'interpolation-roundoff-completion.json',dict(status='CONDITIONAL_PASS',observations=len(result),workers=workers,elapsed_s=time.monotonic()-start,code={str(p):digest(p) for p in (pathlib.Path(__file__),pathlib.Path(__file__).parent/'roundoff.py')}))
 print('PASS arithmetic enclosures for all 19210 reference observations')
