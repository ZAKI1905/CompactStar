#!/usr/bin/env python3
"""Independent piecewise coupled accuracy at every native reconstructed sample.

Each interval starts at its exact retained native left state; this never reruns
or replaces a production trajectory. No Hermite-temperature-path assumption.
"""
import bisect,concurrent.futures,json,pathlib,sys,time
from decimal import Decimal as D,localcontext
from piecewise import MODEL,Coupled,domain,SCALE,step
from run import save,digest
from high_precision import rows,lift
from verify_chemical_reference import corner_uncertainty

def advance(m,box,t,y,u,target):
 L=box['L'];M=box['M']
 while t<target:
  cell=bisect.bisect_right(m.knots,y[0])-1;c=box['cells'][cell]
  distance=min(y[0]-m.knots[cell],m.knots[cell+1]-y[0]);safe=(distance-u)/(2*M)
  bridge=safe<D('.01');h=min(D('1000000'),target-t,D('.01') if bridge else safe)
  assert h>0 and t+h>t
  reach=u+M*h
  assert box['lo']<y[0]-reach and y[0]+reach<box['hi']
  assert all(abs(y[j+1])+reach*SCALE[j+1]<box['Y'][j] for j in range(2))
  q=h/c['r'];assert q<1
  local=L*M*h*h*(2*L*h).exp() if bridge else c['R']*q**5/(1-q)
  y=step(m.coupled_rhs,y,h);u=(L*h).exp()*u+local+D('1e-80');t+=h
  assert u<D('1e-10')
 return t,y,u

def assess_one(task):
 out,role,directory=task;start=time.monotonic()
 with localcontext() as ctx:
  ctx.prec=110;m=Coupled(MODEL,role);box=domain(m)
  pre=json.loads((out/'preflight.json').read_text());forcing=D(pre['roles'][role]['normalization_rhs_bound'])
  histories=rows(directory/'main.accepted_states.tsv');accepted={lift(r['t_after_s']):r for r in histories}
  data=rows(directory/'diagnostics-1920.tsv')
  if len(data)!=1921:return dict(role=role,configuration=directory.name,status='INCOMPLETE')
  t=D(-1);left=None;y=None;u=None;results=[]
  for index,d in enumerate(data):
   target=lift(d['t_s']);p=directory/'checkpoints-1920'/f'point-{index}'
   cp=json.loads((p/'checkpoint.json').read_text())
   if cp['method']=='main-endpoint':
    results.append(dict(index=index,status='EXACT_ENDPOINT'));continue
   a=lift(cp['main_left_time_s']);initial=tuple(lift(v) for v in cp['main_left_state'])
   if a!=left:
    if a:assert initial==tuple(lift(accepted[a][k]) for k in ('x_state','eta_e_MeV','eta_mu_MeV'))
    else:assert initial==(D(0),)*3
    left=a;t=a;y=initial;u=D(0)
   assert target>t
   t,y,u=advance(m,box,t,y,u,target)
   primary=u+forcing*(target-left)*(box['L']*(target-left)).exp()
   sr={r['channel']:r for r in rows(p/'state.tsv')};dr={r['channel']:r for r in rows(p/'diagnostics.tsv')}
   errors=[D(0) if role=='control' and j>0 else primary*SCALE[j] for j in range(3)]
   ratios=[(abs(lift(cp['O2_state'][j])-y[j])+errors[j])/lift(sr[key]['propagated_state_radius']) if lift(sr[key]['propagated_state_radius']) else D(0) for j,key in enumerate(('x','eta_e','eta_mu'))]
   for j,key in enumerate(('x','eta_e','eta_mu')):
    if not lift(sr[key]['propagated_state_radius']):assert lift(cp['O2_state'][j])==y[j] and errors[j]==0
   er=m.energy(y[1:]);ue=corner_uncertainty(m.energy,y[1:],tuple(errors[j] for j in (1,2)))[0]
   numerator=abs(lift(dr['Echem']['O2'])-er)+ue;denominator=lift(dr['Echem']['U'])
   eratio=numerator/denominator if denominator else D(0)
   if not denominator:assert numerator==0
   results.append(dict(index=index,status='CONDITIONAL_PASS' if max(ratios+[eratio])<=1 else 'FAIL',state_ratios=ratios,Echem_ratio=eratio,scaled_reference_error=primary))
  vals=[r for r in results if 'state_ratios' in r]
  result=dict(role=role,configuration=directory.name,status='FAIL' if any(r['status']=='FAIL' for r in results) else 'CONDITIONAL_PASS',points=len(results),worst_state=max(vals,key=lambda r:max(r['state_ratios'])),worst_Echem=max(vals,key=lambda r:r['Echem_ratio']),elapsed_s=time.monotonic()-start,results=results,scope='Coupled ODE endpoint enclosure from each exact native left state, including primary-formula normalization. No assumption about a Hermite thermal path; not a native main-prefix or arithmetic-only audit.')
  target=out/'local-assessment';target.mkdir(exist_ok=True);save(target/f'{role}-{directory.name}.json',result)
  return {k:result[k] for k in ('role','configuration','status','points','worst_state','worst_Echem','elapsed_s')}

if __name__=='__main__':
 out=pathlib.Path(sys.argv[1]);workers=int(sys.argv[2]);tasks=[]
 for role in ('source','control'):
  for p in sorted((out/'native'/role).glob('tier-*')):
   if (p/'diagnostics-1920.tsv').exists() and len(rows(p/'diagnostics-1920.tsv'))==1921:tasks.append((out,role,p))
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:results=list(pool.map(assess_one,tasks))
 save(out/'local-assessment.json',dict(status='CANDIDATE',workers=workers,configurations=results))
