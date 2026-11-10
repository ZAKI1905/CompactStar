"""Actual RK residual propagation with explicit trial dependency boundaries.

Each GSL attempt contributes n RHS records after its cached first derivative.
Only the accepted final attempt feeds its Apply endpoint; discarded event trials
and rejected attempts are retained and counted, never chained as evolution.
"""
import collections,time
from decimal import Decimal as D,localcontext
from residuals import old,jacobian,MODEL
from high_precision import Coupled,rows,lift

def audit(role,name,directory,initial_bound=None,method='rk8pd'):
 start=time.monotonic()
 with localcontext() as ctx:
  ctx.prec=110;m=Coupled(MODEL,role);J,box=jacobian(m);A,W=old.tableau(method);n=len(W)
  trace=rows(directory/'local-rhs.tsv');applies=collections.defaultdict(list)
  for r in rows(directory/'applies.tsv'):applies[int(r['context'])].append(r)
  bound=list(initial_bound or (D(0),)*3);records=[];rhsmax=[D(0)]*3;replaymax=[D(0)]*3;discarded=0;rejections=0
  def key(t,y):return (t,*y)
  def inside(y,u):
   assert box['x_min']<y[0]-u[0] and y[0]+u[0]<box['x_max'],'thermal stage enclosure outside piecewise domain'
   assert all(abs(y[j+1])+u[j+1]<box['eta_abs'][j] for j in range(2)),'chemical stage enclosure outside piecewise domain'
  for mapping in rows(directory/'maps.tsv'):
   group=applies[int(mapping['context'])];assert group,'missing Apply boundaries'
   first=group[0];t0=lift(first['left']);y0=tuple(lift(first[k]) for k in ('x0','e0','mu0'))
   known={key(t0,y0):bound.copy()};used=[]
   for call in group:
    rejections+=int(call['rejected'])
    a=lift(call['left']);b=lift(call['right']);h=b-a;assert h>0
    y0=tuple(lift(call[k]) for k in ('x0','e0','mu0'));y1=tuple(lift(call[k]) for k in ('x1','e1','mu1'))
    assert key(a,y0) in known,'unmapped event-trial initial state'
    incoming=known[key(a,y0)];begin=int(call['row_begin']);end=int(call['row_end']);segment=trace[begin:end]
    expected=n*(1+int(call['rejected']))
    assert len(segment) in (expected,expected+1),('unsupported apply trace',name,mapping['context'],len(segment),expected)
    assert all(int(r['context'])==int(mapping['context']) for r in segment)
    # GSL reuses the initial derivative on rejected retries. The last n records
    # are the accepted attempt's n-1 remaining stages and final derivative.
    tail=segment[-n:];stages=tail[:-1]
    # Cell-selection probes (cell=-1) do not update GSL's cached derivative.
    cached=next((r for r in reversed(trace[:end-n]) if int(r['cell'])>=0 and lift(r['t_s'])==a and tuple(lift(r[k]) for k in ('x','e','mu'))==y0),None)
    assert cached is not None,'cached first derivative not in trace';stages=[cached]+stages
    assert len(stages)==n
    assert tuple(lift(segment[-1][k]) for k in ('x','e','mu'))==y1
    fs=[[lift(r[k]) for k in ('dx','de','dmu')] for r in stages];defects=[]
    for j,r in enumerate(stages):
     yy=tuple(lift(r[k]) for k in ('x','e','mu'));inside(yy,(D(0),)*3)
     real=m.coupled_rhs(yy);defects.append([abs(fs[j][i]-real[i])+D('1e-90') for i in range(3)])
    errors=[]
    for j in range(n):
     yy=tuple(lift(stages[j][k]) for k in ('x','e','mu'))
     residual=[abs(yy[i]-y0[i]-h*sum(A[j][k]*fs[k][i] for k in range(j)))+D('1e-90') for i in range(3)]
     errors.append([incoming[i]+residual[i]+h*sum(abs(A[j][k])*(defects[k][i]+sum(J[i][v]*errors[k][v] for v in range(3))) for k in range(j)) for i in range(3)])
     inside(yy,errors[-1])
     for i in range(3):rhsmax[i]=max(rhsmax[i],defects[j][i])
    residual=[abs(y1[i]-y0[i]-h*sum(W[j]*fs[j][i] for j in range(n)))+D('1e-90') for i in range(3)]
    result=[incoming[i]+residual[i]+h*sum(abs(W[j])*(defects[j][i]+sum(J[i][v]*errors[j][v] for v in range(3))) for j in range(n)) for i in range(3)]
    inside(y1,result)
    previous=known.get(key(b,y1),(D(0),)*3)
    known[key(b,y1)]=[max(v,w) for v,w in zip(previous,result)];used.append(key(b,y1))
   replay=tuple(lift(mapping['replay_'+k]) for k in ('x','e','mu'));target=key(lift(mapping['right']),replay)
   assert target in known,'retained replay does not match an Apply endpoint'
   mismatch=[abs(lift(mapping[k])-replay[i]) for i,k in enumerate(('x','e','mu'))]
   bound=[v+w for v,w in zip(known[target],mismatch)];inside(tuple(lift(mapping[k]) for k in ('x','e','mu')),bound)
   for i in range(3):replaymax[i]=max(replaymax[i],mismatch[i])
   discarded+=sum(1 for r in group if key(lift(r['left']),tuple(lift(r[k]) for k in ('x0','e0','mu0')))==key(t0,y0))-1
   records.append(dict(time_s=mapping['right'],accumulated_arithmetic_enclosure=bound.copy(),retained_replay_difference=mismatch,apply_count=len(group)))
  return dict(status='CONDITIONAL_PIECEWISE_ARITHMETIC_ENCLOSURE',role=role,history=name,method=method,J=J,box=box,records=records,final_enclosure=bound,maximum_rhs_defect=rhsmax,maximum_replay_difference=replaymax,event_trial_resets=discarded,rejected_attempts=rejections,elapsed_s=time.monotonic()-start,assumptions=['Instrumented AuditStepper differs only by class name and passive Apply boundary hooks.','Actual selected Apply stages and retained/replay differences are enclosed against the continuous piecewise real RHS; discarded event trials do not accumulate as serial evolution.','Rational GSL tableaux, 110-digit arithmetic and existing conservative padding retained.','Excludes truncation and static normalization, which are bounded separately.'])
