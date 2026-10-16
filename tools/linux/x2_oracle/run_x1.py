#!/usr/bin/env python3
"""Only preobserved X1 tolerances; diagnostic replays, no rescoring."""
import hashlib,json,os,pathlib,resource,subprocess,time
SRC=pathlib.Path(__file__).resolve().parents[3];OUT=SRC.parents[1]/'x2-fixed-step';DOC=SRC/'docs/validation/x2-fixed-step'
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});assert os.getenv('SLURM_JOB_ID')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cpu():return sum(r.ru_utime+r.ru_stime for r in [resource.getrusage(resource.RUSAGE_SELF),resource.getrusage(resource.RUSAGE_CHILDREN)])
gate=json.loads((DOC/'x1-diagnostic-gate.json').read_text())
for p,h in gate['files'].items():assert sha(pathlib.Path(p))==h,p
prior=sum(json.loads(p.read_text())['cpu_seconds'] for p in OUT.glob('receipt-*.json'));assert prior<750
q=json.loads((SRC/'docs/validation/x1-oracle-anchor/targets.json').read_text())['0']
x1=json.loads((SRC/'docs/validation/x1-oracle-anchor/evidence/runs.json').read_text())
hc=str(next(r['values']['Hc'] for r in x1 if r['target']=='0' and r['variant']=='regularized' and r['tol']=='2e-13'))
start=cpu();wall=time.time();runs=[]
try:
 for variant in ['regularized','unregularized']:
  cases=[]
  for tol in ['2e-10','2e-11','2e-12','2e-13']:
   cases += [(f'baseline-{tol}',tol,{}),(f'fixed-{tol}',tol,{'X2_FIXED_H':'1'})]
  for tag,env in [('center-small',{'X2_CENTER':'1e-12'}),('center-large',{'X2_CENTER':'1e-8'}),('center-correct',{'X2_CORRECT_CENTER':'1'})]:cases.append((tag,'2e-13',env))
  if variant=='regularized':cases.append(('plain','2e-13',{'X2_PLAIN':'1'}))
  for tag,tol,env in cases:
   name=f'x1-{variant}-{tag}';out=OUT/(name+'.out');err=OUT/(name+'.err');assert not out.exists()
   cmd=[str(OUT/'build'/f'x1-{variant}'),q['B'],hc if 'X2_FIXED_H' in env else q['rho_guess'],q['surfaceP_cgs'],tol]
   t=cpu();p=subprocess.run(cmd,env=dict(os.environ,**env),capture_output=True,text=True,timeout=30);out.write_text(p.stdout);err.write_text(p.stderr)
   runs.append(dict(name=name,command=cmd,env=env,cpu_seconds=cpu()-t,returncode=p.returncode,stdout_sha256=sha(out),stderr_sha256=sha(err)))
   if p.returncode:raise RuntimeError(name+': '+p.stderr)
finally:
 rec=dict(phase='x1-diagnostics',job=os.getenv('SLURM_JOB_ID'),affinity=sorted(os.sched_getaffinity(0)),cpu_seconds=cpu()-start,prior_cpu_seconds=prior,wall_seconds=time.time()-wall,runs=runs,gate_sha256=sha(DOC/'x1-diagnostic-gate.json'))
 (OUT/'receipt-x1-diagnostics.json').write_text(json.dumps(rec,indent=2)+'\n')
