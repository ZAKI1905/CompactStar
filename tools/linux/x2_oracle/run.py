#!/usr/bin/env python3
"""Source-bound, single-CPU Slurm runner. Never overwrites a numerical result."""
import hashlib,json,math,os,pathlib,resource,subprocess,sys,time
SRC=pathlib.Path(__file__).resolve().parents[3]
ROOT=SRC.parents[1]; OUT=ROOT/'x2-fixed-step'; DOC=SRC/'docs/validation/x2-fixed-step'
phase=sys.argv[1]
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
assert len(os.sched_getaffinity(0))==1 and os.environ.get('SLURM_JOB_ID')
def cpu():
 return sum(r.ru_utime+r.ru_stime for r in [resource.getrusage(resource.RUSAGE_SELF),resource.getrusage(resource.RUSAGE_CHILDREN)])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
gate=json.loads((DOC/'execution-gate.json').read_text())
for p,h in gate['files'].items():assert sha(pathlib.Path(p))==h,p
prior=sum(json.loads(p.read_text())['cpu_seconds'] for p in OUT.glob('receipt-*.json'))
start=cpu(); wall=time.time();rec=[]
def limits():resource.setrlimit(resource.RLIMIT_CPU,(100,101))
def call(label,args):
 if prior+cpu()-start>750:raise RuntimeError('CPU reserve stop')
 out=OUT/(label+'.json');err=OUT/(label+'.err')
 if out.exists():raise RuntimeError('result already exists: '+label)
 before=cpu();t=time.time()
 p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=110,preexec_fn=limits)
 out.write_text(p.stdout);err.write_text(p.stderr)
 rec.append(dict(label=label,command=list(map(str,args)),returncode=p.returncode,cpu_seconds=cpu()-before,wall_seconds=time.time()-t,stdout_sha256=sha(out),stderr_sha256=sha(err)))
 if p.returncode:raise RuntimeError(label+': '+p.stderr)
 return json.loads(p.stdout)
targets=json.loads((SRC/'docs/validation/x1-oracle-anchor/targets.json').read_text())
x1=json.loads((SRC/'docs/validation/x1-oracle-anchor/evidence/runs.json').read_text())
def solve(t,n,label=None,center='1e-4',eq=96,raw=0,plain=0,mode='extended',fixed=0,hc=None):
 if hc is None:hc=str(next(r['values']['Hc'] for r in x1 if r['target']==t and r['variant']=='regularized' and r['tol']=='2e-13'))
 q=targets[t]
 return call(label or f't{t}-n{n}',[str(OUT/'build'/mode),q['B'],hc,q['surfaceP_cgs'],str(n),center,str(eq),str(raw),str(plain),str(fixed)])
try:
 if phase=='pilot':
  call('checks',[str(OUT/'build/extended'),'--check'])
  for n in [32,64,128,256]:solve('0',n)
 elif phase=='ladder':
  for t in ['0','10','20']:
   for n in [32,64,128,256,512,1024]:
    if t=='0' and n<=256:continue
    solve(t,n)
 elif phase in ['2048','4096']:
  for t in ['0','10','20']:solve(t,int(phase))
 elif phase=='ablations':
  ns=[512,1024]+([2048] if (OUT/'t0-n2048.json').exists() else [])
  for n in ns:
   for cf in ['1e-2','1e-6']:solve('0',n,f'center-{cf}-n{n}',center=cf)
   solve('0',n,f'eos192-n{n}',eq=192)
   solve('0',n,f'plain-n{n}',plain=1)
   solve('0',n,f'raw-n{n}',raw=1)
   # Fixed H binary64 avoids conflating its root-resolution floor with arithmetic.
   hc=json.loads((OUT/f't0-n{n}.json').read_text(),parse_float=str)['values'][8]
   solve('0',n,f'double-n{n}',mode='double',fixed=1,hc=hc)
  for v in ['regularized','unregularized']:
   hc=str(next(r['values']['Hc'] for r in x1 if r['target']=='0' and r['variant']==v and r['tol']=='2e-13'))
   solve('0',max(ns),f'fixed-x1-{v}',fixed=1,hc=hc)
 else:raise ValueError(phase)
finally:
 receipt=dict(phase=phase,job=os.environ['SLURM_JOB_ID'],affinity=sorted(os.sched_getaffinity(0)),cpu_seconds=cpu()-start,prior_cpu_seconds=prior,wall_seconds=time.time()-wall,peak_rss_kib=max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),execution_gate_sha256=sha(DOC/'execution-gate.json'),runs=rec)
 (OUT/f'receipt-{phase}.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt))
