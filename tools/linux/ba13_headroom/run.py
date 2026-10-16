#!/usr/bin/env python3
"""Bounded X2 assessment runner; raw results are write-once and source bound."""
import hashlib,json,os,pathlib,resource,subprocess,sys,time
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1];O=R/'ba13-headroom';TC=R/'toolchain'
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cpu():return sum(x.ru_utime+x.ru_stime for x in (resource.getrusage(resource.RUSAGE_SELF),resource.getrusage(resource.RUSAGE_CHILDREN)))
start=cpu();wall=time.time();records=[];phase=sys.argv[1]
def call(label,args):
 p=O/(label+'.json');err=O/(label+'.err');a=cpu();w=time.time()
 def cap():resource.setrlimit(resource.RLIMIT_CPU,(180,181))
 with p.open('x') as f,err.open('x') as e:
  q=subprocess.run(list(map(str,args)),stdout=f,stderr=e,timeout=200,preexec_fn=cap)
 records.append(dict(label=label,command=list(map(str,args)),returncode=q.returncode,cpu_seconds=cpu()-a,wall_seconds=time.time()-w,sha256=sha(p),stderr_sha256=sha(err)))
 if q.returncode:raise RuntimeError(label+' failed')
 return json.loads(p.read_text()) if not phase.startswith('build') else None
try:
 if phase.startswith('build'):
  for name,source in [('extended','x2_fixed_step.cpp'),('assessment','x2_uncertainty.cpp')]:
   call(phase+'-'+name,[TC/'bin/g++','-O3','-std=c++17','-ffp-contract=off','-Wall','-Wextra','-isystem',TC/'install/Release/Zaki/include',S/'tests/bnv'/source,'-o',O/(name+('-final' if phase=='build-final' else '-v2' if phase=='build-v2' else '')),TC/'install/Release/Zaki/lib/libZaki.a'])
 else:
  t=int(sys.argv[2]);old=R/'x2-fixed-step';target=json.loads((S/'docs/validation/x1-oracle-anchor/targets.json').read_text())[str(t)]
  hc=json.loads((old/f't{t}-n4096.json').read_text(),parse_float=str)['values'][8]
  def solve(tag,n=128,H=hc,cf='1e-4',eq=96):
   return call(f't{t}-{tag}',[O/'assessment',H,target['surfaceP_cgs'],n,cf,eq,0])
  if phase=='pilot-v2':
   # Null observer must reproduce all primary decimal output fields exactly.
   x1=json.loads((S/'docs/validation/x1-oracle-anchor/evidence/runs.json').read_text())
   guess=str(next(r['values']['Hc'] for r in x1 if r['target']==str(t) and r['variant']=='regularized' and r['tol']=='2e-13'))
   replay=call(f't{t}-replay-v2',[O/'extended',target['B'],guess,target['surfaceP_cgs'],128,'1e-4',96,0,0,0])
   assert replay['values']==json.loads((old/f't{t}-n128.json').read_text())['values']
   solve('n128');solve('n512',512)
  elif phase=='regression-final':
   x1=json.loads((S/'docs/validation/x1-oracle-anchor/evidence/runs.json').read_text())
   guess=str(next(r['values']['Hc'] for r in x1 if r['target']==str(t) and r['variant']=='regularized' and r['tol']=='2e-13'))
   call(f't{t}-regression-final',[O/'extended-final',target['B'],guess,target['surfaceP_cgs'],128,'1e-4',96,0,0,0])
   replay=json.loads((O/f't{t}-regression-final.json').read_text(),parse_float=str)
   original=json.loads((old/f't{t}-n128.json').read_text(),parse_float=str)
   assert replay==original, 'full decimal JSON regression'
  elif phase=='derivative-errors':
   from decimal import Decimal as D
   for variant,cf,eq in [('center','1e-2',96),('eos','1e-4',192)]:
    for sign in [-1,1]:solve(f'derivative-{variant}-{sign}',512,str(D(hc)*(1+sign*D('0.00005'))),cf,eq)
  elif phase=='defect512':
   call(f't{t}-defect512',[O/'assessment-v2',hc,target['surfaceP_cgs'],512,'1e-4',96,1])
  elif phase=='defect':
   call(f't{t}-defect',[O/'assessment-v2',hc,target['surfaceP_cgs'],128,'1e-4',96,1])
  elif phase=='sensitivity':
   from decimal import Decimal as D
   solve('center-small',512,cf='1e-6');solve('center-large',512,cf='1e-2');solve('eos192',512,eq=192)
   for frac in ['0.0001','0.00005']:
    for sign in [-1,1]:solve(f'h{frac}-{sign}',512,str(D(hc)*(1+sign*D(frac))))
  else:raise ValueError(phase)
finally:
 data=dict(phase=phase,target=sys.argv[2] if len(sys.argv)>2 else None,job=os.environ['SLURM_JOB_ID'],affinity=sorted(os.sched_getaffinity(0)),cpu_seconds=cpu()-start,wall_seconds=time.time()-wall,peak_rss_kib=max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),runs=records,source={str(p.relative_to(S)):sha(p) for p in [pathlib.Path(__file__),S/'tests/bnv/x2_fixed_step.cpp',S/'tests/bnv/x2_uncertainty.cpp']})
 with (O/f'receipt-{phase}-{sys.argv[2] if len(sys.argv)>2 else "build"}.json').open('x') as f:json.dump(data,f,indent=2)
 print(json.dumps(data))
