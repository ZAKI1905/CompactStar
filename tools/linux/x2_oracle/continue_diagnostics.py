#!/usr/bin/env python3
"""Resume unaffected predeclared diagnostics, preserving the failed binary64 run."""
import hashlib,json,os,pathlib,resource,subprocess,time,sys
SRC=pathlib.Path(__file__).resolve().parents[3];ROOT=SRC.parents[1];OUT=ROOT/'x2-fixed-step';DOC=SRC/'docs/validation/x2-fixed-step'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if len(sys.argv)>1 and sys.argv[1]=='build':
 src=SRC/'tests/bnv/x2_fixed_step.cpp';s=src.read_text();old='residual<32*std::numeric_limits<Real>::epsilon()';assert old in s
 s=s.replace(old,'residual<std::min(32*std::numeric_limits<Real>::epsilon(),Real(5e-16))')
 path=OUT/'build/double-newton.cpp';path.write_text(s)
 tc=ROOT/'toolchain';subprocess.run([str(tc/'bin/g++'),'-O3','-std=c++17','-ffp-contract=off','-DX2_DOUBLE','-isystem',str(tc/'install/Release/Zaki/include'),str(path),'-o',str(OUT/'build/double-newton'),str(tc/'install/Release/Zaki/lib/libZaki.a')],check=True)
 fs=[path,OUT/'build/double-newton',pathlib.Path(__file__),src,OUT/'build/extended']
 (DOC/'continuation-gate.json').write_text(json.dumps({'files':{str(f):sha(f) for f in fs},'purpose':'Correct early Newton stop, not stage gate, in binary64 diagnostic only; preserve first failure.'},indent=2)+'\n');sys.exit()
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});assert os.getenv('SLURM_JOB_ID')
for p,h in json.loads((DOC/'continuation-gate.json').read_text())['files'].items():assert sha(pathlib.Path(p))==h,p
def cpu():return sum(r.ru_utime+r.ru_stime for r in [resource.getrusage(resource.RUSAGE_SELF),resource.getrusage(resource.RUSAGE_CHILDREN)])
prior=sum(json.loads(p.read_text())['cpu_seconds'] for p in OUT.glob('receipt-*.json'));assert prior<650
q=json.loads((SRC/'docs/validation/x1-oracle-anchor/targets.json').read_text())['0'];x1=json.loads((SRC/'docs/validation/x1-oracle-anchor/evidence/runs.json').read_text())
guess=str(next(r['values']['Hc'] for r in x1 if r['target']=='0' and r['variant']=='regularized' and r['tol']=='2e-13'))
tasks=[]
for n in [512,1024,2048]:
 hc=json.loads((OUT/f't0-n{n}.json').read_text(),parse_float=str)['values'][8]
 for tag,cf,eq,raw,plain,mode,fixed,h in [('center-1e-2','1e-2',96,0,0,'extended',0,guess),('center-1e-6','1e-6',96,0,0,'extended',0,guess),('eos192','1e-4',192,0,0,'extended',0,guess),('plain','1e-4',96,0,1,'extended',0,guess),('raw','1e-4',96,1,0,'extended',0,guess),('double-newton','1e-4',96,0,0,'double-newton',1,hc)]:
  label=f'{tag}-n{n}'
  if (OUT/(label+'.json')).exists():continue
  tasks.append((label,[str(OUT/'build'/mode),q['B'],h,q['surfaceP_cgs'],str(n),cf,str(eq),str(raw),str(plain),str(fixed)]))
for v in ['regularized','unregularized']:
 hc=str(next(r['values']['Hc'] for r in x1 if r['target']=='0' and r['variant']==v and r['tol']=='2e-13'))
 tasks.append((f'fixed-x1-{v}',[str(OUT/'build/extended'),q['B'],hc,q['surfaceP_cgs'],'2048','1e-4','96','0','0','1']))
start=cpu();wall=time.time();runs=[]
try:
 for label,cmd in tasks:
  assert prior+cpu()-start<750
  before=cpu();p=subprocess.run(cmd,capture_output=True,text=True,timeout=100)
  out=OUT/(label+'.json');err=OUT/(label+'.err');assert not out.exists();out.write_text(p.stdout);err.write_text(p.stderr)
  runs.append(dict(label=label,command=cmd,returncode=p.returncode,cpu_seconds=cpu()-before,stdout_sha256=sha(out),stderr_sha256=sha(err)))
  # A binary64 diagnostic failure does not discard other independent ablations.
  if p.returncode and not label.startswith('double'):raise RuntimeError(label+': '+p.stderr)
finally:
 rec=dict(phase='continuation',job=os.getenv('SLURM_JOB_ID'),affinity=sorted(os.sched_getaffinity(0)),cpu_seconds=cpu()-start,prior_cpu_seconds=prior,wall_seconds=time.time()-wall,runs=runs,gate_sha256=sha(DOC/'continuation-gate.json'))
 (OUT/'receipt-continuation.json').write_text(json.dumps(rec,indent=2)+'\n')
