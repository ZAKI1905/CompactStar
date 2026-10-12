#!/usr/bin/env python3
"""One isolated predeclared common-B/thermal study case; no trajectory."""
import sys,os,json,time,resource,hashlib,subprocess
from pathlib import Path
sys.dont_write_bytecode=True
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')));RAW=S.parents[1]/'ba13-recovery/run1/raw'
sys.path.insert(0,str(S/'tests/bnv'));from produce_frozen_certificates import thermal_source
case=sys.argv[1];index=int(sys.argv[2])
cases={'nominal':(80000,8192,1,5e-11,1e-14),
 'r20000':(20000,8192,1,5e-11,1e-14),'r40000':(40000,8192,1,5e-11,1e-14),'r160000':(160000,8192,1,5e-11,1e-14),'r320000':(320000,8192,1,5e-11,1e-14),
 'e2048':(80000,2048,1,5e-11,1e-14),'e4096':(80000,4096,1,5e-11,1e-14),'e16384':(80000,16384,1,5e-11,1e-14),'e32768':(80000,32768,1,5e-11,1e-14),
 'center2':(80000,8192,2,5e-11,1e-14),'center05':(80000,8192,.5,5e-11,1e-14),'center025':(80000,8192,.25,5e-11,1e-14),
 'b5e-12':(80000,8192,1,5e-12,1e-14),'b5e-13':(80000,8192,1,5e-13,1e-14),
 'floor15':(80000,8192,1,5e-11,1e-15),'floor16':(80000,8192,1,5e-11,1e-16)}
radial,eos,center,tau,lower=cases[case]
dest=R/'runs'/f'{case}-{index}';partial=dest.with_suffix('.partial');dest.parent.mkdir(exist_ok=True,parents=True)
if dest.exists() or partial.exists():raise RuntimeError('fresh study case required')
exe=R/'build/tests/ba13_uncertainty_pilot';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
args=[str(RAW),str(partial),str(radial),str(eos),str(center),str(tau),str(index),str(lower)]
start=time.monotonic();receipt=dict(case=case,index=index,parameters=cases[case],job=os.getenv('SLURM_JOB_ID'),node=os.uname().nodename,
 executable_sha256=sha(exe),predeclaration_sha256=sha(S/'docs/validation/ba13-uncertainty/PREDECLARATION.md'),predeclaration_commit='b48c642',commands=[])
try:
 with dest.with_suffix('.log').open('x') as log:
  for mode in ['star','thermal']:
   cmd=[str(exe),mode]+args;receipt['commands'].append(cmd)
   subprocess.run(cmd,cwd=S,stdout=log,stderr=subprocess.STDOUT,check=True)
   if mode=='star':thermal_source(partial/'profile')
 receipt['status']='COMPLETE_MEASUREMENTS_NOT_BOUNDS'
except Exception as e:
 receipt['status']='REFUSED';receipt['error']=str(e)
finally:
 ru=resource.getrusage(resource.RUSAGE_CHILDREN);own=resource.getrusage(resource.RUSAGE_SELF)
 receipt.update(wall_seconds=time.monotonic()-start,cpu_seconds=ru.ru_utime+ru.ru_stime+own.ru_utime+own.ru_stime,max_rss_kib=max(ru.ru_maxrss,own.ru_maxrss))
 dest.with_suffix('.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
if receipt['status']=='REFUSED':raise SystemExit(1)
files={str(p.relative_to(partial)):sha(p) for p in partial.rglob('*') if p.is_file()}
(partial/'COMPLETE.json').write_text(json.dumps(dict(files=files,stage='uncertainty measurement'),indent=2)+'\n')
for p in partial.rglob('*'):
 if p.is_file():
  with p.open('rb') as f:os.fsync(f.fileno())
partial.rename(dest)
print(receipt)
