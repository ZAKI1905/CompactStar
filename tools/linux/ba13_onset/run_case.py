"""One new immutable target-B case through unchanged scientific owners."""
import sys,os,json,time,resource,subprocess,hashlib
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=S.parents[1];OUT=R/'ba13-onset-pilot';RAW=R/'ba13-recovery/run1/raw'
sys.path.insert(0,str(S/'tests/bnv'));from produce_frozen_certificates import thermal_source
name=sys.argv[1];i=int(sys.argv[2]);rad=int(sys.argv[3]);eo=int(sys.argv[4]);cap=int(sys.argv[5]);dest=OUT/'runs'/f'{name}-{i}';partial=dest.with_suffix('.partial');dest.parent.mkdir(exist_ok=True)
assert not dest.exists() and not partial.exists()
exe=OUT/'build/onset_pilot';args=[str(RAW),str(partial),str(rad),str(eo),'1','5e-11',str(i),'1e-14',str(cap)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rec=dict(case=name,target=i,radial=rad,eos=eo,cap=cap,job=os.getenv('SLURM_JOB_ID'),executable_sha256=sha(exe),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=S,text=True).strip(),predeclaration=sha(S/'docs/validation/ba13-onset-pilot/PREDECLARATION.md'),predictions=sha(S/'docs/validation/ba13-onset-pilot/evidence/predictions.json'),eos_predictions=sha(S/'docs/validation/ba13-onset-pilot/evidence/eos-predictions.json'),commands=[]);start=time.monotonic()
try:
 for mode in ['star','thermal']:
  cmd=[str(exe),mode]+args;rec['commands'].append(cmd);subprocess.run(cmd,cwd=S,check=True)
  if mode=='star':thermal_source(partial/'profile')
 files={str(p.relative_to(partial)):sha(p) for p in partial.rglob('*') if p.is_file()}
 (partial/'COMPLETE.json').write_text(json.dumps(dict(files=files),indent=2)+'\n');partial.rename(dest);rec['status']='COMPLETE_NOT_QUALIFIED'
except Exception as e:rec['status']='REFUSED';rec['error']=str(e)
finally:
 u=resource.getrusage(resource.RUSAGE_SELF);v=resource.getrusage(resource.RUSAGE_CHILDREN);rec.update(cpu_seconds=u.ru_utime+u.ru_stime+v.ru_utime+v.ru_stime,wall_seconds=time.monotonic()-start,max_rss_kib=max(u.ru_maxrss,v.ru_maxrss))
 dest.with_suffix('.receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
if rec['status']=='REFUSED':sys.exit(1)
