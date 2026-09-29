from pathlib import Path
import subprocess,os,json,time
q=Path(__file__).resolve().parent
src=Path('/Users/keeper/Documents/CompactStar/repo/CompactStar')
w=Path('/Users/keeper/Documents/CompactStar/worktrees')
old=w/'CompactStar-phase6a1-controlled-bnv-implementation';base=old/'build/phase6a1-controlled-bnv-debug';recon=w/'CompactStar-phase6a1-checkpoint-reconstruction-validation'
inputs=dict(oracle_root=recon/'build/phase6a1-checkpoint-reconstruction-validation/run/oracle',oracle_result=recon/'build/phase6a1-checkpoint-reconstruction-validation/oracle-result.json',matrix=recon/'docs/validation/phase6a1_checkpoint_reconstruction_solve_matrix.tsv',profile=base/'phase6a1-fresh-qualification/chemical-characterization/run-dhp6aivf/t8192-r80000',certificate=base/'phase6a1-fresh-qualification/certificate-t8192-r80000.txt',thermal=base/'phase6a1-frozen-certificate-final/chemical-profiles/target-0/thermal',frozen=base/'phase6a1-frozen-certificate-final/certificate.tsv',coefficients=base/'phase6a1-frozen-certificate-final/coefficients.tsv',entry=old/'docs/validation/phase6a1_controlled_bnv_entry_hashes.json')
run=q/'adr0017-old-Debug';run.mkdir()
env={**os.environ,'PATH':'/usr/bin:/bin:/opt/homebrew/bin:/opt/local/bin:'+os.environ['PATH'],'PYTHONDONTWRITEBYTECODE':'1'}
verify=src/'tests/bnv/adr0017_production_verify.py'
commands=[]
def call(cmd,name):
 commands.append([str(x) for x in cmd]);(run/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
 start=time.monotonic()
 with (run/(name+'.log')).open('w') as f:r=subprocess.run(cmd,cwd=src,env=env,stdout=f,stderr=f)
 print(name,r.returncode,time.monotonic()-start,flush=True)
 if r.returncode:raise SystemExit(r.returncode)
cmd=['/Users/keeper/miniforge3/bin/python3',verify,'authenticate']
for k,v in inputs.items():cmd+=['--'+k.replace('_','-'),v]
call(cmd+['--output',run/'auth'],'authentication')
call([q/'old-Debug/tests/adr0017_production_qualification',inputs['matrix'],inputs['profile'],inputs['certificate'],inputs['thermal'],run/'context',inputs['entry'],inputs['frozen'],inputs['coefficients'],run/'output',run/'auth/oracle-qualified.flag'],'execution')
call(['/Users/keeper/miniforge3/bin/python3',verify,'verify','--output',run/'output','--oracle-root',inputs['oracle_root'],'--oracle-result',inputs['oracle_result'],'--result',run/'result.json'],'verification')
