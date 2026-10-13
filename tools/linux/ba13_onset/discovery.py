"""Reconstruction implementation cross-check and EOS discovery prediction freeze."""
import json,resource,time
from pathlib import Path
from equilibrium import analyze
from high_precision import constants
from freeze import predict
S=Path(__file__).resolve().parents[3];R=S.parents[1];out=R/'ba13-onset-pilot';old=R/'ba13-uncertainty/runs';start=time.monotonic()
hp=constants();(out/'constants.json').write_text(json.dumps(hp,indent=2)+'\n');d={}
for case in ['r20000','r40000','nominal','r160000','r320000','e2048','e4096']:
 for i in [0,10,20]:
  d[f'{i}:{case}']=analyze(old/f'{case}-{i}',hp);print(case,i,flush=True)
(out/'discovery.json').write_text(json.dumps(d,indent=2)+'\n')
p={}
for i in [0,10,20]:
 for key in ['L_nu','C_nu','r_nu']:
  p[f'{i}:{key}']=predict([d[f'{i}:{c}'][key] for c in ['e2048','e4096','nominal']])
(S/'docs/validation/ba13-onset-pilot/evidence/eos-predictions.json').write_text(json.dumps(p,indent=2)+'\n')
u=resource.getrusage(resource.RUSAGE_SELF)
(out/'discovery-resource.json').write_text(json.dumps(dict(cpu_seconds=u.ru_utime+u.ru_stime,wall_seconds=time.monotonic()-start,max_rss_kib=u.ru_maxrss,job=__import__('os').getenv('SLURM_JOB_ID')),indent=2)+'\n')
