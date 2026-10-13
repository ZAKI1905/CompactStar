"""Sequential Slurm numerical budget controller, fixed cases, no adaptive tuning."""
import csv,json,os,resource,subprocess,sys,time,signal
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=S.parents[1];OUT=R/'ba13-onset-pilot';PY=sys.executable;HERE=Path(__file__).resolve().parent
assert os.getenv('SLURM_JOB_ID'), 'numerical work requires Slurm'
assert not subprocess.check_output(['git','status','--porcelain'],cwd=S),'commit tools/predictions before outcomes'
# Conservative failed-job/startup charge; discovery counted exactly from measured process.
prior=json.loads((OUT/'discovery-resource.json').read_text())['cpu_seconds']+2
records=[];start=time.monotonic();limit=3500-prior

def execute(name,cmd,ceiling=300):
 remaining=limit-(time.monotonic()-start)
 if remaining<ceiling:
  records.append(dict(name=name,status='RESOURCE_STOP',required_wall=ceiling,remaining_wall=remaining));save();raise SystemExit('resource stop; no additional numerical work')
 ru=resource.getrusage(resource.RUSAGE_CHILDREN);before=ru.ru_utime+ru.ru_stime;t=time.monotonic()
 try:
  with (OUT/(name+'.log')).open('x') as log:
   p=subprocess.Popen(list(map(str,cmd)),stdout=log,stderr=subprocess.STDOUT,cwd=S,start_new_session=True)
   try:code=p.wait(timeout=min(600,remaining))
   except subprocess.TimeoutExpired:
    os.killpg(p.pid,signal.SIGKILL);p.wait();code=124
 except OSError as error:
  code=125
 ru=resource.getrusage(resource.RUSAGE_CHILDREN);records.append(dict(name=name,command=list(map(str,cmd)),exit_code=code,cpu_seconds=ru.ru_utime+ru.ru_stime-before,wall_seconds=time.monotonic()-t,max_rss_kib=ru.ru_maxrss,job=os.getenv('SLURM_JOB_ID')));save();print(name,code,records[-1]['cpu_seconds'],flush=True);return code==0

def save():
 u=resource.getrusage(resource.RUSAGE_SELF);v=resource.getrusage(resource.RUSAGE_CHILDREN)
 (OUT/'resources.json').write_text(json.dumps(dict(prior_cpu_charged=prior,controller_cpu=u.ru_utime+u.ru_stime,children_cpu=v.ru_utime+v.ru_stime,total_cpu=prior+u.ru_utime+u.ru_stime+v.ru_utime+v.ru_stime,wall_seconds=time.monotonic()-start,ceiling=3600,reserve=100,records=records),indent=2)+'\n')

admissible={}
for n,cap in [(2048,11),(4096,11),(8192,11),(16384,11),(32768,11),(65536,11),(32768,12)]:
 ok=execute(f'table-{n}-{cap}',[OUT/'build/onset_table_audit',OUT/'tables'/f'{n}-{cap}',n,cap],60);admissible[n,cap]=ok
 if not ok and n<=8192:raise SystemExit('STOP P5 identity/admissibility failed')
execute('high-precision',[PY,'-B',HERE/'high_precision.py',OUT],120)
# Refuse new stars if independent binary64 ordering cannot represent their table.
hp=json.loads((OUT/'high-precision.json').read_text())
for row in hp['pairs']:
 if row['case']!='historical32768' and not (row['high_precision_order'] and row['rounded_order']):
  n,cap=map(int,row['case'].split('-'));admissible[n,cap]=False
cases=[('r640000',640000,8192,11,i) for i in [0,10,20]]
cases += [(f'e{n}c',80000,n,11,i) for n in [16384,32768,65536] for i in [0,10,20]]
cases += [('e32768d12',80000,32768,12,i) for i in [0,10,20]]
cases += [(f'r{n}',n,8192,11,i) for n in range(80001,80006) for i in [0,10,20]]
cases += [('cross',640000,16384,11,20)]
for name,rad,n,cap,i in cases:
 if not admissible[n,cap]:records.append(dict(name=f'{name}-{i}',status='EOS_REFUSED'));save();continue
 if execute(f'run-{name}-{i}',[PY,'-B',HERE/'run_case.py',name,i,rad,n,cap],200 if rad==640000 else 90):
  execute(f'analyze-{name}-{i}',[PY,'-B',HERE/'analyze_case.py',name,i],250 if name=='r640000' else 60)
for i in [5,15]:
 for c in ['r20000','r40000','nominal','r160000','r320000']:
  execute(f'secondary-{c}-{i}',[PY,'-B',HERE/'analyze_case.py',c,i,'old'],60)
# Independent stars use original fixed B, own analytic EOS, and same nominal surface pressure.
for i in [0,10,20]:
 path=R/'ba13-uncertainty/runs'/f'nominal-{i}'
 row=next(csv.DictReader((path/'star.tsv').open(),delimiter='\t'))
 # Read final profile line, pressure conversion only (no integration).
 last=(path/'profile/profile.tsv').read_text().splitlines()[-1].split();ps=float(last[4])/1.3235754613110692e-6*1.602176634e33
 for tol in ['2e-10','2e-12']:
  execute(f'oracle-{i}-{tol}',[OUT/'build/onset_structure_oracle',row['B_target'],row['rho_c'],ps,tol],120)
save()
