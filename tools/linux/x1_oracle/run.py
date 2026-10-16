"""Sequential Slurm-only X1 runner, measured CPU budget and authenticated sources."""
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

SRC=Path(__file__).resolve().parents[3]
DOC=SRC/'docs/validation/x1-oracle-anchor'
ROOT=SRC.parents[1]
OUT=ROOT/'x1-oracle-anchor'
start=time.monotonic()
if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Slurm required')
if int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=1:raise SystemExit('one CPU required')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
gate=json.loads((DOC/'execution-gate.json').read_text())
for path,h in gate['files'].items():
    if sha(Path(path))!=h:raise SystemExit('authentication failed: '+path)
for path,h in gate['review_files'].items():
    if sha(Path(path))!=h:raise SystemExit('review authentication failed: '+path)
if (OUT/'runs.json').exists():raise SystemExit('refuse overwrite or repeat numerical run')
head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=SRC,text=True).strip()
def cpu():
    s=resource.getrusage(resource.RUSAGE_SELF);c=resource.getrusage(resource.RUSAGE_CHILDREN)
    return s.ru_utime+s.ru_stime+c.ru_utime+c.ru_stime
receipts=[];rows=[]
def execute(name,cmd):
    remaining=120-cpu()
    if remaining<62:raise RuntimeError('budget reserve prevents next numerical process')
    before=cpu();wall=time.monotonic()
    def limits():
        resource.setrlimit(resource.RLIMIT_CPU,(60,60))
        resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    receipt=OUT/(name+'.time')
    with (OUT/(name+'.out')).open('w') as out,(OUT/(name+'.err')).open('w') as err:
        p=subprocess.run(['/usr/bin/time','-v','-o',str(receipt)]+cmd,stdout=out,stderr=err,preexec_fn=limits,timeout=61)
    used=cpu()-before
    record=dict(name=name,command=cmd,exit_code=p.returncode,cpu_seconds=used,wall_seconds=time.monotonic()-wall,slurm_job=os.environ['SLURM_JOB_ID'],head=head)
    receipts.append(record)
    (OUT/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
    if p.returncode or used>60 or cpu()>120:raise RuntimeError('solve/check/resource failure: '+name)
    return (OUT/(name+'.out')).read_text()
status='INCONCLUSIVE'
try:
    checks=execute('mathematical-checks',[str(OUT/'build/x1_regularized_oracle'),'--check'])
    if not json.loads(checks)['pass']:raise RuntimeError('mathematical checks failed')
    targets=json.loads((DOC/'targets.json').read_text())
    for variant,tols,exe in [('unregularized',['2e-11','2e-13'],'onset_structure_oracle'),('regularized',['2e-10','2e-11','2e-12','2e-13'],'x1_regularized_oracle')]:
        for i in ['0','10','20']:
            for tol in tols:
                target=targets[i];name=f'{variant}-{i}-{tol}'
                cmd=[str(OUT/'build'/exe),target['B'],target['rho_guess'],target['surfaceP_cgs'],tol]
                text=execute(name,cmd);values=list(map(float,text.split()))
                if len(values)!=9 or not all(map(math.isfinite,values)):raise RuntimeError('invalid output')
                v=dict(zip(['B','R','M','r_e','L','C','rb','mu_inf','Hc'],values))
                row=dict(variant=variant,target=i,tol=tol,values=v)
                if variant=='regularized':row['checks']=json.loads((OUT/(name+'.err')).read_text())
                rows.append(row);(OUT/'runs.json').write_text(json.dumps(rows,indent=2)+'\n')
                if abs(v['B']/float(target['B'])-1)>5e-13:raise RuntimeError('target alignment failed: '+name)
    execute('score',[sys.executable,'-I','-B',str(SRC/'tools/linux/x1_oracle/score.py'),str(OUT)])
    status='COMPLETE_NOT_QUALIFIED'
except Exception as e:
    (OUT/'STOP.txt').write_text(str(e)+'\n');print(str(e),file=sys.stderr)
finally:
    s=resource.getrusage(resource.RUSAGE_SELF);c=resource.getrusage(resource.RUSAGE_CHILDREN)
    ledger=dict(status=status,total_cpu_seconds=cpu(),controller_cpu=s.ru_utime+s.ru_stime,children_cpu=c.ru_utime+c.ru_stime,wall_seconds=time.monotonic()-start,max_child_rss_kib=c.ru_maxrss,max_self_rss_kib=s.ru_maxrss,job=os.environ['SLURM_JOB_ID'],cpus=os.environ['SLURM_CPUS_PER_TASK'],node=os.environ.get('SLURMD_NODENAME'),head=head,gate_sha256=sha(DOC/'execution-gate.json'),runs=len(rows))
    (OUT/'resources.json').write_text(json.dumps(ledger,indent=2)+'\n');print(json.dumps(ledger))
if status!='COMPLETE_NOT_QUALIFIED':sys.exit(1)
