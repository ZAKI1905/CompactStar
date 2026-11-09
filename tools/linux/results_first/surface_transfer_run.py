#!/usr/bin/env python3
"""Bounded Slurm wrapper; derived from the existing static assessment runner."""
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time


def cpu():
    r=resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+r.ru_utime+r.ru_stime


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
    source=Path(os.environ['RF_SOURCE']).resolve();output=Path(os.environ['RF_OUTPUT']).resolve()
    commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    status=subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True)
    (output/'source-preflight.json').write_text(json.dumps(dict(actual_commit=commit,requested_commit=os.environ['RF_COMMIT'],status=status),indent=2)+'\n')
    if commit!=os.environ['RF_COMMIT'] or status:raise RuntimeError('clean committed source required')
    if source==output or source in output.parents:raise RuntimeError('external output required')
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    start,wall=cpu(),time.monotonic();steps=[];code=1
    try:
        audit=output/'historical-ledger';audit.mkdir()
        for n in ('ledger.json.gz','amended-expectation-violations.json'):
            shutil.copyfile(source/'docs/validation/results-first-bnv-pilot/static-completion/evidence'/n,audit/n)
        calls=[(name,[source/'tests/bnv'/name]) for name in ('test_results_first_cut.py','test_surface_transfer.py','test_results_first_static.py')]
        calls += [('surface-transfer',[source/'tools/linux/results_first/surface_transfer.py','--output',output/'assessment']),
                  ('historical-ledger-audit',[source/'tools/linux/results_first/static_validate.py','--assessment',audit])]
        for label,args in calls:
            remaining=max(1,int(120-(cpu()-start)))
            def limit():resource.setrlimit(resource.RLIMIT_CPU,(remaining,remaining))
            begin=cpu()
            with (output/(label+'.out')).open('x') as out,(output/(label+'.err')).open('x') as err:
                run=subprocess.run([sys.executable,'-I','-B',*map(str,args)],cwd=source,
                                   stdout=out,stderr=err,timeout=150,preexec_fn=limit)
            steps.append(dict(label=label,arguments=list(map(str,args)),exit_code=run.returncode,cpu_seconds=cpu()-begin))
            if run.returncode:raise RuntimeError(label+' failed')
        if cpu()-start>120:raise RuntimeError('aggregate CPU ceiling')
        code=0
    finally:
        receipt=dict(execution_commit=commit,job=os.environ['SLURM_JOB_ID'],host=os.uname().nodename,
            affinity=sorted(os.sched_getaffinity(0)),cpu_seconds=cpu()-start,wall_seconds=time.monotonic()-wall,
            max_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,steps=steps,exit_code=code,
            global_stellar_integrations=0,trajectories=0)
        (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        files={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in output.rglob('*')
               if p.is_file() and not p.name.startswith(('slurm-','allocation-'))}
        (output/'manifest.json').write_text(json.dumps(files,indent=2)+'\n')
    return code


if __name__=='__main__':raise SystemExit(main())
