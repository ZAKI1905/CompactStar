#!/usr/bin/env python3
"""Bounded endpoint ledger and historical tests on EKU; no stellar calculations."""
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
    root=Path(os.environ['RF_SOURCE']).resolve();out=Path(os.environ['RF_OUTPUT']).resolve()
    commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
    dirty=subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True)
    (out/'source-preflight.json').write_text(json.dumps(dict(actual_commit=commit,
        requested_commit=os.environ['RF_COMMIT'],status=dirty),indent=2)+'\n')
    if dirty or commit!=os.environ['RF_COMMIT']:raise RuntimeError('clean exact source required')
    if root==out or root in out.parents:raise RuntimeError('external output required')
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    start=cpu();wall=time.monotonic();steps=[];code=1
    try:
        # Run old authority tests against their original committed ADR, without
        # rewriting either the historical manifest or the newly amended ADR.
        history=out/'historical-source';history.mkdir()
        # This checkout is sparse/partial: broad git archive requests omitted
        # unrelated blobs. Copy materialized inputs and restore only the pinned
        # historical ADR. Each historical consumer still checks its own hashes.
        for name in ('tools', 'tests/bnv', 'docs/validation'):
            shutil.copytree(root/name,history/name)
        adr='docs/adr/ADR-0022-restricted-bnv-common-reference-residual.md'
        (history/adr).parent.mkdir(parents=True)
        (history/adr).write_bytes(subprocess.check_output([
            'git','-C',str(root),'show','afd4b06074a5f597ef4cbe38ce3c9cb04386f556:'+adr]))
        oldaudit=out/'historical-ledger';oldaudit.mkdir()
        for name in ('ledger.json.gz','amended-expectation-violations.json'):
            shutil.copyfile(root/'docs/validation/results-first-bnv-pilot/static-completion/evidence'/name,oldaudit/name)
        transfer=out/'historical-transfer'
        shutil.copytree(root/'docs/validation/results-first-bnv-pilot/surface-transfer/evidence/assessment',transfer)
        calls=[(name,[history/'tests/bnv'/name]) for name in ('test_results_first_cut.py','test_surface_transfer.py','test_results_first_static.py')]
        calls += [('endpoint-tests',[root/'tests/bnv/test_endpoint_acceptance.py']),
                  ('historical-ledger-audit',[history/'tools/linux/results_first/static_validate.py','--assessment',oldaudit]),
                  ('historical-transfer-audit',[history/'tools/linux/results_first/surface_transfer_validate.py','--assessment',transfer]),
                  ('endpoint-assessment',[root/'tools/linux/results_first/endpoint_acceptance.py','--output',out/'assessment']),
                  ('endpoint-audit',[root/'tools/linux/results_first/endpoint_validate.py','--assessment',out/'assessment'])]
        for name,args in calls:
            remain=max(1,int(120-(cpu()-start)))
            def limit():resource.setrlimit(resource.RLIMIT_CPU,(remain,remain))
            begin=cpu()
            with (out/(name+'.out')).open('x') as stdout,(out/(name+'.err')).open('x') as stderr:
                run=subprocess.run([sys.executable,'-I','-B',*map(str,args)],cwd=root,stdout=stdout,stderr=stderr,
                                   timeout=150,preexec_fn=limit)
            steps.append(dict(name=name,exit_code=run.returncode,cpu_seconds=cpu()-begin))
            if run.returncode:raise RuntimeError(name+' failed')
        if cpu()-start>120:raise RuntimeError('CPU limit')
        code=0
    finally:
        receipt=dict(execution_commit=commit,job=os.environ['SLURM_JOB_ID'],host=os.uname().nodename,
            cpu_seconds=cpu()-start,wall_seconds=time.monotonic()-wall,max_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            affinity=sorted(os.sched_getaffinity(0)),steps=steps,exit_code=code,stellar_calculations=0,trajectories=0,
            historical_source_commit='afd4b06074a5f597ef4cbe38ce3c9cb04386f556')
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        files={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*')
               if p.is_file() and 'historical-source' not in p.parts and not p.name.startswith(('slurm-','allocation-'))}
        (out/'manifest.json').write_text(json.dumps(files,indent=2)+'\n')
    return code
if __name__=='__main__':raise SystemExit(main())
