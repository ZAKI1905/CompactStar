#!/usr/bin/env python3
"""Bounded Slurm runner: algebra tests and retained-array assessment only."""
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


def cpu():
    return time.process_time()+resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Slurm required')
    source = Path(os.environ['RF_SOURCE']).resolve()
    output = Path(os.environ['RF_OUTPUT']).resolve()
    commit = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != os.environ['RF_COMMIT']:
        raise RuntimeError('wrong execution commit')
    if subprocess.check_output(['git', '-C', str(source), 'status', '--porcelain'], text=True):
        raise RuntimeError('source must be clean')
    if output == source or source in output.parents:
        raise RuntimeError('outputs must be external')
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    begin, wall = cpu(), time.monotonic()
    receipts = []
    code = 1
    def call(label, args):
        remaining = max(1, int(60-(cpu()-begin)))
        def limit():
            resource.setrlimit(resource.RLIMIT_CPU, (remaining, remaining))
        c0, t0 = cpu(), time.monotonic()
        with (output/(label+'.out')).open('x') as stdout, (output/(label+'.err')).open('x') as stderr:
            run = subprocess.run([sys.executable, '-I', '-B', *map(str, args)], cwd=source,
                                 stdout=stdout, stderr=stderr, preexec_fn=limit, timeout=120)
        receipts.append(dict(label=label, arguments=list(map(str, args)), exit_code=run.returncode,
                             cpu_seconds=cpu()-c0, wall_seconds=time.monotonic()-t0))
        return run.returncode
    try:
        for test in ('test_results_first_cref', 'test_ba13_final_decision', 'test_results_first_static'):
            if call(test, [source/f'tests/bnv/{test}.py']):
                raise RuntimeError('test failed: '+test)
        code = call('static-assessment', [source/'tools/linux/results_first/static_completion.py', '--output', output/'assessment'])
        if code not in (0, 2):
            raise RuntimeError('assessment execution failed')
        if call('independent-ledger-audit', [source/'tools/linux/results_first/static_validate.py',
                                           '--assessment', output/'assessment']):
            code = 1
            raise RuntimeError('serialized-ledger audit failed')
        if cpu()-begin > 60:
            raise RuntimeError('aggregate CPU ceiling exceeded')
    finally:
        receipt = dict(execution_commit=commit, slurm_job_id=os.environ['SLURM_JOB_ID'],
            host=os.uname().nodename, cpu_seconds=cpu()-begin, wall_seconds=time.monotonic()-wall,
            cpu_ceiling_seconds=60, max_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            exit_code=code, steps=receipts, stellar_calculations=0, trajectories_run=0)
        (output/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
        manifest = {str(p.relative_to(output)): dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size)
                    for p in sorted(output.rglob('*')) if p.is_file() and p.name not in ('external-manifest.json',)
                    and not p.name.startswith(('slurm-', 'allocation-'))}
        (output/'external-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
