#!/usr/bin/env python3
"""Necessary static C-REF check; cannot issue an operational certificate."""
import argparse
import copy
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import time

REPO = Path(__file__).resolve().parents[3]
DOC = REPO / 'docs/validation/results-first-bnv-pilot'
AUTHORITY = 'owner:RF-BNV-P2-20261108:4509e13217ec442a2649ab0e2dcf950f3c0f5c97'

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def authenticate(registry, document_root=DOC):
    if registry['scope'] != 'RF-BNV-P2-NATIVE-BOUNDED-v1' or registry['status'] != 'ACCEPTED_FOR_PREPARATION_ONLY':
        raise ValueError('wrong preparation authority scope')
    if registry['paired_trajectories_authorized'] or registry['historical_certificate_promoted']:
        raise ValueError('preparation authority cannot launch or promote')
    for key, name in [('approval_sha256','OWNER_APPROVAL.md'),('decision_sha256','OWNER_DECISION.md'),
                      ('review_sha256','evidence/INDEPENDENT_REVIEW.md')]:
        if digest(document_root/name) != registry[key]:
            raise ValueError('changed approval evidence: '+name)
    for key in ('C-REF','C-TAN','C-CONT','C-INPUT','C-SUPPORT','C-CUT','C-MUN'):
        if registry['accepted_conditions'].get(key) != AUTHORITY:
            raise ValueError('unaccepted condition '+key)
    if registry['accepted_standards'].get('S') != AUTHORITY:
        raise ValueError('unaccepted numerical standard')

def partial_budget(A, b, f, e, *, absolute=False, origin=False):
    if not all(math.isfinite(x) and x >= 0 for x in (A,b,f,e)) or b+e >= 1:
        raise ValueError('invalid budget or unseparated denominator')
    coefficient = f if origin else 2+f
    return A+coefficient*e if absolute else (A*(1-b)+coefficient*e)/(1-b-e)

def assess(matrix, comparison):
    result=copy.deepcopy(matrix)
    failed=[];expectation_failed=[];maxima={}
    for row in result['rows']:
        n=row['quantity'];i=row['index']
        if n not in comparison:
            row['restricted_preparation']={'status':'CONSUMER_SCOPE_ONLY' if n=='metric_structure' else 'PENDING_FULL_COMPLETION'}
            continue
        detail=comparison[n]
        agreement=max(abs(v) for component in detail['signed_error'] for v in component)
        e=10*(9.1e-12 if n=='Ltilde_Me' else agreement)
        T=row['historical']['threshold']
        A=row['historical']['drift_bound']+row['finite_added_estimates']
        total=partial_budget(A,detail['base_error_estimate'],detail['fit_operator_norm'],e,
                             absolute=n.startswith('t_'),origin=i==0)
        u=total/T;change=u-row['partial_utilization']
        candidate=dict(status='FAIL_CREF_POINTWISE_COMPLETION' if u>1 else 'PENDING_FULL_COMPLETION',
            interpretation='pointwise absolute production-X2 agreement; LMe accepted exception',
            max_absolute_agreement=agreement,common_residual_allowance=e,
            budget_lower_bound=total,utilization_lower_bound=u,
            increase_from_finite_ledger=change,full_conditional_utilization=None,
            unresolved_terms_remain=True,operational_certificate=False)
        row['restricted_preparation']=candidate
        if u>1:failed.append({'index':i,'quantity':n,'utilization_lower_bound':u})
        if change>1e-3:expectation_failed.append({'index':i,'quantity':n,'increase':change})
        if n not in maxima or u>maxima[n]['utilization_lower_bound']:maxima[n]=dict(index=i,**candidate)
    # Existing rows, fields and dispositions must remain exactly equal.
    restored=copy.deepcopy(result)
    for row in restored['rows']:del row['restricted_preparation']
    if restored != matrix:raise RuntimeError('historical evidence was altered')
    result['restricted_summary']=dict(status='BLOCKED_CREF_INTERPRETATION' if failed else 'PENDING_FULL_COMPLETION',
        failed_rows=failed,review_expectation_failures=expectation_failed,group_maxima=maxima,
        historical_rows_preserved=len(matrix['rows']),operational_certificate_issued=False,
        trajectories_authorized=False,trajectories_run=0)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required for assessment')
    resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    begin=time.process_time();wall=time.monotonic()
    registry=json.loads((DOC/'preparation-registry.json').read_text());authenticate(registry)
    files=[REPO/'docs/validation/ba13-final-decision/evidence/matrix.json.gz',
           REPO/'docs/validation/ba13-evidence/evidence/comparison.json.gz']
    inputs=[]
    for p in files:
        with gzip.open(p,'rt') as f:inputs.append(json.load(f))
    result=assess(*inputs)
    result['preparation_provenance']={str(p.relative_to(REPO)):digest(p) for p in
        files+[DOC/'preparation-registry.json',DOC/'CREF_COMPLETION_MODEL.md',Path(__file__)]}
    result['preparation_provenance']['execution_commit']=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip()
    summary=dict(result['restricted_summary'],provenance=result['preparation_provenance'],
                 slurm_job_id=os.environ['SLURM_JOB_ID'],cpu_seconds=time.process_time()-begin,
                 wall_seconds=time.monotonic()-wall)
    a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'sidecar.json.gz').open('xb') as f:
        f.write(gzip.compress(json.dumps(result,sort_keys=True,allow_nan=False).encode(),mtime=0))
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:summary[k] for k in ('status','historical_rows_preserved','operational_certificate_issued','cpu_seconds')}))
    print('Failing rows:',len(summary['failed_rows']))
    # Scientific failure is represented by nonzero exit, after evidence retention.
    return 2 if summary['failed_rows'] else 3

if __name__=='__main__':raise SystemExit(main())
