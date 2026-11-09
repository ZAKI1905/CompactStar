#!/usr/bin/env python3
"""Independent mpmath recomputation and mutation audit of conditional accounting."""
import argparse
import copy
import gzip
import json
import os
from pathlib import Path
import mpmath as mp

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'docs/validation/results-first-bnv-pilot'

def load(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return json.load(f)

def validate(new, old, bounds, comp, m4):
    checks=0
    def check(ok,label):
        nonlocal checks
        checks+=1
        if not ok:raise ValueError(label)
    def eq(a,b,label):
        check(abs(mp.mpf(str(a))-b)<=mp.mpf('1e-95')*max(abs(b),mp.mpf('1e-100')),label)
    def exact(a,b,label):
        check(json.dumps(a,sort_keys=True)==json.dumps(b,sort_keys=True),label)
    exact({k:v for k,v in new.items() if k not in ('rows','endpoint_acceptance_summary')},
          {k:v for k,v in old.items() if k!='rows'},'historical metadata changed')
    check(len(new['rows'])==len(old['rows'])==462,'row count')
    required=['S','C-REF','C-TAN','C-CONT','C-INPUT','C-SUPPORT','C-CUT','C-MUN','ADR-0022','C-ENDPOINT','METRIC-CONSUMER-SCOPE']
    maximum=mp.mpf(0);increase=mp.mpf(0);groups={};counts={}
    for row,prior in zip(new['rows'],old['rows']):
        exact({k:v for k,v in row.items() if k!='endpoint_acceptance'},prior,'historical typed record changed')
        s=row['endpoint_acceptance'];n=row['quantity'];i=row['index']
        check(s['conditions']==required and s['rigorous_bound'] is False,'condition or evidence class lost')
        counts[s['status']]=counts.get(s['status'],0)+1
        if n not in comp:
            check(s['status']==prior['surface_transfer']['status'] and s['full_utilization']==prior['surface_transfer']['full_utilization'],'scope-only or exact row changed')
            continue
        U=mp.mpf(bounds['allowances'][n]['interior_point_bound'])
        fit=mp.mpf(str(comp[n]['fit_operator_norm']));N=mp.mpf(str(prior['static_completion']['known_numerator']));b=mp.mpf(str(prior['static_completion']['known_baseline']))
        point=(fit+(0 if i==0 else 2))*U
        end=mp.mpf(0)
        if i==20:
            x=list(map(lambda v:mp.mpf(str(v)),m4['achieved_x']))
            end=2*U*max(mp.mpf(0),mp.mpf('1e-6')+x[20])/abs(x[20]-x[10])
        base=mp.mpf(0) if n.startswith('t_') else U
        denominator=1-b-base
        # Independently expanded total, including nonlinear denominator change.
        old_budget=N/(1-b)
        delta_budget=(point+end+old_budget*base)/denominator
        total=old_budget+delta_budget
        threshold=mp.mpf(str(prior['historical']['threshold']))
        util=total/threshold;delta=delta_budget/threshold
        for key,value in [('full_interval_point_bound',U),('additional_numerator',point),
                          ('additional_endpoint',end),('additional_baseline',base),
                          ('full_denominator',denominator),('full_budget',total),
                          ('full_utilization',util),('utilization_increase',delta)]:eq(s[key],value,key)
        if n.startswith('t_'):check(s['baseline_cut_bound']=='NOT_APPLICABLE','absolute tangent normalization')
        else:eq(s['baseline_cut_bound'],U,'relative baseline')
        check(s['status']==('CONDITIONALLY_QUALIFIED' if util<=1 else 'FAIL_CONDITIONAL_BUDGET'),'budget status')
        # The earlier scenario remains unauthorized historical evidence.
        scenario=prior['surface_transfer']['endpoint_extension_scenario']
        check(scenario['authorized'] is False and scenario['qualification'] is False,'historical scenario promoted')
        check(abs(mp.mpf(scenario['utilization'])-util)<mp.mpf('1e-65'),'old scenario numerical disagreement')
        maximum=max(maximum,util);groups[n]=max(groups.get(n,mp.mpf(0)),util)
        if delta>increase:increase=delta;impact=dict(quantity=n,index=i)
    summary=new['endpoint_acceptance_summary']
    check(summary['row_statuses']==counts,'status counts')
    check(summary['status']=='CONDITIONAL_STATIC_PASS','static status')
    eq(summary['maximum_utilization'],maximum,'maximum')
    eq(summary['maximum_cut_utilization_increase'],increase,'maximum increment')
    check(summary['maximum_cut_increase_row']==impact,'impact identity')
    for n,value in groups.items():eq(summary['group_max_utilization'][n],value,'group maximum')
    check(summary['conditions']==required and summary['rigorous_bound'] is False,'summary conditions')
    check(summary['historical_CREF_failures']==42 and summary['original_expectation_violations']==355 and summary['amended_expectation_violations']==353,'historical failures')
    check(summary['actual_path_uncertainty'] is None and summary['all_restricted_pilot_gates_pass'] is False and
          summary['nonlinear_trajectory_qualified'] is False and summary['trajectories_authorized'] is False and
          summary['historical_certificate_promoted'] is False,'unearned readiness promotion')
    return checks

def main():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('EKU Slurm required')
    p=argparse.ArgumentParser();p.add_argument('--assessment',type=Path,required=True);args=p.parse_args()
    new=load(args.assessment/'ledger.json.gz');old=load(BASE/'surface-transfer/evidence/assessment/ledger.json.gz')
    inputs=[load(BASE/'surface-transfer/evidence/assessment/bounds.json'),
            load(ROOT/'docs/validation/ba13-evidence/evidence/comparison.json.gz'),
            load(ROOT/'docs/validation/ba13-closure/evidence/m4.json')]
    mp.mp.dps=130
    checks=validate(new,old,*inputs);mutations={}
    for label in ('baseline','endpoint','cut-zero','threshold','historical-scenario','condition','rigorous','metric','nonlinear-ready','unknown-zero'):
        bad=copy.deepcopy(new);row=next(r for r in bad['rows'] if r['quantity']=='envelope' and r['index']==20);s=row['endpoint_acceptance']
        if label=='baseline':s['additional_baseline']='0'
        elif label=='endpoint':s['additional_endpoint']='0'
        elif label=='cut-zero':s['full_interval_point_bound']='0'
        elif label=='threshold':row['historical']['threshold']*=2
        elif label=='historical-scenario':row['surface_transfer']['endpoint_extension_scenario']['authorized']=True
        elif label=='condition':s['conditions'].remove('C-ENDPOINT')
        elif label=='rigorous':s['rigorous_bound']=True
        elif label=='metric':next(r for r in bad['rows'] if r['quantity']=='metric_structure')['endpoint_acceptance']['status']='CONDITIONALLY_QUALIFIED'
        elif label=='nonlinear-ready':bad['endpoint_acceptance_summary']['all_restricted_pilot_gates_pass']=True
        else:bad['endpoint_acceptance_summary']['actual_path_uncertainty']=0
        try:validate(bad,old,*inputs)
        except ValueError as e:mutations[label]=str(e)
        else:raise RuntimeError('mutation escaped: '+label)
    result=dict(status='PASS_CONDITIONAL_ACCOUNTING_NOT_TRAJECTORY_QUALIFICATION',checks=checks,
        mutations_refused=mutations,arithmetic='mpmath 130 digits, expanded denominator correction, no assessor imports',
        historical_records_preserved=True,stellar_calculations=0,trajectories=0)
    (args.assessment/'independent-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
