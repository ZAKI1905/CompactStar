#!/usr/bin/env python3
"""Independent Decimal accounting and refusal checks for this one local witness."""
import argparse
from collections import Counter
import copy
from decimal import Decimal, localcontext
import gzip
import json
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]


def load(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:
        return json.load(f)


def D(x):return Decimal(str(x))


def validate(new,old,comparison,m4,bounds):
    count=0
    def check(value,label):
        nonlocal count
        count+=1
        if not value:raise ValueError(label)
    def close(a,b,label):
        check(abs(D(a)-b)<=D('1e-60')*max(abs(b),D('1e-100')),label)
    check({k:v for k,v in new.items() if k not in ('rows','surface_transfer_summary')}==
          {k:v for k,v in old.items() if k!='rows'},'original metadata changed')
    check(len(new['rows'])==len(old['rows'])==462,'row count')
    summary=new['surface_transfer_summary']
    check(summary['qualified'] is False and summary['full_C_CUT_budget_effect'] is None,'unsupported qualification')
    check(summary['historical_CREF_failures']==42 and summary['expectation_violations']==353,'historical counts')
    check(bounds['endpoint_extension_supported'] is False and D(bounds['eta'])>0,'endcap dependency erased')
    maximum=Decimal(0);increase=Decimal(0);impact_row=None
    for row,prior in zip(new['rows'],old['rows']):
        check({k:v for k,v in row.items() if k!='surface_transfer'}==prior,'original typed record changed')
        n,i=row['quantity'],row['index'];s=row['surface_transfer'];p=prior['static_completion']
        if n not in comparison:
            check(s['status']==p['status'] and s['full_utilization']==p.get('full_utilization'),'independent gate changed')
            continue
        check(s['status']=='UNQUALIFIED_ENDPOINT_EXTENSION','unsupported status')
        check(s['full_budget'] is None and s['full_utilization'] is None and
              s['full_interval_point_bound'] is None,'endcap promoted')
        check(s['baseline_cut_bound']==('NOT_APPLICABLE' if n.startswith('t_') else None),'baseline promoted')
        U=D(bounds['allowances'][n]['interior_point_bound'])
        check(U>0 and D(s['interior_point_bound'])==U,'interior bound changed')
        e=s['endpoint_extension_scenario']
        check(e['authorized'] is False and e['qualification'] is False,'scenario promoted')
        add_n=(D(comparison[n]['fit_operator_norm'])+(0 if i==0 else 2))*U
        add_e=Decimal(0)
        if i==20:
            xx=list(map(D,m4['achieved_x']))
            add_e=2*U*max(Decimal(0),D('1e-6')+xx[20])/abs(xx[20]-xx[10])
        add_b=Decimal(0) if n.startswith('t_') else U
        N=D(p['known_numerator']);b=D(p['known_baseline']);threshold=D(row['historical']['threshold'])
        total=(N+add_n+add_e)/(1-b-add_b);util=total/threshold
        delta=(total-N/(1-b))/threshold
        for key,value in (('additional_numerator',add_n),('additional_endpoint',add_e),
                          ('additional_baseline',add_b),('full_budget',total),
                          ('utilization',util),('utilization_increase',delta)):
            close(e[key],value,key+' arithmetic')
        check(util<1,'conditional scenario exceeds historical cap')
        maximum=max(maximum,util)
        if delta>increase:increase=delta;impact_row=dict(quantity=n,index=i)
    close(summary['endpoint_extension_scenario_max_utilization'],maximum,'maximum utilization')
    close(summary['endpoint_extension_scenario_max_utilization_increase'],increase,'maximum increment')
    check(summary['endpoint_extension_scenario_max_increase_row']==impact_row,'maximum increment identity')
    return count


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('EKU Slurm required')
    p=argparse.ArgumentParser();p.add_argument('--assessment',type=Path,required=True);args=p.parse_args()
    new=load(args.assessment/'ledger.json.gz')
    old=load(ROOT/'docs/validation/results-first-bnv-pilot/static-completion/evidence/ledger.json.gz')
    comp=load(ROOT/'docs/validation/ba13-evidence/evidence/comparison.json.gz')
    m4=load(ROOT/'docs/validation/ba13-closure/evidence/m4.json')
    bounds=load(args.assessment/'bounds.json')
    with localcontext() as context:
        context.prec=110
        count=validate(new,old,comp,m4,bounds);mutations={}
        for label in ('drop-baseline','drop-final-endpoint','promote-scenario','promote-endcap','change-historical-threshold'):
            bad=copy.deepcopy(new)
            row=next(r for r in bad['rows'] if r['index']==20 and r['quantity']=='envelope')
            s=row['surface_transfer'];e=s['endpoint_extension_scenario']
            if label=='drop-baseline':e['additional_baseline']='0'
            elif label=='drop-final-endpoint':e['additional_endpoint']='0'
            elif label=='promote-scenario':e['authorized']=True
            elif label=='promote-endcap':s['full_budget']=e['full_budget']
            else:row['historical']['threshold']*=2
            try:validate(bad,old,comp,m4,bounds)
            except ValueError as error:mutations[label]=str(error)
            else:raise RuntimeError('mutation escaped: '+label)
    result=dict(status='PASS_ACCOUNTING_AND_REFUSAL_NOT_QUALIFICATION',checks=count,
                mutations_refused=mutations,row_statuses=dict(Counter(r['surface_transfer']['status'] for r in new['rows'])),
                historical_ledger_preserved=True,arithmetic='Decimal 110 digits; no assessor imports',
                global_stellar_integrations=0,trajectories=0)
    (args.assessment/'independent-transfer-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
