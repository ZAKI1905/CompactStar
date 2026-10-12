#!/usr/bin/env python3
"""Compare Mac/native BA13 components without promoting inherited error budgets."""
import argparse
import json
import math
import sys
from pathlib import Path
from ba13_refit import load, NAMES

def main():
    if sys.flags.optimize:raise RuntimeError('platform comparison requires active assertions')
    p=argparse.ArgumentParser();p.add_argument('--native',type=Path,required=True);p.add_argument('--mac',type=Path,required=True)
    p.add_argument('--native-certificate',type=Path,required=True);p.add_argument('--mac-certificate',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();nx,ns=load(a.native);mx,ms=load(a.mac)
    nc=json.loads(a.native_certificate.read_text());mc=json.loads(a.mac_certificate.read_text())
    records={}
    for name in NAMES:
        nv,ne,scale,kind=ns[name];mv,me,_,_=ms[name]
        rows=[];maxscaled=0.;outside=0;missing=0
        for i in range(21):
            for j,s in enumerate(scale):
                delta=abs(nv[i][j]-mv[i][j]);bound=ne[i][j]+me[i][j]
                # Report abscissa mismatch separately, and a linear alignment
                # diagnostic; it is not a new tolerance or an error certificate.
                dx=nx[i]-mx[i];slope=mc['fits'][name][j]['slope']
                aligned=abs(nv[i][j]-mv[i][j]-slope*dx)
                maxscaled=max(maxscaled,delta/s)
                if delta>bound:outside+=1
                if bound==0 and name not in ('P2_average','species_support'):missing+=1
                rows.append(dict(index=i,component=j,absolute_difference=delta,scaled_difference=delta/s,
                    sum_declared_errors=bound,within_declared_errors=delta<=bound,
                    achieved_abscissa_difference=dx,linear_alignment_diagnostic=aligned,
                    agreement_status='DIAGNOSTIC_ONLY_ERROR_VALIDATION_REQUIRED'))
        fields={key:max(abs(nc['samples'][i][key][name]-mc['samples'][i][key][name]) for i in range(21))
                for key in ('drift_bound','threshold','utilization')}
        records[name]=dict(max_scaled_difference=maxscaled,comparisons=len(rows),outside_declared_error_sum=outside,
            zero_error_sum_without_exact_zero_authority=missing,certificate_max_absolute_differences=fields,components=rows)
    blockers=[
        'Historical fitter sets radius, envelope, mu_n_profile and metric_structure numerical errors to zero; no per-target refinement bounds are serialized.',
        'Historical thermal Cstar, surface gravity and redshifted temperature errors are 64 epsilon times value, not resolution/refinement bounds.',
        'Ltilde/direct-reference errors are order-16 versus order-32 differences; zero differences do not certify zero quadrature/background/EOS error.',
        'Native and historical Phase-5C target adapters inherit B0 background/center characterization; per-target independent refinement is not demonstrated.',
        'Historical bound is cumulative maximum drift plus largest fit residual. Pointwise numerical uncertainty enters thresholds/residuals, not an added envelope term.',
        'Regression uses a 64 epsilon max(1,abs(value)) weight floor. For dimensional tiny Z/Ltilde values this is a conditioning floor, not physical uncertainty.'
    ]
    report=dict(status='UNQUALIFIED_PLATFORM_AGREEMENT',contract='PROPOSED ADR-0020 B: independently validated summed componentwise errors',
        fit_semantics='Unchanged authenticated historical assembler; discrepancies documented without threshold/window changes',
        max_achieved_abscissa_difference=max(abs(a-b) for a,b in zip(nx,mx)),quantities=records,blockers=blockers,
        owner_approval=False,full_campaign_authorized=False)
    with a.output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(report['status'])
    for name,r in records.items():print(name,'max scaled difference',r['max_scaled_difference'],'outside declared errors',r['outside_declared_error_sum'],'of',r['comparisons'])

if __name__=='__main__':main()
