#!/usr/bin/env python3
"""Bounded static reassessment of authenticated BA13 data; never invokes a solver."""
import copy,csv,decimal,gzip,hashlib,json,os,pathlib,sys
import numpy as np
from model import correlated_increment,feasible_cap,complete_total
S=pathlib.Path(__file__).resolve().parents[3]; R=S.parents[1]
O=pathlib.Path(os.environ['BA13_FINAL_OUTPUT_ROOT'])
E=S/'docs/validation/ba13-evidence/evidence'; C=S/'docs/validation/ba13-closure/evidence'
D=decimal.Decimal;decimal.getcontext().prec=60
sources={}
def read(p,dec=False):
    p=pathlib.Path(p);sources[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt') as f:return json.load(f,parse_float=D if dec else float)
def write(name,d):
    with (O/name).open('x') as f:json.dump(d,f,indent=1,allow_nan=False,default=str);f.write('\n')
def main():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
    old=read(E/'matrix.json.gz'); comp=read(R/'ba13-closure-20261019/assessment-final/component-data.json.gz')
    details=read(E/'comparison.json.gz'); m3=read(C/'m3.json');m4=read(C/'m4.json')
    cert=read(R/'ba13-recovery/run1/fit/certificate.json')
    registry=read(S/'docs/validation/ba13-final-decision/acceptance-registry.json')
    checks={}
    def check(n,b):
        checks[n]=bool(b)
        if not b:raise AssertionError(n)
    x=np.array(m4['achieved_x']); knots=x[[20,10,0]]
    comparison={}
    for n,d in details.items():
        P=np.array(comp['P'][n]);X=np.array(d['reference_values'])
        ix=np.column_stack([np.interp(x,knots,X[::-1,c]) for c in range(X.shape[1])])
        err=P-ix;pair=np.maximum.accumulate(np.max(abs(err-err[0]),axis=1))
        check(n+' independent linear interpolator',np.max(abs(ix-np.array(d['interpolated'])))<5e-16)
        check(n+' paired discrepancy',np.max(abs(pair-np.array(d['paired'])))<5e-16)
        check(n+' curvature remains conditional',d['validated_supremum'] is False)
        comparison[n]=dict(max_paired=float(max(pair)),fit_norm=d['fit_operator_norm'],K=d['K'])
    # Recompute every retained tangent ratio from raw decimal endpoint text.
    retained=read(E/'tangents.json'); tangents={}
    for target in (0,10,20):
        raw={}
        for grid in (512,1024,2048):
            for sign in (-1,1):
                p=(R/f'ba13-closure-20261019/t{target}-h0.000025-{sign}.json' if grid==512
                   else E/f't{target}-n{grid}-offset-{sign}.json')
                raw[grid,sign]=read(p,True)
        for n,k in [('t_n',0),('t_e',2),('t_mu',3)]:
            A={g:raw[g,1]['N'][k]-raw[g,-1]['N'][k] for g in (512,1024,2048)}
            B={g:raw[g,1]['primary'][0]-raw[g,-1]['primary'][0] for g in A}
            t={g:A[g]/B[g] for g in A}; d1=t[1024]-t[512];d2=t[2048]-t[1024]
            p=(abs(d1/d2)).ln()/D(2).ln(); oldt=retained['records'][f'{target}:{n}']
            joint=correlated_increment(A[2048]-A[1024],B[2048]-B[1024],t[1024],B[2048])
            check(f'{target}:{n} raw ratio identity',abs(joint-d2)<D('1e-55'))
            check(f'{target}:{n} stored ratio identity',all(abs(t[g]-D(oldt['t'][str(g)]))<D('1e-55') for g in t))
            check(f'{target}:{n} order gate',d1*d2>0 and D('3.5')<=p<=D('4.5'))
            check(f'{target}:{n} width gate',D(3)<=D(str(oldt['stencil_ratio']))<=D(5))
            tangents[f'{target}:{n}']=dict(t={str(g):str(v) for g,v in t.items()},order=str(p),fine_tail=str(abs(d2)/15),
                joint_refinement_identity_residual=str(joint-d2),stencil_tail_observation=oldt['stencil_tail_observation'],
                units='absolute dN_i/dB',evidence_type='validated_numerical_estimate',
                scope='correlated GL4 discretization only; asymptotic error model, not full tangent uncertainty',
                remaining_joint_reference_error=None,remaining_stencil_and_surface_error=None)
    hold=read(E/'holdouts.json'); frozen=read(E/'frozen-holdouts.json')
    check('holdout criteria hash',hold['criterion_sha256']==sources[str(E/'frozen-holdouts.json')])
    for i in ('5','15'):
        check('holdout '+i+' own surface law',hold['targets'][i]['pressure']==frozen['targets'][i]['surfaceP_cgs'])
        check('holdout '+i+' all frozen comparisons',all(v['discrepancy']<=v['criterion'] for v in hold['targets'][i]['groups'].values()))
    # Every gap gets its own inherited coordinate uncertainty, never a global sampled supremum.
    gaps={}
    for n,d in details.items():
        base='mu_B' if n=='mu_n_profile' else n
        gaps[n]=[dict(left=v['left'],right=v['right'],achieved_spacing=v['achieved_spacing'],
                      coordinate_estimate=v['coordinate_estimate'],tested_spacing=v['tested_spacing'],
                      conditional_remainder=d['K']*v['tested_spacing']**2/8,
                      K=d['K'],evidence_type='conditional_estimate',authority=None)
                 for v in m3['curvatures'][base]['intervals']]
        check(n+' all 20 gaps',len(gaps[n])==20)
    rows=[]
    for oldrow in old['rows']:
        row=copy.deepcopy(oldrow);n=row['quantity'];i=row['index'];h=row['historical']
        check(f'{i}:{n} historical identity',h=={k:cert['samples'][i][k][n] for k in h})
        total=complete_total(row['terms'],registry)
        check(f'{i}:{n} complete total refuses missing inputs',total is None if n!='P2_average' else total==0)
        row['prior_evidence_assessment_status']=oldrow['status']
        row['scientific_disposition']='EXACT_IDENTITY_ONLY' if n=='P2_average' else 'UNQUALIFIED'
        row['complete_reference_uncertainty']=0 if n=='P2_average' else None
        row['nonlinear_and_inter_sample_remainder']=0 if n=='P2_average' else None
        row['joint_missing_accuracy_requirement']=None
        row['common_blockers']=([] if n=='P2_average' else ['C-CONT'] if n=='species_support'
            else ['NODEWISE','C-REF','C-CONT','C-INPUT'] if n=='metric_structure'
            else ['C-REF','C-CONT','C-INPUT']+(['C-TAN'] if n.startswith('t_') else []))
        row['production_reference_discrepancy']=next((t['value'] for t in row['terms'] if t['name']=='direct_P_minus_interpolated_X2'),None)
        if n in details:
            d=details[n];row['joint_missing_accuracy_requirement']=feasible_cap(h['threshold'],h['drift_bound'],
                row['finite_added_estimates'],d['fit_operator_norm'],d['base_error_estimate'],absolute=n.startswith('t_'))
        # Do not overwrite prior evidence or totals. Conditional requirements never confer PASS.
        rows.append(row)
    check('462 unique rows',len(rows)==len({(r['index'],r['quantity']) for r in rows})==462)
    group_summary={}
    for n in sorted({r['quantity'] for r in rows}):
        rr=[r for r in rows if r['quantity']==n]; worst=max(rr,key=lambda r:r['finite_added_estimates']/r['headroom'])
        caps=[r['joint_missing_accuracy_requirement']['single_axis_pointwise_cap'] for r in rr if r['joint_missing_accuracy_requirement']]
        group_summary[n]=dict(worst_finite_fraction=worst['finite_added_estimates']/worst['headroom'],
            finite_added=worst['finite_added_estimates'],historical_utilization=worst['historical']['utilization'],
            historical_headroom=min(r['headroom'] for r in rr),partial_remaining=min(r['remaining_headroom_after_finite_terms'] for r in rr),
            min_required_pointwise_cap=min(caps) if caps else None,units=worst['units'],blockers=worst['common_blockers'],
            proposed_total_utilization=None if n!='P2_average' else 0)
    summary=dict(disposition='E',full_qualification=False,owner_accepted=False,
        counts={v:sum(r['status']==v for r in rows) for v in ('PASS','FAIL','UNQUALIFIED')},
        actual_x_endpoint=str(x[-1]),declared_x_endpoint='-0.000001',
        endpoint_shortfall=max(0,float(x[-1]+1e-6)),
        endpoint_sliver_qualified=False,groups=group_summary,
        reason='joint reference/derivative remainder and quantitative continuous/input model incomplete; nodewise metric reference gated',
        stellar_runs=0,checks=len(checks))
    write('matrix.json',dict(summary=summary,rows=rows));write('summary.json',summary)
    write('tangent-reassessment.json',tangents);write('conditional-intervals.json',gaps)
    write('checks.json',checks);write('source-manifest.json',sources)
    cols=['index','quantity','historical_drift','threshold','headroom','production_reference_discrepancy','finite_added_estimates',
          'complete_reference_uncertainty','nonlinear_and_inter_sample_remainder','proposed_total_utilization','status','common_blockers','required_pointwise_cap']
    with (O/'matrix.tsv').open('x') as f:
        w=csv.DictWriter(f,cols,delimiter='\t');w.writeheader()
        for r in rows:
            d={k:r.get(k) for k in cols};d.update(historical_drift=r['historical']['drift_bound'],threshold=r['historical']['threshold'],
                common_blockers=','.join(r['common_blockers']),required_pointwise_cap=(r['joint_missing_accuracy_requirement'] or {}).get('single_axis_pointwise_cap'))
            w.writerow(d)
    print(json.dumps(summary,indent=1))
if __name__=='__main__':main()
