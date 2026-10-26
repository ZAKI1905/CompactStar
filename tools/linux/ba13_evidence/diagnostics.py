#!/usr/bin/env python3
"""Retained independent-method and correlated derivative consistency evidence."""
import decimal
import json
import os
from reference import S,R,O,F,sha

if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
D=decimal.Decimal;decimal.getcontext().prec=60
def read(p):return json.loads(p.read_text(),parse_float=D)
def ratio(v,k):return (v[1]['N'][k]-v[-1]['N'][k])/(v[1]['primary'][0]-v[-1]['primary'][0])
tg=read(O/'tangents.json');m3=read(F/'m3.json');m1=read(F/'m1.json')
out={}
for k,n in ((0,'t_n'),(2,'t_e'),(3,'t_mu')):
    records={}
    for i,j in ((0,10),(10,20),(0,20)):
        tvals=[]
        for target in (i,j):
            wide={s:read(R/f'ba13-headroom/t{target}-h0.00005-{s}.json') for s in (-1,1)}
            narrow={s:read(R/f'ba13-closure-20261019/t{target}-h0.000025-{s}.json') for s in (-1,1)}
            correction=(ratio(wide,k)-ratio(narrow,k))/3
            tvals.append(D(tg['records'][f'{target}:{n}']['t']['2048'])-correction)
        sec={}
        for resolution in (512,1024,2048):
            base=R/('ba13-headroom' if resolution==512 else 'ba13-closure-20261019')
            a,b=[read(base/f't{target}-n{resolution}.json') for target in (i,j)]
            sec[str(resolution)]=(b['N'][k]-a['N'][k])/(b['primary'][0]-a['primary'][0])
        m4=read(F/'m4.json');dx=m4['achieved_x'][j]-m4['achieved_x'][i]
        # Sign comes from the retained analytic chain-rule stencil curvature.
        curv=[m3['curvatures'][n]['width_checks'][str(target)][-1]['curvature'][0] for target in (i,j)]
        predicted=sum(curv)/2*dx*dx/12
        delta=sum(tvals)/2-sec['2048']
        endpoint_tail=sum(D(tg['records'][f'{target}:{n}']['fine_tail']) for target in (i,j))/2
        sec_increment=abs(sec['2048']-sec['1024'])
        records[f'{i}-{j}']=dict(mean_extrapolated_offset_minus_central_secant=delta,
            predicted_curvature_mean_value_term=predicted,residual=delta-predicted,
            measured_secant_last_increment=sec_increment,offset_mean_tail=endpoint_tail,
            central_secants=sec,
            interpretation='independent stencil geometry on same X2 method; not an independent solver derivative bound',
            unresolved='fourth-width term, higher mean-value derivatives, moving-cut tangent response and arithmetic cancellation')
    roots={str(t):m1['root_alignment'][str(t)] for t in (0,10,20)}
    out[n]=dict(cross_geometry=records,root_alignment=roots,
        independent_method_tangent_reference_available=False,
        classification='offset discretization experimentally supported; full derivative uncertainty remains conditional/unresolved')
abl=S/'docs/validation/x2-fixed-step/evidence/ablations.json'
ab=read(abl)
result=dict(tangents=out,retained_reference_ablations=ab,
    ablation_source=dict(path=str(abl.relative_to(R)),sha256=sha(abl)),
    ablation_scope='primary X2 outputs at retained cases only; no blanket transfer to every passive/differentiated quantity',
    retained_failures_preserved=['double-newton N512 FAIL Newton convergence',
        'historical M1 INCONCLUSIVE','historical M2 INCONCLUSIVE','historical X1/X2 failures'],
    common_bias_enclosure=None,
    mu_n_profile=dict(field='mu_n_inf_MeV',nodes=2049,
        contract='PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT sections 8,9,11; compare_candidate.py profile parser',
        identity='at equilibrium eta=0: mu_n_inf=mu_B_inf; X2 H=log(mu_B/h0), nu=K-H, hence mu_B exp(nu)=h0 exp(K)',
        covers='all source-support nodes for the equilibrium/reference field; current eta-dependent correction remains separate',
        owner_definition_acceptance=False),
    metric_structure=dict(fields=['r_km','m_km','nu','lambda'],nodes=2049,
        independent_nodewise_reference_required=True,observer_benchmarked=False,
        stop='reference uncertainty/continuation/authority prerequisites unclosed; no observer developed',
        dense_output='quadratic collocation interior output is not granted endpoint fourth order'))
print(json.dumps(result,indent=2,default=lambda v:str(v) if isinstance(v,D) else v,allow_nan=False))
