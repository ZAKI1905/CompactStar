"""Frozen X1 scoring; numerical results cannot change any acceptance constant."""
import json
import math
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[3]
DOC = SRC / 'docs/validation/x1-oracle-anchor'
OUT = Path(sys.argv[1])
read = lambda p: json.loads(p.read_text())
anchor = read(DOC / 'frozen/anchor.json')['targets']
old_raw = read(DOC / 'frozen/oracle-retained.json')
production = read(DOC / 'frozen/paired-production.json')['values']
targets = read(DOC / 'targets.json')
rows = read(OUT / 'runs.json')
qs = ['L', 'C', 'R', 'M', 'r_e']
allq = qs + ['mu_inf', 'B']
mapping = dict(L='L_nu', C='C_nu', r_e='r_nu')
old = {key: {q: value[mapping.get(q, q)] for q in allq} for key, value in old_raw.items()}
data = {(r['variant'], r['target'], r['tol']): r['values'] for r in rows}
def relative(a, b):
    return (a-b)/b
def ladder(seq, q):
    vals = [s[q] for s in seq]
    inc = [(b-a)/vals[-1] for a,b in zip(vals, vals[1:])]
    ratios = [abs(b/a) if a else None for a,b in zip(inc, inc[1:])]
    monotone = all(v > 0 for v in inc) or all(v < 0 for v in inc)
    return dict(values=vals, increments_relative=inc, ratios=ratios, monotone=monotone)

report = dict(status='INCONCLUSIVE', rigorous_global_bound=None, A1={}, A2={}, A2s={}, A3=[], independence={}, ladders={}, alignment=[])
complete = len(rows)==18
alignment = True
for r in rows:
    residual=relative(r['values']['B'], float(targets[r['target']]['B']))
    alignment &= abs(residual)<=5e-13
    report['alignment'].append(dict(variant=r['variant'],target=r['target'],tol=r['tol'],residual=residual,pass_alignment=abs(residual)<=5e-13))

a1 = complete
for i in ['0','10','20']:
    reg=[data.get(('regularized',i,t)) for t in ['2e-10','2e-11','2e-12','2e-13']]
    unreg=[old[i+':2e-10'],data.get(('unregularized',i,'2e-11')),old[i+':2e-12'],data.get(('unregularized',i,'2e-13'))]
    report['ladders'][i]={}
    for variant,seq in [('regularized',reg),('unregularized',unreg)]:
        if all(seq):report['ladders'][i][variant]={q:ladder(seq,q) for q in allq}
    report['A1'][i]={}
    for q in qs:
        row=report['ladders'][i].get('regularized',{}).get(q)
        passed=bool(row and all(v is not None and v<=1/3.3 for v in row['ratios']) and row['ratios'][-1]<=.3)
        u=None
        if row:
            rho=row['ratios'][-1]
            if rho is not None and rho<1:u=abs(row['increments_relative'][-1])*rho/(1-rho)
        report['A1'][i][q]=dict(status='PASS' if passed else ('FAIL' if row else 'INCONCLUSIVE'),U_emp=u,conditional_only=True)
        a1 &= passed

independence=complete
for i in ['0','10','20']:
    r=data.get(('regularized',i,'2e-13'));u=data.get(('unregularized',i,'2e-13'))
    report['independence'][i]={}
    if r and u:
        for q in allq:
            diff=relative(u[q],r[q]);emp=report['A1'][i].get(q,{}).get('U_emp')
            passed=emp is not None and abs(diff)<=10*emp
            report['independence'][i][q]=dict(unregularized_minus_regularized_rel=diff,limit_10U=None if emp is None else 10*emp,status=('PASS' if passed else 'FAIL') if q in qs and emp is not None else 'INCONCLUSIVE')
            if q in qs:independence &= passed

gate=a1 and alignment and independence
for i in ['0','10','20']:
    final=data.get(('regularized',i,'2e-13'))
    report['A2'][i]={};report['A2s'][i]={}
    if not final:continue
    for q in qs:
        pred=anchor[i][q]['predicted_analytic_limit'];band=anchor[i][q]['band_total_excluding_oracle']
        diff=relative(final[q],pred);emp=report['A1'][i][q]['U_emp']
        eligible=gate and emp is not None and emp<=band
        report['A2'][i][q]=dict(status=('PASS' if abs(diff)<=band+emp else 'FAIL') if eligible else 'INCONCLUSIVE',residual_rel=diff,band=band,U_emp=emp,reason='scored empirically' if eligible else 'A1/completion/independence gate or excessive/undefined U_emp')
        if q!='R':
            move=relative(final[q],old[i+':2e-12'][q]);sign=1 if q=='r_e' else -1
            report['A2s'][i][q]=dict(status='PASS' if move*sign>0 else ('FAIL' if move*sign<0 else 'INCONCLUSIVE'),change_from_retained_2e12_rel=move,predicted_sign=sign,decisive_absolute_test=eligible)

for variant,tols in [('unregularized',['2e-11','2e-13']),('regularized',['2e-10','2e-11','2e-12','2e-13'])]:
    for t in tols:
        base=data.get((variant,'0',t))
        for i in ['10','20']:
            val=data.get((variant,i,t))
            if not base or not val:continue
            for q in qs+['mu_inf']:
                drift=relative(val[q],base[q]);ref=relative(old[i+':2e-12'][q],old['0:2e-12'][q])
                phase_diffs={p:drift-relative(v[i][q],v['0'][q]) for p,v in production.items()}
                report['A3'].append(dict(variant=variant,tol=t,target=i,quantity=q,drift=drift,retained_drift=ref,solver_difference=drift-ref,oracle_status='PASS' if abs(drift-ref)<=1e-12 else 'FAIL',phase_differences=phase_diffs,production_status='PASS' if len(phase_diffs)==6 and max(map(abs,phase_diffs.values()))<=1e-11 else 'FAIL'))
report['gates']=dict(complete=complete,alignment=alignment,A1=a1,independence=independence,A2_eligible_global=gate)
if gate:
    statuses=[v['status'] for t in report['A2'].values() for v in t.values()]
    report['status']='FAIL' if 'FAIL' in statuses else ('INCONCLUSIVE' if 'INCONCLUSIVE' in statuses else 'EMPIRICAL_PASS_ONLY')
# Distinguish target-zero common offsets from paired differential residuals.
report['anchor_error_decomposition']={}
if all(data.get(('regularized',i,'2e-13')) for i in ['0','10','20']):
    for q in qs:
        e={i:report['A2'][i][q]['residual_rel'] for i in ['0','10','20']}
        report['anchor_error_decomposition'][q]=dict(absolute_relative_to_frozen_anchor=e,common_target0=e['0'],target_minus_target0={i:e[i]-e['0'] for i in ['10','20']},note='anchor discrepancy is not a known exact numerical error')
(OUT/'scores.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps(report['gates']))
