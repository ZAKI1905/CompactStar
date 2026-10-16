"""Read-only independent Decimal audit of frozen scoring and raw-output provenance.

No stellar solves, new tolerances, fitted bands or writes to the frozen scorer.
Run under Slurm; charge resource usage alongside the original job.
"""
import hashlib
import json
import os
import resource
from collections import Counter
from decimal import Decimal, localcontext
from pathlib import Path

S=Path(__file__).resolve().parents[3];D=S/'docs/validation/x1-oracle-anchor';O=S.parents[1]/'x1-oracle-anchor'
assert os.environ.get('SLURM_JOB_ID')
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
resource.setrlimit(resource.RLIMIT_CPU,(30,30))
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
gate=read(D/'execution-gate.json')
for path,h in {**gate['files'],**gate['review_files']}.items():assert sha(Path(path))==h,path
rows=read(O/'runs.json');scores=read(O/'scores.json');counts=Counter()
raw={}
fields=['B','R','M','r_e','L','C','rb','mu_inf','Hc']
for row in rows:
    key=(row['variant'],row['target'],row['tol'])
    text=(O/('-'.join(key)+'.out')).read_text();xs=text.split()
    assert list(map(float,xs))==[row['values'][q] for q in fields]
    raw[key]=dict(zip(fields,map(Decimal,xs)))
anchor=json.loads((D/'frozen/anchor.json').read_text(),parse_float=Decimal)['targets']
old=json.loads((D/'frozen/oracle-retained.json').read_text(),parse_float=Decimal)
prod=json.loads((D/'frozen/paired-production.json').read_text(),parse_float=Decimal)['values']
targets=json.loads((D/'targets.json').read_text())
names=dict(L='L_nu',C='C_nu',r_e='r_nu')
with localcontext() as ctx:
    ctx.prec=60
    for i in ['0','10','20']:
        for q in ['L','C','R','M','r_e']:
            vals=[raw[('regularized',i,t)][q] for t in ['2e-10','2e-11','2e-12','2e-13']]
            inc=[b-a for a,b in zip(vals,vals[1:])]
            ratios=[abs(b/a) if a else None for a,b in zip(inc,inc[1:])]
            passed=all(r is not None and r<=1/Decimal('3.3') for r in ratios) and ratios[-1]<=Decimal('.3')
            expected='PASS' if passed else 'FAIL'
            assert scores['A1'][i][q]['status']==expected;counts['A1_'+expected]+=1
            assert scores['A2'][i][q]['status']=='INCONCLUSIVE'
            assert Decimal(str(scores['A2'][i][q]['band']))==anchor[i][q]['band_total_excluding_oracle']
            residual=vals[-1]/anchor[i][q]['predicted_analytic_limit']-1
            counts['anchor_residual_within_original_band']+=abs(residual)<=anchor[i][q]['band_total_excluding_oracle']
            if q!='R':
                move=vals[-1]-old[i+':2e-12'][names.get(q,q)]
                predicted=1 if q=='r_e' else -1
                status='PASS' if move*predicted>0 else 'FAIL' if move*predicted<0 else 'INCONCLUSIVE'
                assert scores['A2s'][i][q]['status']==status;counts['direction_'+status]+=1
            u=abs(inc[-1]/vals[-1])*ratios[-1]/(1-ratios[-1])
            disagreement=abs(raw[('unregularized',i,'2e-13')][q]/vals[-1]-1)
            status='PASS' if disagreement<=10*u else 'FAIL'
            assert scores['independence'][i][q]['status']==status;counts['agreement_'+status]+=1
    for row in rows:
        k=(row['variant'],row['target'],row['tol'])
        assert abs(raw[k]['B']/Decimal(targets[row['target']]['B'])-1)<=Decimal('5e-13')
    for r in scores['A3']:
        q=r['quantity'];i=r['target'];v=r['variant'];t=r['tol'];oq=names.get(q,q)
        drift=raw[(v,i,t)][q]/raw[(v,'0',t)][q]-1
        delta=drift-(old[i+':2e-12'][oq]/old['0:2e-12'][oq]-1)
        status='PASS' if abs(delta)<=Decimal('1e-12') else 'FAIL'
        assert r['oracle_status']==status;counts['paired_oracle_'+status]+=1
        phasepass=True
        for p,pv in prod.items():
            good=abs(drift-(pv[i][q]/pv['0'][q]-1))<=Decimal('1e-11')
            counts['paired_phase_PASS' if good else 'paired_phase_FAIL']+=1;phasepass &= good
        assert r['production_status']==('PASS' if phasepass else 'FAIL')
assert len(rows)==18
u=resource.getrusage(resource.RUSAGE_SELF);used=u.ru_utime+u.ru_stime
prior=read(O/'resources.json')['total_cpu_seconds']
assert prior+used<120
result=dict(status='PASS',purpose='independent arithmetic and provenance audit, not a global error bound',counts=dict(counts),cpu_seconds=used,combined_measured_cpu_seconds=prior+used,max_rss_kib=u.ru_maxrss,affinity=sorted(os.sched_getaffinity(0)),job=os.environ['SLURM_JOB_ID'],source_sha256=sha(Path(__file__)),numerical_solver_runs_added=0)
(O/'independent-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
