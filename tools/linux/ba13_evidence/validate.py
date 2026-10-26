#!/usr/bin/env python3
"""Independent successor recomputation and executable source mutation controls."""
import gzip
import hashlib
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import numpy as np

S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1]
O=pathlib.Path(os.environ.get('BA13_EVIDENCE_OUTPUT_ROOT',R/'ba13-evidence-20261026'))
F=R/'ba13-closure-20261019/assessment-final'


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
    checks={};tests=[];mutations={}
    def check(name,condition):
        checks[name]=bool(condition)
        if not condition:raise AssertionError(name)
    out=O/'tests';out.mkdir()
    for script in sorted((S/'tests/bnv').glob('test_ba13_*.py')):
        q=subprocess.run([sys.executable,'-B',str(script)],capture_output=True,text=True)
        (out/(script.stem+'.out')).write_text(q.stdout+q.stderr)
        tests.append(dict(script=script.name,returncode=q.returncode,output_sha256=sha(out/(script.stem+'.out'))))
        check(script.name,q.returncode==0)
    mx=read(O/'matrix.json');cert=read(R/'ba13-recovery/run1/fit/certificate.json')
    prior=read(F/'matrix.json');old={(r['index'],r['quantity']):r for r in prior['rows']}
    rows={(r['index'],r['quantity']):r for r in mx['rows']}
    check('462 unique rows',len(rows)==len(mx['rows'])==462)
    check('historical values and statuses preserved',all(r['historical']==old[k]['historical'] and
          r['historical_assessment_status']==old[k]['status'] for k,r in rows.items()))
    check('certificate identities preserved',all(r['historical']=={n:cert['samples'][k[0]][n][k[1]]
          for n in ('drift_bound','threshold','utilization')} for k,r in rows.items()))
    check('only original exact P2 passes',all((r['status']=='PASS')==(k[1]=='P2_average') for k,r in rows.items()))
    check('incomplete scientific totals remain null',all(r['proposed_total_utilization'] is None and
          r['remaining_physical_headroom'] is None and r['added_justified_error'] is None
          for r in rows.values() if r['status']=='UNQUALIFIED'))
    check('no implicit nonzero analytic labels',all(t['evidence_type']!='analytic_bound' or t['value']==0
          for r in rows.values() for t in r['terms']))
    check('all nonzero sources explicit',all(t['provenance'] and t['evidence_type']
          for r in rows.values() for t in r['terms'] if t['value']))
    check('all missing terms have reasons',all(t['missing'] for r in rows.values() for t in r['terms'] if t['value'] is None))
    check('no invented authority',all(t['authority'] is None for r in rows.values() for t in r['terms']))
    with gzip.open(F/'component-data.json.gz','rt') as f:comp=json.load(f)
    details=read(O/'comparison.json');tg=read(O/'tangents.json');m4=read(F/'m4.json')
    x=np.array(m4['achieved_x']);knots=x[[20,10,0]]
    for n,d in details.items():
        if n.startswith('t_'):
            X=np.array([[float(tg['records'][f'{t}:{n}']['t']['2048'])] for t in (20,10,0)])
        elif n=='mu_n_profile':
            scale=comp['original_units']['mu_B']['scales'][0]/np.array(comp['original_units'][n]['scales'])
            X=np.array([comp['X'][str(t)]['mu_B'] for t in (20,10,0)])*scale
        else:X=np.array([comp['X'][str(t)][n] for t in (20,10,0)])
        expected=np.column_stack([np.interp(x,knots,X[:,c]) for c in range(X.shape[1])])
        check(n+' reference is independently bound to X2',np.max(abs(expected-np.array(d['interpolated'])))<4e-16)
        P=np.array(comp['P'][n]);e=P-expected
        pair=np.maximum.accumulate(np.max(abs(e-e[0]),axis=1))
        check(n+' all21 direct comparison',np.max(abs(pair-np.array(d['paired'])))<5e-16)
        check(n+' nonzero root charged',d['root_pair']>0 and all(any(t['name']=='X2_common_B_alignment' and
              t['value']==d['root_pair'] for t in rows[(i,n)]['terms']) for i in range(1,21)))
        check(n+' full numerator normalization',all(
            abs(next(t['value'] for t in rows[(i,n)]['terms'] if t['name']=='full_numerator_normalization')-
                (rows[(i,n)]['historical']['drift_bound']+sum(t['value'] or 0 for t in rows[(i,n)]['terms']
                    if t['name']!='full_numerator_normalization'))*d['base_error_estimate']/(1-d['base_error_estimate']))<1e-25
            for i in range(21)) if not n.startswith('t_') else True)
        check(n+' all21 fit norm includes inner half',d['fit_operator_norm']>5.)
    hold=read(O/'holdouts.json');frozen=read(O/'frozen-holdouts.json')
    check('both withheld targets scored',set(hold['targets'])=={'5','15'})
    check('frozen holdout criteria unchanged',hold['criterion_sha256']==sha(O/'frozen-holdouts.json') and
          sha(O/'frozen-holdouts.json')==sha(S/'docs/validation/ba13-evidence/evidence/frozen-holdouts.json'))
    check('all 30 holdout group comparisons pass',all(v['discrepancy']<=v['criterion'] for t in hold['targets'].values() for v in t['groups'].values()))
    check('holdouts use own cutoffs',all(hold['targets'][str(i)]['pressure']==frozen['targets'][str(i)]['surfaceP_cgs'] for i in (5,15)))
    check('all 9 tangent estimates meet frozen order gate',len(tg['records'])==9 and all(v['order_ok'] and v['width_ok'] for v in tg['records'].values()))
    check('tangent correlated identity',max(abs(float(v['correlated_identity_residual'])) for v in tg['records'].values())<1e-55)
    check('old M1 M2 bytes remain unchanged',all(subprocess.run(['git','-C',str(S),'diff','--quiet',
        '093d121135e35bc50ca2a807470906bcb28ef5c3','--',f'docs/validation/ba13-closure/evidence/{n}.json']).returncode==0 for n in ('m1','m2')))
    # True mutations: alter executable assessor modules, run the same test suite,
    # and retain each failing receipt. These are intended failures, not gate tuning.
    variants={
      'wrong-reference-accepted':('reference.py',"if reference_kind != 'X2':","if False:"),
      'comparison-reference-omitted':('reference.py','np.asarray(production,dtype=LD)-np.asarray(reference,dtype=LD)',
                                     'np.asarray(production,dtype=LD)-np.asarray(production,dtype=LD)'),
      'root-charge-omitted':('reference.py','return error, paired, float(root_band)','return error, paired, 0.0'),
      'acceptance-bypass':('evidence.py',"if not t['authority'] or not acceptance or acceptance!=t['authority']:",'if False:'),
      'evidence-promoted':('evidence.py','value=value,evidence_type=evidence_type,',"value=value,evidence_type='analytic_bound',"),
      'historical-resolution-rule':('evidence.py','if reference_band is None or not complete or not accepted:',
                                  'if reference_band is None or reference_band>=discrepancy or not complete or not accepted:')}
    for name,(file,needle,replacement) in variants.items():
        with tempfile.TemporaryDirectory(prefix='ba13-mutation-',dir=O) as tmp:
            tmp=pathlib.Path(tmp)
            for mod in ('evidence.py','reference.py'):
                code=(S/'tools/linux/ba13_evidence'/mod).read_text()
                if mod==file:
                    check(name+' mutation applied',needle in code)
                    code=code.replace(needle,replacement)
                (tmp/mod).write_text(code)
            env=dict(os.environ,BA13_EVIDENCE_MODULES=str(tmp))
            p=subprocess.run([sys.executable,'-B',str(S/'tests/bnv/test_ba13_evidence.py')],capture_output=True,text=True,env=env)
            (out/(name+'.out')).write_text(p.stdout+p.stderr)
            mutations[name]=dict(returncode=p.returncode,detected=p.returncode!=0,
                                output_sha256=sha(out/(name+'.out')),
                                source_sha256=sha(tmp/file))
            check(name+' detected',p.returncode!=0 and 'FAIL' in p.stderr)
    result=dict(status='PASS',checks=checks,check_count=len(checks),tests=tests,mutations=mutations,
                source_commit=subprocess.check_output(['git','-C',str(S),'rev-parse','HEAD'],text=True).strip(),
                assessor_sha256=sha(S/'tools/linux/ba13_evidence/assemble.py'),
                decision_sha256=sha(S/'tools/linux/ba13_evidence/evidence.py'),
                reference_sha256=sha(S/'tools/linux/ba13_evidence/reference.py'),
                matrix_sha256=sha(O/'matrix.json'),scientific_qualification=False)
    with (O/'validation.json').open('x') as f:json.dump(result,f,indent=1);f.write('\n')
    print(json.dumps(dict(status='PASS',checks=len(checks),tests=43,mutations=len(mutations))))


if __name__=='__main__':main()
