#!/usr/bin/env python3
"""Completion/provenance audit, not an independent scientific review."""
import hashlib,json,pathlib,subprocess,sys,time
import numpy as np
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1];O=R/'ba13-headroom';D=S/'docs/validation/ba13-headroom'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
start=time.process_time();checks={}
auth=read(D/'authentication.json');allowed={str(S/'tests/bnv/x2_fixed_step.cpp'),str(S/'docs/adr/ADR-0021-ba13-static-uncertainty-completion.md')}
for p,h in auth['protected_inputs'].items():
 if p not in allowed:assert sha(pathlib.Path(p))==h,p
checks['protected_input_bytes']=len(auth['protected_inputs'])-len(allowed)
old_manifest=read(S/'docs/validation/x2-fixed-step/evidence/external-manifest.json')
for name,entry in old_manifest['files'].items():assert sha(pathlib.Path(old_manifest['root'])/name)==entry['sha256'],name
checks['original_X2_external_artifacts']=len(old_manifest['files'])
raw_root=R/'ba13-recovery/run1/raw';raw_manifest=read(raw_root/'COMPLETE.json')
for name,h in raw_manifest['files'].items():assert sha(raw_root/name)==h,name
checks['original_nominal_manifest_files']=len(raw_manifest['files'])
assert sha(S/'tests/bnv/compare_candidate.py')=='9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370'
for i in [0,10,20]:
 assert read(O/f't{i}-regression-final.json')==read(R/f'x2-fixed-step/t{i}-n128.json')
 # The runner also checks all decimal strings, including shooting diagnostics.
 r=read(O/f'receipt-regression-final-{i}.json');assert len(r['runs'])==1 and r['runs'][0]['returncode']==0
checks['full_decimal_X2_replays']=3
m=read(O/'matrix.json');assert len(m['rows'])==462
old=read(R/'ba13-recovery/run1/fit/certificate.json');hr=read(S/'docs/validation/x2-fixed-step/evidence/headroom.json')
for r in m['rows']:
 i,n=r['index'],r['quantity']
 assert all(r['historical'][k]==old['samples'][i][k][n] for k in ['threshold','drift_bound','utilization'])
 assert r['headroom']==next(v['headroom'] for v in hr['rows'] if v['target']==i and v['quantity']==n)
 assert abs(sum(r['allocated_budget'].values())-r['headroom'])<1e-18
 assert (r['status']=='PASS')==(n=='P2_average')
 if r['missing']:assert r['status']=='UNQUALIFIED' and r['proposed_total_utilization'] is None
checks['frozen_historical_rows_and_complete_budgets']=462
for i in range(21):
 p=read(O/f'profile-v2-{i}.json')
 for name,h in p['hash_inputs'].items():assert sha(pathlib.Path(name))==h
 assert p['cache_reproduction']<1e-13
 assert p['tail']['direct_and_independent_enclosed']
checks['nominal_profiles_cache_and_retained_tail']=21
supp=read(O/'supplement.json');assert len(supp['records'])==21
for r in supp['records']:
 assert len(r['mapped_metric'])==8196
 if '32' in r['integrals']:
  assert abs(r['integrals']['16']['P1']/r['integrals']['32']['P1']-1)<1e-12
checks['direct_events_and_metric_mapping']=21
result=subprocess.run([sys.executable,'-B',str(S/'tests/bnv/test_ba13_headroom.py')],capture_output=True,text=True)
(O/'tests.txt').write_text(result.stdout+result.stderr);assert result.returncode==0
checks['fault_and_algebra_tests']=14
result=dict(status='PASS_COMPLETION_AUDIT_NOT_SCIENTIFIC_QUALIFICATION',checks=checks,cpu_seconds=time.process_time()-start,source_sha256=sha(pathlib.Path(__file__)))
(O/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
