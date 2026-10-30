#!/usr/bin/env python3
"""Verify immutable inputs, complete inventory, provenance and packaged outputs."""
import gzip,hashlib,json,pathlib,subprocess,sys
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1]
P=S/'docs/validation/ba13-final-decision';E=P/'evidence'
BASE='84e1062db78b91e4fb0988074941b1aea6019917'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return json.load(f)
def git(*a):return subprocess.check_output(['git','-C',str(S),*a],text=True).strip()
checks={}
def check(n,b):
    checks[n]=bool(b)
    if not b:raise AssertionError(n)
obj=subprocess.check_output(['git','-C',str(S),'cat-file','commit',BASE])
h=subprocess.check_output(['git','-C',str(S),'hash-object','-t','commit','--stdin'],input=obj,text=False).decode().strip()
check('full starting SHA recomputed from commit object',h==BASE)
check('entry branch remains unchanged',git('rev-parse','candidate/linux-ba13-evidence-20261026')==BASE)
check('descendant ancestry',subprocess.run(['git','-C',str(S),'merge-base','--is-ancestor',BASE,'HEAD']).returncode==0)
auth=read(E/'authentication.json')
master=next(b['sha'] for b in auth['branches'] if b['branch']=='master')
check('canonical master unchanged',git('rev-parse','master')==master)
paths=git('diff',BASE,'--name-only').splitlines()
allowed=('docs/validation/ba13-final-decision/','tools/linux/ba13_final_decision/',
         'tests/bnv/test_ba13_final_decision.py','docs/adr/ADR-0021-ba13-static-uncertainty-completion.md')
check('only authorized new validation and ADR paths changed',all(p.startswith(allowed) for p in paths))
oldadr=subprocess.check_output(['git','-C',str(S),'show',BASE+':docs/adr/ADR-0021-ba13-static-uncertainty-completion.md'])
check('ADR history append only',(S/'docs/adr/ADR-0021-ba13-static-uncertainty-completion.md').read_bytes().startswith(oldadr))
check('predeclaration immutable',subprocess.check_output(['git','-C',str(S),'show','d5b40b3:docs/validation/ba13-final-decision/PREDECLARATION.md'])==(P/'PREDECLARATION.md').read_bytes())
for p,h in auth['identities'].items():
    if p.endswith('ADR-0021-ba13-static-uncertainty-completion.md'):continue
    check('authenticated identity '+p,sha(pathlib.Path(p))==h)
old=read(S/'docs/validation/ba13-evidence/evidence/matrix.json.gz');new=read(E/'matrix.json.gz')
check('462 rows',len(new['rows'])==462)
provenance={}
for a,b in zip(old['rows'],new['rows']):
    check('historical and prior ledger preserved '+str((a['index'],a['quantity'])),all(b[k]==v for k,v in a.items()))
    for t in b['terms']:
        for v in t['provenance']:provenance[v['path']]=v['sha256']
    if b['quantity']!='P2_average':
        check('missing cannot become completed '+str((b['index'],b['quantity'])),b['proposed_total_utilization'] is None and b['complete_reference_uncertainty'] is None)
for p,h in provenance.items():check('term source '+p,sha(R/p)==h)
for p,h in read(E/'source-manifest.json').items():check('assessment source '+p,sha(pathlib.Path(p))==h)
for n in ('external-manifest.json',):
    m=read(E/n)
    for p,v in m['files'].items():check('external '+p,sha(pathlib.Path(m['root'])/p)==v['sha256'])
check('independent reference replay matches retained review',read(E/'independent-reference-replay.json')==read(R/'ba13-closure-independent-review-20261020/outputs/reference-diagnostics.json'))
check('owner acceptance registry empty',read(P/'acceptance-registry.json')['accepted_conditions']=={} and read(P/'acceptance-registry.json')['accepted_standards']=={})
check('no stellar calculations',read(E/'resources.json')['stellar_evaluations']==0)
check('budget ceiling honored',read(E/'resources.json')['measured_numerical_cpu_seconds']<600)
result=dict(status='PASS',check_count=len(checks),checks=checks,source_commit=git('rev-parse','HEAD'),
            scope='delivery/content verification, not independent scientific acceptance',all_term_sources=len(provenance),changed_paths=paths)
if len(sys.argv)>1:pathlib.Path(sys.argv[1]).write_text(json.dumps(result,indent=1)+'\n')
print(json.dumps(dict(status='PASS',checks=len(checks),term_sources=len(provenance))))
