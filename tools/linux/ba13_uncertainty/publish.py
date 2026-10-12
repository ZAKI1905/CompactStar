#!/usr/bin/env python3
"""Fail-closed publication of completed pilot evidence, never a certificate."""
import csv,hashlib,json,os,subprocess
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=Path(os.environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
required=['summary.json','disposition-matrix.json','disposition-matrix.tsv','resources.json','measurements.tsv','paired-integrals.json','holdout-validation.json','holdout-predictions.json','target-holdout-predictions.json','thermal-diagnostic.json','fit-units.json','chemical-diagnostic.json','integral-diagnostic.json','profile-diagnostic.json','p1-validation.json','eos-ordering.tsv']
for n in required:assert (R/n).is_file(),n
for n,count in [('integral-diagnostic.json',27),('profile-diagnostic.json',55),('p1-validation.json',27),('holdout-validation.json',352)]:
 assert len(json.loads((R/n).read_text()))==count,n
t=json.loads((R/'thermal-diagnostic.json').read_text());assert len(t['records'])==546 and len(t['high_precision'])==27
c=json.loads((R/'chemical-diagnostic.json').read_text())['records'];assert [v['index'] for v in c]==list(range(21))
assert all(len(v['integrals'])==5 for v in c if v['index'] in (0,10,20))
u=json.loads((R/'fit-units.json').read_text())['quantities'];assert len(u)==22 and sum(v['components'] for v in u.values())==10303
matrix=json.loads((R/'disposition-matrix.json').read_text());assert len(matrix)==462
assert len({(v['index'],v['quantity']) for v in matrix})==462
assert all(v['proposed_qualified_uncertainty'] is None or v['quantity']=='P2_average' for v in matrix)
assert all(v['qualified_utilization'] is None or v['quantity']=='P2_average' for v in matrix)
summary=json.loads((R/'summary.json').read_text());assert not summary['qualified']
complete=[];failed=[]
for p in sorted((R/'runs').glob('*.receipt.json')):
 v=json.loads(p.read_text());dest=R/'runs'/p.name.replace('.receipt.json','')
 if v['status']=='REFUSED':assert not (dest/'COMPLETE.json').exists();failed.append(dest.name);continue
 m=json.loads((dest/'COMPLETE.json').read_text())
 for n,h in m['files'].items():assert sha(dest/n)==h,(dest,n)
 assert set(m['files'])=={str(q.relative_to(dest)) for q in dest.rglob('*') if q.is_file() and q.name!='COMPLETE.json'}
 complete.append(dest.name)
assert len(complete)==55 and len(failed)==3
raw=S.parents[1]/'ba13-recovery/run1/raw'
def table(p):return list(csv.DictReader(p.open(),delimiter='\t'))
old=table(raw/'target_stars.tsv');heat=table(raw/'thermal_surface.tsv')
nominal_checks=[]
for i in (0,5,10,15,20):
 p=R/'runs'/f'nominal-{i}';st=table(p/'star.tsv')[0];t=table(p/'thermal.tsv')
 check=dict(index=i,B_exact=float(st['B_solved'])==float(old[i]['B_solved']),Cstar_exact=all(float(a['Cstar'])==float(b['Cstar_erg_K']) for a,b in zip(t,heat[i*13:(i+1)*13])))
 check['B_target_exact']=float(st['B_target'])==float(old[i]['B_target'])
 check['achieved_B_change_relative']=(float(st['B_solved'])-float(old[i]['B_solved']))/float(old[i]['B_solved'])
 assert check['B_target_exact']
 if i in (0,5,10,20):assert check['B_exact'] and check['Cstar_exact'],check
 check['binding_status']='EXACT_ORIGINAL' if check['B_exact'] and check['Cstar_exact'] else 'DISTINCT_COMMON_TARGET_RESOLVE_NO_CERTIFICATE_BINDING'
 nominal_checks.append(check)
files={}
for p in R.rglob('*'):
 if not p.is_file() or p.name in ('STUDY_COMPLETE.json','STUDY_COMPLETE.partial.json'):continue
 rel=p.relative_to(R)
 if rel.parts[0]=='build' and str(rel) not in ('build/dependency-provenance.json','build/CMakeCache.txt','build/tests/ba13_uncertainty_pilot','build/tests/ba13_uncertainty_eos_audit','build/tests/ba13_original_helpers.hpp'):continue
 files[str(rel)]=dict(sha256=sha(p),size=p.stat().st_size)
sources=list((S/'tools/linux/ba13_uncertainty').glob('*'))+[S/'tests/bnv/uncertainty_pilot.cpp',S/'tests/bnv/uncertainty_eos_audit.cpp',S/'tests/CMakeLists.txt',S/'docs/validation/ba13-uncertainty/PREDECLARATION.md']
record=dict(status='COMPLETE_BOUNDED_PILOT_NOT_QUALIFIED',root=str(R),source_commit=subprocess.check_output(['git','-C',str(S),'rev-parse','HEAD'],text=True).strip(),
 source_files={str(p.relative_to(S)):sha(p) for p in sources if p.is_file()},completed_cases=complete,refused_cases=failed,nominal_identity_checks=nominal_checks,files=files)
p=R/'STUDY_COMPLETE.partial.json';p.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
with p.open('rb') as f:os.fsync(f.fileno())
p.replace(R/'STUDY_COMPLETE.json');print('verified publication',len(files),'files',len(complete),'complete',len(failed),'refused')
