#!/usr/bin/env python3
"""Decimal analysis of frozen raw outputs. No stellar solves or X1 rescoring."""
import decimal,hashlib,json,math,pathlib
D=decimal.Decimal;decimal.getcontext().prec=60
SRC=pathlib.Path(__file__).resolve().parents[3];OUT=SRC.parents[1]/'x2-fixed-step';DOC=SRC/'docs/validation/x2-fixed-step';E=DOC/'evidence';E.mkdir(exist_ok=True)
keys=['B','R','M','r_e','L','C','rb','mu_inf','Hc','mass_stage_residual','mass_audit_residual','B_audit_residual','mu_drift','stage_residual','max_u_defect']
ns=[32,64,128,256,512,1024,2048,4096]
def read(p):return json.loads(p.read_text(),parse_float=D)
def convert(x):
 if isinstance(x,D):return float(x)
 if isinstance(x,dict):return {k:convert(v) for k,v in x.items()}
 if isinstance(x,list):return [convert(v) for v in x]
 return x
def save(name,x):(E/name).write_text(json.dumps(convert(x),indent=2,allow_nan=False)+'\n')
def sequence(rows,i):
 ds=[rows[j+1][i]-rows[j][i] for j in range(len(rows)-1)]
 ps=[math.log2(abs(float(ds[j]/ds[j+1]))) if ds[j] and ds[j+1] else None for j in range(len(ds)-1)]
 # This arithmetic-resolution screen is diagnostic only; no X1 allowance.
 floor=D(64)*D(2)**-63*abs(rows[-1][i])
 gate=len(ds)>=3 and all(p is not None and 3.5<=p<=4.5 for p in ps[-2:]) and ds[-3]*ds[-2]>0 and ds[-2]*ds[-1]>0 and min(abs(x) for x in ds[-3:])>floor
 return dict(increments=ds,orders=ps,status='PASS' if gate else 'INCONCLUSIVE',arithmetic_screen_absolute=floor,conditional_Richardson_relative=abs(ds[-1])/15/abs(rows[-1][i]) if gate else None)
x1=read(SRC/'docs/validation/x1-oracle-anchor/evidence/runs.json');anchor=read(SRC/'docs/validation/x1-oracle-anchor/frozen/anchor.json');targets=read(SRC/'docs/validation/x1-oracle-anchor/targets.json')
allrows={};scores={};raw={}
for t in ['0','10','20']:
 rows=[]
 for n in ns:
  p=OUT/f't{t}-n{n}.json';r=read(p);raw[p.name]=r;rows.append(r['values'])
 allrows[t]=rows;scores[t]={k:sequence(rows,i) for i,k in enumerate(keys[:9])}
 # B is a shooting residual, not an order claim, even if it accidentally contracts.
 scores[t]['B']['status']='INCONCLUSIVE';scores[t]['B']['conditional_Richardson_relative']=None
 scores[t]['alignment_max']=max(abs(r[0]/D(targets[t]['B'])-1) for r in rows)
 scores[t]['anchor']={q:dict(residual=rows[-1][keys.index(q)]/r['predicted_analytic_limit']-1,original_band=r['band_total_excluding_oracle'],classification='DESCRIPTIVE ONLY; X1 A2 stays INCONCLUSIVE') for q,r in anchor['targets'][t].items()}
 scores[t]['cross_method']={v:{q:rows[-1][keys.index(q)]/next(r['values'][q] for r in x1 if r['target']==t and r['tol']=='2e-13' and r['variant']==v)-1 for q in keys[:9]} for v in ['regularized','unregularized']}
 scores[t]['cross_method_status']='INCONCLUSIVE: no eligible global X1 uncertainty; residuals reported without an invented allowance'
 scores[t]['paired']={q:rows[-1][i]/allrows['0'][-1][i]-1 for i,q in enumerate(keys[:9])}
 scores[t]['identities']=dict(zip(keys[9:],rows[-1][9:]))
ablations={}
for tag in ['center-1e-2','center-1e-6','eos192','plain','raw','double-newton']:
 rows=[];diff={}
 for n in [512,1024,2048]:
  p=OUT/f'{tag}-n{n}.json'
  if not p.read_text().strip():diff[str(n)]={'status':'FAIL','error':p.with_suffix('.err').read_text().strip()};continue
  r=read(p);raw[p.name]=r;rows.append(r['values']);b=read(OUT/f't0-n{n}.json')['values'];diff[str(n)]={q:r['values'][i]/b[i]-1 for i,q in enumerate(keys[:9])}
 ablations[tag]={'differences':diff,'C_sequence':sequence(rows,5) if len(rows)==3 else None}
for v in ['regularized','unregularized']:
 p=OUT/f'fixed-x1-{v}.json';r=read(p);raw[p.name]=r
 old=next(r['values'] for r in x1 if r['target']=='0' and r['variant']==v and r['tol']=='2e-13')
 ablations['same_H_'+v]={q:old[q]/r['values'][i]-1 for i,q in enumerate(keys[:9])}

def trace(p):
 text=p.read_text();rows=[[D(x) for x in l.split()[1:]] for l in text.splitlines() if l.startswith('SEG')]
 return rows[-4:]
diagnostics={}
for p in sorted(list(OUT.glob('x1-*.out'))+list(OUT.glob('shell-*.out'))):
 vals=[D(x) for x in p.read_text().split()];seg=trace(p.with_suffix('.err'));last=seg[-1]
 diagnostics[p.stem]=dict(values=dict(zip(keys[:9],vals)),segments=seg,mass_closure=last[7]/last[3]-1,sub_neutron_JL=seg[-1][5]-seg[-2][5],sub_neutron_JC=seg[-1][6]-seg[-2][6])
# Preserve exact decimal raw numeric text in compact evidence, not only binary64 JSON.
raw_text={p.name:p.read_text() for p in sorted(OUT.glob('*.json')) if p.name.startswith(('t0-','t10-','t20-','fixed-x1-','double','center','eos192','plain','raw','checks'))}
save('runs.json',raw);save('raw-decimal-text.json',raw_text);save('scores.json',scores);save('ablations.json',ablations);save('x1-diagnostics.json',diagnostics)
receipts=[json.loads(p.read_text()) for p in sorted(OUT.glob('receipt-*.json'))];save('resources.json',{'cpu_seconds':sum(r['cpu_seconds'] for r in receipts),'ceiling':900,'receipts':receipts})
checks=read(OUT/'checks.json');save('checks.json',checks)
# Authentication is rerun at delivery, with contents bound independently of timestamps.
manifest={'root':str(OUT),'source_commit_at_analysis':__import__('subprocess').check_output(['git','rev-parse','HEAD'],cwd=SRC,text=True).strip(),'files':{str(p.relative_to(OUT)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(OUT.rglob('*')) if p.is_file()}}
save('external-manifest.json',manifest)
print('CPU',sum(r['cpu_seconds'] for r in receipts))
for t in scores:
 print(t,[(q,scores[t][q]['status'],scores[t][q]['orders'][-1],scores[t][q]['conditional_Richardson_relative']) for q in ['M','R','r_e','C','L']])
 print('anchor',scores[t]['anchor'])
