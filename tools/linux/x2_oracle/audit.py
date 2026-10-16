#!/usr/bin/env python3
"""Independent Decimal/provenance checks; no ODE solve, no score changes."""
import decimal,hashlib,json,os,pathlib,resource,subprocess
D=decimal.Decimal;decimal.getcontext().prec=80
SRC=pathlib.Path(__file__).resolve().parents[3];DOC=SRC/'docs/validation/x2-fixed-step';OUT=SRC.parents[1]/'x2-fixed-step'
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});assert os.getenv('SLURM_JOB_ID')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(),parse_float=D)
counts={};
for gate in ['execution-gate.json','x1-diagnostic-gate.json','continuation-gate.json','shell-gate.json']:
 x=json.loads((DOC/gate).read_text())
 for p,h in x['files'].items():assert sha(pathlib.Path(p))==h,p
 counts[gate]=len(x['files'])
x1=SRC/'docs/validation/x1-oracle-anchor'
for manifest,base in [(x1/'SHA256SUMS',x1),(x1/'frozen/SHA256SUMS',x1/'frozen')]:
 for line in manifest.read_text().splitlines():
  h,f=line.split(None,1);assert sha(base/f)==h,f
x=json.loads((x1/'execution-gate.json').read_text())
for group in ['files','review_files']:
 for p,h in x[group].items():assert sha(pathlib.Path(p))==h,p
for doc in [x1/'evidence/external-manifest.json',DOC/'evidence/external-manifest.json']:
 x=json.loads(doc.read_text())
 for p,r in x['files'].items():
  f=pathlib.Path(x['root'])/p;assert sha(f)==r['sha256'] and f.stat().st_size==r['bytes'],str(f)
 counts[str(doc.relative_to(SRC))]=len(x['files'])
old=json.loads((x1/'evidence/runs.json').read_text());retained=json.loads((x1/'frozen/oracle-retained.json').read_text());keys=['B','R','M','r_e','L','C','rb','mu_inf','Hc']
reproduced=0
for v in ['regularized','unregularized']:
 for tol in ['2e-10','2e-11','2e-12','2e-13']:
  rows=[r['values'] for r in old if r['target']=='0' and r['variant']==v and r['tol']==tol]
  if rows:ref=[rows[0][k] for k in keys]
  else:ref=[retained['0:'+tol][{'r_e':'r_nu','L':'L_nu','C':'C_nu'}.get(k,k)] for k in keys]
  actual=list(map(float,(OUT/f'x1-{v}-baseline-{tol}.out').read_text().split()))
  assert actual==ref,(v,tol,actual,ref);reproduced+=1
counts['bit_identical_X1_baseline_replays']=reproduced
scores=json.loads((DOC/'evidence/scores.json').read_text());order_pass=0;targets=read(x1/'targets.json')
for t in ['0','10','20']:
 vals=[read(OUT/f't{t}-n{n}.json')['values'] for n in [32,64,128,256,512,1024,2048,4096]]
 for row in vals:assert abs(row[0]/D(targets[t]['B'])-1)<=D('2e-16')
 for i,k in enumerate(keys[1:],1):
  diffs=[vals[j+1][i]-vals[j][i] for j in range(7)]
  orders=[(abs(diffs[j]/diffs[j+1])).ln()/D(2).ln() for j in range(6)]
  assert all(D('3.5')<=o<=D('4.5') for o in orders[-2:])
  assert diffs[-3]*diffs[-2]>0 and diffs[-2]*diffs[-1]>0
  assert scores[t][k]['status']=='PASS'
  estimate=abs(diffs[-1])/15/abs(vals[-1][i]);reported=D(str(scores[t][k]['conditional_Richardson_relative']))
  assert abs(estimate/reported-1)<D('1e-14');order_pass+=1
  assert abs(orders[-1]-D(str(scores[t][k]['orders'][-1])))<D('1e-12')
 final=vals[-1];assert abs(final[10])<D('1e-7') and abs(final[11])<D('1e-7') and abs(final[12])<D('1e-12')
 assert scores[t]['B']['status']=='INCONCLUSIVE'
counts['physical_scalar_and_Hc_order_PASS']=order_pass
# Source support correction must eliminate lower-shell JL exactly in raw decimal text.
for p in OUT.glob('shell-*.err'):
 rows=[[D(x) for x in l.split()[1:]] for l in p.read_text().splitlines() if l.startswith('SEG')][-4:]
 assert rows[-1][5]==rows[-2][5],p.name
counts['zero_sub_neutron_JL_tests']=len(list(OUT.glob('shell-*.err')))
headroom=json.loads((DOC/'evidence/headroom.json').read_text());hist=read(SRC/'docs/validation/ba13-recovery/evidence/fits/run1/certificate-summary.json')
for row in headroom['rows']:
 s=hist['samples'][row['target']];q=row['quantity'];expected=s['threshold'][q]-s['drift_bound'][q]
 assert abs(D(str(row['headroom']))-expected)<D('1e-19')
counts['historical_headrooms']=len(headroom['rows'])
changed=subprocess.check_output(['git','diff','--name-only','6cc7e20','HEAD'],cwd=SRC,text=True).splitlines()
assert all(p.startswith(('docs/validation/x2-fixed-step/','tools/linux/x2_oracle/')) or p=='tests/bnv/x2_fixed_step.cpp' for p in changed)
prior=sum(json.loads(p.read_text())['cpu_seconds'] for p in OUT.glob('receipt-*.json'))
r=resource.getrusage(resource.RUSAGE_SELF);cpu=r.ru_utime+r.ru_stime
result={'status':'PASS','purpose':'Independent arithmetic, raw replay and authentication audit; no global numerical enclosure','counts':counts,'cpu_seconds':cpu,'prior_numerical_cpu_seconds':prior,'combined_cpu_seconds':prior+cpu,'affinity':sorted(os.sched_getaffinity(0)),'job':os.getenv('SLURM_JOB_ID'),'source_sha256':sha(pathlib.Path(__file__)),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=SRC,text=True).strip(),'changed_files':changed}
assert prior+cpu<900
(OUT/'independent-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
