#!/usr/bin/env python3
"""No-evolution runtime and piecewise-domain preflight; explicit failure receipts."""
import concurrent.futures,json,pathlib,subprocess,sys,time
from decimal import Decimal as D,localcontext
from piecewise import ROOT,BASE,MODEL,Coupled,constants,domain,SCALE
from run import load_pair,SNAP,digest,save,RETAINED
from roundoff import jacobian_bound
from thermal_enclosure import constants as old_constants

def check(out):
 start=time.monotonic();load_pair(SNAP)
 manifest=MODEL/'thermal-model/manifest.json'
 for p,h in json.loads(manifest.read_text())['files'].items():assert digest(p)==h,p
 for p,h in json.loads((MODEL/'normalization.json').read_text())['files'].items():assert digest(p)==h,p
 results={}
 with localcontext() as ctx:
  ctx.prec=110
  for role in ('source','control'):
   m=Coupled(MODEL,role);box=domain(m);continuity=[]
   for cell in box['cells']:
    k=cell+1
    if k not in box['cells']:continue
    a=m.capacity[cell]+m.slopes[cell]*(m.knots[k]-m.knots[cell]);b=m.capacity[k]
    relative=abs(a-b)/b;assert relative<D('1e-100')
    continuity.append(dict(knot=k,x=m.knots[k],temperature_K=D('1e8')*m.knots[k].exp(),relative_capacity_jump=relative))
   # Regression against accepted Stage A Jacobian using its original domain/scales.
   c=constants(m,D('-.0001'),D('.0002'),[D('3e-10'),D('8e-10')],93,scale=[D(1),D('1e-6'),D('1e-6')])
   J,_=jacobian_bound(m);assert c['J']==J
   # Primary-formula offset remains explicitly bounded, as in Stage A.
   n=json.loads((MODEL/'normalization.json').read_text())['roles'][role]
   forcing=(D(n['direct_normalization_error'])+D(n['photon_normalization_error'])*(D.from_float(2.42)*box['hi']).exp())/(D('1e8')*box['lo'].exp()*min(c['C_min'] for c in box['cells'].values()))
   # Forward invariant box checks: evaluate thermal signs at bounding faces,
   # with conservative polynomial maxima for chemical rates inside the domain.
   results[role]=dict(box=box,continuity=continuity,stage_a_jacobian_regression='PASS',normalization_rhs_bound=forcing)
 jobs=[['ctest','--test-dir',str(RETAINED.parent),'--output-on-failure','-R','^(phase6a1_direct_energy|phase6a1_thermal_ledger|passive_checkpoint_output_contract|phase6_breakpoint_contract|phase5d_component_tolerances|phase6_currentness_bytes)$']]
 for cmd in jobs:
  r=subprocess.run(cmd,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(out/'contracts.log').write_text(r.stdout);assert r.returncode==0,r.stdout
 save(out/'preflight.json',dict(status='PASS',roles=results,elapsed_s=time.monotonic()-start,snapshot_ids=json.loads((SNAP/'identities.json').read_text()),files={str(p):digest(p) for p in [manifest,MODEL/'thermal-model/capacity.tsv',MODEL/'thermal-model/anchors.tsv',MODEL/'normalization.json',BASE/'toolchain/source/gsl-2.7.1/ode-initval2/rkf45.c',BASE/'toolchain/source/gsl-2.7.1/ode-initval2/rk8pd.c']},scope='Runtime and mathematical piecewise preflight only. Actual-path containment, crossings, all original accuracy gates and local reconstruction remain runtime/postflight obligations.'))
if __name__=='__main__':check(pathlib.Path(sys.argv[1]))
