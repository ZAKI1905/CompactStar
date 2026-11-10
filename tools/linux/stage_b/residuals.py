#!/usr/bin/env python3
"""Extend Stage A actual-stage residual propagation to the Stage B piecewise RHS.

The global positive Jacobian is valid across continuous Cstar knots. Native
smooth-extension discrepancies enter measured RHS defects explicitly. This is
a conservative arithmetic/method-evaluation allowance, not a truncation bound.
"""
import concurrent.futures,pathlib,sys
from piecewise import MODEL,domain,Coupled
import importlib.util
spec=importlib.util.spec_from_file_location('stage_a_roundoff',pathlib.Path(__file__).parents[1]/'preliminary_pair/roundoff.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
from run import save

def jacobian(model):
 b=domain(model)
 J=[[max(c['J'][i][j] for c in b['cells'].values()) for j in range(3)] for i in range(3)]
 return J,dict(x_min=b['lo'],x_max=b['hi'],eta_abs=b['Y'],cells=list(b['cells']),piecewise_continuous=True)

def one(task):
 out,role,k=task;old.jacobian_bound=jacobian
 r=old.audit(MODEL,role,f'reference-{k}',directory_override=out/f'{role}-supplement'/f'reference-{k}',write=False)
 r['assumptions'][0]='Frozen binary64 Cstar table with piecewise interpolation and evaluated photon/direct anchors defines the real RHS; normalization uncertainty remains separate.'
 r['assumptions'].append('The measured RHS residual includes smooth native cell-extension versus real piecewise-field differences at event strips; no single-cell derivative assumption is made.')
 target=out/'roundoff';target.mkdir(exist_ok=True);save(target/f'{role}-reference-{k}.json',r)
 return dict(role=role,reference=k,status=r['status'],final_enclosure=r['final_enclosure'],maximum_rhs_defect=r['maximum_rhs_defect'],maximum_replay_difference=r['maximum_replay_difference'])
if __name__=='__main__':
 out=pathlib.Path(sys.argv[1]);workers=int(sys.argv[2]);tasks=[(out,r,k) for r in ('source','control') for k in range(5)]
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:results=list(pool.map(one,tasks))
 save(out/'roundoff.json',dict(status='CONDITIONAL_PIECEWISE_RESIDUAL_ENCLOSURE',histories=results))
