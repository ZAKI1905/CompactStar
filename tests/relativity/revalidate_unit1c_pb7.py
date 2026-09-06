#!/usr/bin/env python3
"""Test-side PB7 prerequisite: complete stars vs converged homogeneous solution.
No production particle response API. Bounds are ADR-0012 U12's 2e-4, applied
also to each independently estimated numerical uncertainty (subdominance).
"""
from pathlib import Path
import sys,json
sys.dont_write_bytecode = True
import numpy as np
import revalidate_background as b
b.D=Path(sys.argv[1]);out=Path(sys.argv[2]);a=b.load('g8192r40000')
z={'species':['n','p','e','mu'],'status':'IN PROGRESS','sequence':{}}
for g in [4096,8192]:
    counts={j:b.count(b.load(f'seq{g}_{j}')) for j in [-4,-2,-1,1,2,4]}
    B=(8*(counts[1]-counts[-1])-(counts[2]-counts[-2]))/.006
    B2=(8*(counts[2]-counts[-2])-(counts[4]-counts[-4]))/.012
    z['sequence'][str(g)]={'h0005':B.tolist(),'h001':B2.tolist(),'step_relative':(B2/B-1).tolist()}
    assert max(abs(B2/B-1))<2e-4,'PB7 sequence step criterion failure'
out.write_text(json.dumps(z,indent=2))
Hs=[]
for rtol,atol,N in [(2e-10,1e-13,20000),(2e-11,1e-14,40000)]:
    H=b.homogeneous(a,rtol,atol,N);Hs.append(H)
    err=H/B-1;row={'rtol':rtol,'atol':atol,'partitions':N,'H':H.tolist(),'signed_relative':err.tolist(),'absolute_relative':abs(err).tolist()}
    z.setdefault('homogeneous',[]).append(row);out.write_text(json.dumps(z,indent=2));print(row,flush=True)
    assert max(abs(err))<2e-4,'PB7 coherent disagreement >= 2e-4'
floor=abs(Hs[1]/Hs[0]-1);table=abs(np.array(z['sequence']['4096']['h0005'])/B-1);step=abs(B2/B-1)
z['error_estimates']={'homogeneous_refinement_relative':floor.tolist(),'table_relative':table.tolist(),'sequence_step_relative':step.tolist(),'conservative_sum':(floor+table+step).tolist()}
assert max(floor+table+step)<1e-4,'PB7 error estimate not subdominant to 2e-4'
z['status']='ADR-0011 Q4 UNIT PREREQUISITE NUMERICALLY SATISFIED IN UNIT-1 REVALIDATION CANDIDATE — PENDING INDEPENDENT REVIEW AND HUMAN RATIFICATION'
z['boundary']='PB7 itself incomplete; ADR-0011 Q4 remains pending; no Phase-5B production APIs.'
out.write_text(json.dumps(z,indent=2));print(z['status'],flush=True)
