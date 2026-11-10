#!/usr/bin/env python3
"""Independent coupled reference over the authorized horizon; no native replacement."""
import argparse,json,os,pathlib,time
from decimal import Decimal as D,localcontext
from piecewise import MODEL,Coupled,integrate,END,ROOT,BASE
from run import save,digest
p=argparse.ArgumentParser();p.add_argument('out',type=pathlib.Path);p.add_argument('role',choices=['source','control']);p.add_argument('--step',type=int,default=1000000);a=p.parse_args()
assert json.loads((a.out/'preflight.json').read_text())['status']=='PASS'
start=time.monotonic()
with localcontext() as ctx:
 ctx.prec=110;m=Coupled(MODEL,a.role)
 result=integrate(m,[END*i/1920 for i in range(1921)],D(a.step))
 result.update(role=a.role,precision=110,max_step_s=a.step,elapsed_s=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID'),status='CANDIDATE_PIECEWISE_REFERENCE',inputs={str(p):digest(p) for p in [MODEL/'thermal-model/capacity.tsv',MODEL/'thermal-model/anchors.tsv',ROOT/'tools/linux/stage_b/piecewise.py',pathlib.Path(__file__)]})
 save(a.out/f'hp-{a.role}-{a.step}.json',result)
