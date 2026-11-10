#!/usr/bin/env python3
"""Collect compact, source-bound report evidence; retain full records externally."""
import csv,json,pathlib,shutil,sys
from decimal import Decimal as D,InvalidOperation
from piecewise import ROOT  # Establish the retained Stage A helper import path.
from run import save,digest
out=pathlib.Path(sys.argv[1]);sealed=pathlib.Path(sys.argv[2]);target=ROOT/'docs/validation/results-first-bnv-pilot/stage-b/evidence'
target.mkdir(exist_ok=True)
def compact(v):
 if isinstance(v,list):return [compact(x) for x in v]
 if isinstance(v,dict):return {k:compact(x) for k,x in v.items()}
 if isinstance(v,str):
  try:
   if any(c in v for c in '.eE') and D(v).is_finite():return format(D(v),'.17g')
  except InvalidOperation:pass
 return v
bindings={}
for name in ('assessment.json','qualification.json','final-summary.json','chemical-response.json','events-assessment.json','ledger-assessment.json','accuracy.json','matched.json','coverage.json','roundoff.json','interpolation-roundoff-completion.json','local-assessment.json','resources.json','finalization-provenance.json','job-status-final.json'):
 p=out/name
 if not p.exists():raise FileNotFoundError(p)
 r=json.loads(p.read_text())
 if name=='assessment.json':
  for role in r['roles'].values():
   for config in role['configurations']:
    config.pop('R20',None)
 if name=='resources.json':
  for build in r['builds']:build.pop('bindings',None)
 save(target/name,compact(r));bindings[name]=dict(source=str(p),sha256=digest(p),report_sha256=digest(target/name),format='17 significant digit display copy; original strings retained externally')
for name in ('temperature.csv','curve-source.csv','curve-control.csv'):
 p=out/name
 with p.open() as f:data=list(csv.DictReader(f))
 with (target/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(compact(data))
 bindings[name]=dict(source=str(p),sha256=digest(p),report_sha256=digest(target/name),format='17 significant digit display copy; original strings retained externally')
for name in ('temperature.png','physics.png','numerics.png','temperature.pdf','physics.pdf','numerics.pdf','test-piecewise.log','test-piecewise-final.log','test-apply-residuals.log','test-apply-residuals-final.log','test-audit-source.log','report-import-after-fix.log'):
 p=out/name;shutil.copy2(p,target/name);bindings[name]=dict(source=str(p),sha256=digest(p))
shutil.copy2(sealed/'receipt.json',target/'manifest-receipt.json')
save(target/'report-bindings.json',dict(status='DISPLAY_COPIES_BOUND_TO_FULL_EXTERNAL_EVIDENCE',files=bindings,receipt_sha256=digest(sealed/'receipt.json')))
