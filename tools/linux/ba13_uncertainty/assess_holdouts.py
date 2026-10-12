#!/usr/bin/env python3
import csv,json
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')))
r=json.loads((R/'holdout-predictions.json').read_text())+json.loads((R/'target-holdout-predictions.json').read_text())
for v in r:
 case={'radial':'r320000','eos':'e32768'}[v['axis']]
 path=R/'runs'/f'{case}-{v["index"]}'/v['table']
 if not path.exists():
  v['heldout_status']='INCONCLUSIVE_EOS_GENERATION_REFUSED';v['heldout_value']=None;continue
 rows=list(csv.DictReader(path.open(),delimiter='\t'))
 y=float(rows[v['knot']][v['quantity']]);v['heldout_value']=y;v['heldout_increment']=abs(y-v['values'][-1])
 v['heldout_status']=('PASS_CONDITIONAL_PREDICTION' if v['heldout_increment']<=v['finest_remainder_prediction'] else 'FAIL_PREDICTION') if v['estimator_eligible'] else 'INCONCLUSIVE'
(R/'holdout-validation.json').write_text(json.dumps(r,indent=2)+'\n')
from collections import Counter
print(Counter(v['heldout_status'] for v in r))
for q in sorted(set(v['quantity'] for v in r)):
 print(q,Counter(v['heldout_status'] for v in r if v['quantity']==q))
