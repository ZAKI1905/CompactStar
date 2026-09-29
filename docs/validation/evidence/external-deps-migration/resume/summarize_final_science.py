from pathlib import Path
import csv,hashlib,json
q=Path(__file__).resolve().parent
hist=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6-adr0017-production-implementation/build/phase6-adr0017-production-qualification/result.json')
roots={'old':q/'adr0017-old-Debug','new':q/'adr0017-new-Debug-attempt2'}
results={k:json.loads((r/'result.json').read_text()) for k,r in roots.items()}
historical=json.loads(hist.read_text())
science=lambda x:{k:v for k,v in x.items() if k!='performance'}
assert all(x['pass'] for x in results.values())
assert science(results['old'])==science(results['new'])==science(historical)
non_timing=lambda x:{k:v for k,v in x['performance'].items() if not k.endswith('_s')}
assert non_timing(results['old'])==non_timing(results['new'])==non_timing(historical)
files={};checkpoints={};ignored=[]
for key,root in roots.items():
 files[key]={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (root/'output').iterdir() if f.is_file() and f.name not in ['performance.tsv','checkpoints.tsv']}
 with (root/'output/checkpoints.tsv').open() as f:
  reader=csv.DictReader(f,delimiter='\t');ignored=[k for k in reader.fieldnames if k.endswith(('_wall_s','_cpu_s'))]
  checkpoints[key]=[{k:v for k,v in row.items() if k not in ignored} for row in reader]
assert files['old']==files['new']
assert checkpoints['old']==checkpoints['new']
assert len(checkpoints['old'])==241
summary=dict(result='PASS',scientific_result_equal_to_historical=True,main_artifact_hashes=files,checkpoint_records=241,checkpoint_scientific_record_sha256=hashlib.sha256(json.dumps(checkpoints['old'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),excluded_checkpoint_columns=ignored,scientific_result=science(results['new']),non_timing_performance_fields=non_timing(results['new']))
(q/'evidence/ADR0017-same-mode-comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
rows=[]
for stack in ['old','new']:
 for mode in ['Debug','Release']:
  receipts=sorted((q/f'{stack}-{mode}').rglob('regression-evidence.json'))
  assert receipts,(stack,mode)
  for receipt in receipts:
   data=json.loads(receipt.read_text());side=json.loads(Path(data['execution_sidecar']).read_text())
   assert data['producer_raw_rc']==0
   rows.append(dict(stack=stack,mode=mode,receipt=str(receipt),receipt_sha256=hashlib.sha256(receipt.read_bytes()).hexdigest(),generated_sha256=data['generated_sha256'],baseline_sha256=data['baseline_sha256'],controls=data['controls'],sidecar=side))
(q/'evidence/Phase5D1-complete-matrix.json').write_text(json.dumps(rows,indent=2)+'\n')
print('ADR0017 OLD/NEW/historical science exact; four complete Phase5D1 same-baseline receipts authenticated')
