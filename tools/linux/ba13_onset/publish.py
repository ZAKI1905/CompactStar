"""Authenticate retained pilot evidence and copy only compact records into Git."""
import gzip,hashlib,json,subprocess
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=S.parents[1];OUT=R/'ba13-onset-pilot';DOC=S/'docs/validation/ba13-onset-pilot/evidence'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
# Fail closed on protected owner/fitter/certificate drift.
subprocess.run(['python3','-B',str(Path(__file__).parent/'authenticate.py')],check=True,cwd=S)
for name,commit in [('PREDECLARATION.md','8bb72e2'),('evidence/predictions.json','e657efb'),('evidence/eos-predictions.json','650871e')]:
 current=S/'docs/validation/ba13-onset-pilot'/name
 original=subprocess.check_output(['git','show',commit+':docs/validation/ba13-onset-pilot/'+name],cwd=S)
 assert original==current.read_bytes(), 'frozen file changed: '+name
verified=[]
for p in (OUT/'runs').glob('*/COMPLETE.json'):
 for n,h in json.loads(p.read_text())['files'].items():assert sha(p.parent/n)==h,(p,n)
 verified.append(p.parent.name)
resources=json.loads((OUT/'resources.json').read_text());score=json.loads((OUT/'score-resource.json').read_text());total=resources['total_cpu']+score['cpu_seconds']
assert total<=3600, 'numerical resource ceiling exceeded'
# Build is reported separately, including all build wrapper records.
build=[]
for p in OUT.glob('build-*.log'):
 for line in p.read_text().splitlines():
  if line.startswith('RESOURCE_USAGE_JSON '):build.append(json.loads(line[len('RESOURCE_USAGE_JSON '):]))
record=dict(numerical_cpu_seconds_including_conservative_startup_charge=total,build_cpu_seconds=sum(v['user_cpu_seconds']+v['system_cpu_seconds'] for v in build),ceiling_seconds=3600,completed_new_stars=len(verified),numerical_job_ids=['14687','14689','14690','14691'],failed_setup_job='14686 (no numerical work)',note='Publication/authentication hashing is I/O bookkeeping, not numerical work; no Slurm accounting database. Process receipts and conservative sequential wall guard retained.')
(OUT/'resource-summary.json').write_text(json.dumps(record,indent=2)+'\n')
files={str(p.relative_to(OUT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='PILOT_COMPLETE.json'}
source={str(p.relative_to(S)):sha(p) for folder in [S/'tools/linux/ba13_onset',S/'docs/validation/ba13-onset-pilot'] for p in sorted(folder.rglob('*')) if p.is_file() and 'evidence' not in p.parts}
source.update({str(p.relative_to(S)):sha(p) for p in (S/'tests/bnv').glob('onset*')})
protected={name:sha(S/name) for name in ['tests/eos/structure1/table.hpp','tests/bnv/compare_candidate.py','tests/bnv/frozen_sensitivity.cpp','CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp','CompactStar/Core/src/TOVSolver.cpp','CompactStar/Core/src/NStar.cpp','CompactStar/Physics/Evolution/src/StarContext.cpp','CompactStar/Physics/Driver/Thermal/Boundary/src/TbDefinition.cpp']}
manifest=dict(status='VALIDATION_ONLY_NOT_QUALIFIED',entry='4f3cf97462ca030f715e26808938f781a2a40945',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=S,text=True).strip(),predeclaration_commit='8bb72e2',radial_prediction_commit='e657efb',final_blind_gate_commit='650871e',root=str(OUT),files=files,source_files=source,protected=protected,completed_cases=verified,resources=record)
(OUT/'PILOT_COMPLETE.json').write_text(json.dumps(manifest,indent=2)+'\n')
for name in ['scores.json','high-precision.json','resource-summary.json','constants.json']:(DOC/name).write_bytes((OUT/name).read_bytes())
for name in ['metrics.json','discovery.json','resources.json','PILOT_COMPLETE.json']:
 with (DOC/(name+'.gz')).open('wb') as stream:
  with gzip.GzipFile(filename='',mode='wb',fileobj=stream,mtime=0) as g:g.write((OUT/name).read_bytes())
receipt=dict(manifest=str(OUT/'PILOT_COMPLETE.json'),sha256=sha(OUT/'PILOT_COMPLETE.json'),files=len(files),retained_bytes=sum(v['bytes'] for v in files.values()),resources=record)
(DOC/'artifact-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
