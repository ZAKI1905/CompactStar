#!/usr/bin/env python3
"""Fail closed on retained input or unchanged producer drift; no numerical work."""
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(__file__).resolve().parents[3]
inv=json.loads((root/'docs/validation/results-first-bnv-pilot/native-assembly/evidence/input-inventory.json').read_text())
checks=[]
for x in inv['checks']:
 p=pathlib.Path(x['path']); actual=hashlib.sha256(p.read_bytes()).hexdigest()
 if actual!=x['expected_sha256']: raise SystemExit(f'STOP input mismatch: {p}')
 checks.append({'path':str(p),'sha256':actual})
files=['CompactStar/Analysis/src/ParticleNumberResponse.cpp','CompactStar/Core/src/TOVSolver.cpp',
'CompactStar/Core/src/NStar.cpp','CompactStar/Core/src/RotationSolver.cpp',
'CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp','CompactStar/EOS/src/LocalThermodynamics.cpp',
'tests/analysis/phase5b_freegas_validation.cpp','tests/analysis/particle_number_homogeneous.hpp',
'tests/rotochemical/fixture.hpp','tests/bnv/campaign_fixture.hpp']
for name in files:
 data=(root/name).read_bytes()
 original=subprocess.check_output(['git','show','aa6239b963b7a48946f585af158482659d8fc59f:'+name],cwd=root)
 if data!=original: raise SystemExit(f'STOP producer changed: {name}')
 checks.append({'path':name,'sha256':hashlib.sha256(data).hexdigest()})
print(json.dumps({'status':'PASS','checks':checks},indent=2))
