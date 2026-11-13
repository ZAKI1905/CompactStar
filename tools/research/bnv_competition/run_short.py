#!/usr/bin/env python3
"""Bind native S2 invocation to retained owner archive and current source bytes."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
out=Path(sys.argv[1]);assert os.environ.get('SLURM_JOB_ID')
snapshot=ROOT/'docs/validation/results-first-bnv-pilot/parallel-preparation/snapshots/source.json'
s=json.loads(snapshot.read_text())
binary=out/'build/tests/bnv_competition_short'
args=[str(binary),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(out/'short-assembly'),str(out/'short')]
digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources=[p for p in (ROOT/'CompactStar').rglob('*') if p.suffix in ['.cpp','.hpp']]
sources += [ROOT/'tests/native/competition.cpp',ROOT/'tests/bnv/competition_contract.cpp']
receipt=dict(args=args,snapshot_sha256=digest(snapshot),binary_sha256=digest(binary),job=os.environ['SLURM_JOB_ID'],
             source_sha256={str(p.relative_to(ROOT)):digest(p) for p in sources},
             entry_commit=os.environ['COMP_COMMIT'])
(out/'short-invocation.json').write_text(json.dumps(receipt,indent=2)+'\n')
start=time.monotonic()
with (out/'short.log').open('w') as f:r=subprocess.run(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
(out/'short-completion.json').write_text(json.dumps(dict(exit_code=r.returncode,wall_seconds=time.monotonic()-start,job=os.environ['SLURM_JOB_ID']),indent=2)+'\n')
raise SystemExit(r.returncode)
