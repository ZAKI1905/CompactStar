#!/usr/bin/env python3
"""Merge only complete, hash-verified, isolated BA13 coefficient shards."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
from ba13_stage import verify

p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
if a.output.exists() or a.output.with_suffix('.partial').exists():raise RuntimeError('fresh collection root required')
for i in range(21):verify(a.run/f'coefficient-{i}')
verify(a.run/'stars');verify(a.run/'native-adapters')
for i in range(21):
    shard=a.run/f'coefficient-{i}'
    for name in ('target_stars.tsv','tangent.tsv'):
        if (shard/name).read_bytes() != (a.run/'stars'/name).read_bytes():
            raise RuntimeError('shard target/tangent identity mismatch')
    original=a.run/'native-adapters'/'chemical-profiles'/f'target-{i}'
    for p in original.rglob('*'):
        if p.is_file() and p.read_bytes() != (shard/'chemical-profiles'/f'target-{i}'/p.relative_to(original)).read_bytes():
            raise RuntimeError('shard profile/adapter identity mismatch')
out=a.output.with_suffix('.partial');out.mkdir(parents=True)
for name in ('target_stars.tsv','tangent.tsv'):shutil.copyfile(a.run/'stars'/name,out/name)
for stem,count in (('coefficients',1),('thermal_surface',13),('profile_metrics',2049)):
    header=None
    with (out/(stem+'.tsv')).open('w') as f:
        for i in range(21):
            path=a.run/f'coefficient-{i}'/f'{stem}-shard-{i}.tsv'
            lines=path.read_text().splitlines(keepends=True)
            if len(lines)!=count+1:raise RuntimeError('incomplete shard')
            if header is None:header=lines[0];f.write(header)
            if lines[0]!=header:raise RuntimeError('header mismatch')
            rows=list(csv.DictReader(lines,delimiter='\t'))
            if any(int(row['index'])!=i for row in rows):raise RuntimeError('shard identity mismatch')
            f.writelines(lines[1:])
shutil.copytree(a.run/'native-adapters'/'chemical-profiles',out/'chemical-profiles')
files={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}
(out/'COMPLETE.json').write_text(json.dumps(dict(files=files,stage='collected static raw samples'),indent=2)+'\n')
for p in out.rglob('*'):
    if p.is_file():
        with p.open('rb') as f:os.fsync(f.fileno())
fd=os.open(out,os.O_RDONLY);os.fsync(fd);os.close(fd)
out.rename(a.output)
fd=os.open(a.output.parent,os.O_RDONLY);os.fsync(fd);os.close(fd)
print('COMPLETE',a.output)
