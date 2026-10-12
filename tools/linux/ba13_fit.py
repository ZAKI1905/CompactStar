#!/usr/bin/env python3
"""Run the byte-authenticated historical static fitter with atomic publication.

Publication denotes a complete candidate, never scientific qualification.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from ba13_stage import digest, verify, rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--historical-package',type=Path)
    a=p.parse_args()
    source=Path(__file__).resolve().parents[2]
    fitter=source/'tests/bnv/compare_candidate.py'
    if digest(fitter)!='9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370':
        raise RuntimeError('historical fitter bytes changed')
    if a.historical_package:
        root=a.historical_package.resolve()
        manifest={line.split('  ',1)[1]:line.split('  ',1)[0] for line in (root/'SHA256SUMS').read_text().splitlines()}
        for name,sha in manifest.items():
            if digest(root/name)!=sha:raise RuntimeError('historical archive input changed: '+name)
        if a.work.resolve()!=root/'original/phase6a1-frozen-certificate-final':
            raise RuntimeError('historical input identity mismatch')
    else:
        verify(a.work)
    inputs={str(f.relative_to(a.work)):digest(f) for f in a.work.rglob('*') if f.is_file()}
    for name,count in [('target_stars',21),('tangent',21),('coefficients',21),('thermal_surface',273),('profile_metrics',43029)]:
        rows(a.work/(name+'.tsv'),count)
    partial=a.output.with_suffix('.partial')
    if a.output.exists() or partial.exists():raise RuntimeError('fresh fit output required')
    partial.mkdir(parents=True)
    cmd=[sys.executable,'-B',str(fitter),'frozen','--work',str(a.work.resolve()),
         '--output-json',str(partial/'certificate.json'),'--output-tsv',str(partial/'certificate.tsv')]
    with (partial/'fit.log').open('x') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    certificate=json.loads((partial/'certificate.json').read_text())
    if len(certificate['samples'])!=21 or len(certificate['fits'])!=22:
        raise RuntimeError('incomplete fitted certificate')
    rows_tsv=list(__import__('csv').DictReader((partial/'certificate.tsv').open(),delimiter='\t'))
    if len(rows_tsv)!=462:raise RuntimeError('incomplete serialized certificate')
    for name,sha in inputs.items():
        if digest(a.work/name)!=sha:raise RuntimeError('input mutated during fit')
    receipt=dict(status='HISTORICAL_SEMANTICS_CANDIDATE_NOT_PROMOTED',command=cmd,
        fitter_sha256=digest(fitter),input_root=str(a.work.resolve()),inputs=inputs,
        slurm_job=os.environ.get('SLURM_JOB_ID'),node=os.uname().nodename,
        limitations=['Historical zero/absent error assignments are preserved, not independently validated',
                     'A passing original fitter does not establish platform agreement or owner acceptance'])
    (partial/'binding.json').write_text(json.dumps(receipt,indent=2)+'\n')
    files={str(f.relative_to(partial)):digest(f) for f in partial.rglob('*') if f.is_file()}
    (partial/'COMPLETE.json').write_text(json.dumps(dict(stage='historical frozen fit candidate',files=files),indent=2)+'\n')
    for f in partial.iterdir():
        with f.open('rb') as stream:os.fsync(stream.fileno())
    fd=os.open(partial,os.O_RDONLY);os.fsync(fd);os.close(fd)
    partial.rename(a.output)
    fd=os.open(a.output.parent,os.O_RDONLY);os.fsync(fd);os.close(fd)
    print((a.output/'fit.log').read_text(),end='')
    print('COMPLETE candidate',a.output)

if __name__=='__main__':main()
