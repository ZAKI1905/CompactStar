#!/usr/bin/env python3
"""Capture T1/T2 numerical artifacts and optional exact reference comparison."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

EOS_SHA='5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd'
SEQ='NStar/Dark_Core/0.8/0.8_19x19_Sequence.tsv'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def artifacts(root,mode):
    return {str(p.relative_to(root)):sha(p) for p in sorted(root.rglob('*'))
            if p.is_file() and p.suffix in ('.tsv','.eos')
            and p.name!='DS(CMF)-1_with_crust.eos'
            and (mode=='T2' or str(p.relative_to(root)) not in (SEQ,'EOS/Fermi_Gas_0.8mn.eos'))}
def main():
    p=argparse.ArgumentParser()
    p.add_argument('executable',type=Path);p.add_argument('root',type=Path)
    p.add_argument('eos',type=Path);p.add_argument('--mode',choices=['T1','T2'],default='T2')
    p.add_argument('--sequence-root',type=Path);p.add_argument('--reference',type=Path)
    a=p.parse_args();root=a.root.resolve()
    assert len(str(root))<=100 and not root.exists()
    assert sha(a.eos)==EOS_SHA
    (root/'EOS').mkdir(parents=True)
    shutil.copyfile(a.eos,root/'EOS/DS(CMF)-1_with_crust.eos')
    if a.mode=='T1':
        for rel in (SEQ,'EOS/Fermi_Gas_0.8mn.eos'):
            (root/rel).parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(a.sequence_root/rel,root/rel)
    start=time.monotonic()
    with root.with_suffix('.log').open('w') as log:
        result=subprocess.run([str(a.executable),str(root),a.mode],stdout=log,stderr=log,
                              env={**os.environ,'MPLBACKEND':'Agg'})
    output=artifacts(root,a.mode)
    record=dict(mode=a.mode,returncode=result.returncode,wall_seconds=time.monotonic()-start,
                files=output,count=len(output))
    if a.reference:
        reference=artifacts(a.reference,a.mode)
        record['different_files']=[f for f in sorted(output.keys()|reference.keys()) if output.get(f)!=reference.get(f)]
    root.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='files'}))
    assert result.returncode==0
    assert len(output)==(38 if a.mode=='T2' else 13)
    assert not record.get('different_files')
if __name__=='__main__':main()
