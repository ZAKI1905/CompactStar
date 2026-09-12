"""Negative controls alter scratch manifests, never a protected source."""
import sys
sys.dont_write_bytecode=True
import copy,json,tempfile
from pathlib import Path
from qualified_suite import ROOT,check_manifest
entry=ROOT/'build/phase5d1c2-audit/entry-hashes.json'
original=check_manifest(entry)
with tempfile.TemporaryDirectory(prefix='manifest-controls-',dir=ROOT/'build') as scratch:
    for name in ['missing_path','changed_hash','duplicate_path','missing_baseline']:
        data=copy.deepcopy(original)
        if name=='missing_path':data['protected'].pop()
        if name=='changed_hash':data['protected'][0]['actual']='0'*64
        if name=='duplicate_path':data['protected'][-1]=data['protected'][0]
        if name=='missing_baseline':data['baselines'].pop()
        p=Path(scratch)/f'{name}.json';p.write_text(json.dumps(data))
        refused=False
        try:check_manifest(p)
        except RuntimeError:refused=True
        if not refused:raise RuntimeError('manifest mutation escaped: '+name)
        print(name,'REFUSED')
check_manifest(entry)
print('PASS exact 33 source paths, 10 baselines, 3 special hashes and scratch negative controls')
