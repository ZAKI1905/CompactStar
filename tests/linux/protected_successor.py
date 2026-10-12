#!/usr/bin/env python3
"""Exact opt-in candidate authority; never overrides the historical checker."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tests/rotochemical'))
import manifest
AUTHORITY=ROOT/'docs/validation/linux-corrections/protected-successor.json'
AUTHORITY_SHA256='bb7c0c7e3882cb9b9e16047a17a2bc8c6939541704006ef416f1e727d7d1cbd0'


def check(root, record, authority):
    if record != authority:
        raise RuntimeError('record differs from exact candidate authority')
    for group,count in [('protected',33),('baselines',10),('special',3)]:
        if len(record[group])!=count:
            raise RuntimeError('wrong topology')
        for path,expected in record[group].items():
            if manifest.digest(root/path)!=expected:
                raise RuntimeError('source/artifact mismatch: '+path)


def unique_object(pairs):
    result={}
    for key,value in pairs:
        if key in result: raise RuntimeError('duplicate manifest key: '+key)
        result[key]=value
    return result


def main():
    if manifest.digest(AUTHORITY)!=AUTHORITY_SHA256: raise RuntimeError('candidate authority changed')
    authority=json.loads(AUTHORITY.read_text(),object_pairs_hook=unique_object)
    p=argparse.ArgumentParser();p.add_argument('--entry',type=Path);p.add_argument('--write-entry',type=Path);args=p.parse_args()
    if args.entry: check(ROOT,json.loads(args.entry.read_text(),object_pairs_hook=unique_object),authority)
    old=manifest.protected_authority(ROOT)
    if set(old)!=set(authority['protected']): raise RuntimeError('protected topology changed')
    if manifest.digest(ROOT/'docs/validation/PHASE5D1_FROZEN_CONTEXT_PROVENANCE_PLAN.md')!=authority['historical_plan_sha256']: raise RuntimeError('historical authority changed')
    changed={p for p in old if old[p]!=authority['protected'][p]}
    if changed!={t['path'] for t in authority['transitions']}: raise RuntimeError('transition inventory changed')
    for t in authority['transitions']:
        if t['old']!=old[t['path']] or t['new']!=authority['protected'][t['path']]: raise RuntimeError('transition hash mismatch')
    historical=[dict(path=p,expected=h,actual=manifest.digest(ROOT/p)) for p,h in old.items() if manifest.digest(ROOT/p)!=h]
    check(ROOT,authority,authority)
    if args.write_entry:
        with args.write_entry.open('x') as f: f.write(json.dumps(authority,indent=2)+'\n')
    refused=[]
    try: json.loads('{"protected":{},"protected":{}}',object_pairs_hook=unique_object)
    except RuntimeError: refused.append('duplicate_path_json')
    else: raise RuntimeError('duplicate JSON key escaped')
    with tempfile.TemporaryDirectory() as tmp:
        scratch=Path(tmp)
        for group in ('protected','baselines','special'):
            for path in authority[group]:
                dest=scratch/path;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/path,dest)
        check(scratch,authority,authority)
        mutations={
            'missing_protected': lambda d: d['protected'].pop(next(iter(d['protected']))),
            'extra_protected': lambda d: d['protected'].update({'extra':'0'*64}),
            'forged_hash': lambda d: d['protected'].update({next(iter(d['protected'])):'0'*64}),
            'missing_baseline': lambda d: d['baselines'].pop(next(iter(d['baselines']))),
            'missing_special': lambda d: d['special'].pop(next(iter(d['special']))),
            'forged_status': lambda d: d.update(status='ACCEPTED'),
            'forged_transition': lambda d: d['transitions'][0].update(old='0'*64),
        }
        for name,mutate in mutations.items():
            d=copy.deepcopy(authority);mutate(d)
            try: check(scratch,d,authority)
            except RuntimeError: refused.append(name)
            else: raise RuntimeError('escaped mutation: '+name)
        # Every protected path: real byte mutation after entry, then restoration.
        for path in authority['protected']:
            p=scratch/path;original=p.read_bytes();p.write_bytes(original+b'\n')
            try: check(scratch,authority,authority)
            except RuntimeError: refused.append('source:'+path)
            else: raise RuntimeError('escaped source mutation: '+path)
            p.write_bytes(original)
        check(scratch,authority,authority)
    print(json.dumps(dict(historical_status='FAIL',historical_mismatches=historical,
        proposed_status='PASS_CANDIDATE_NOT_OWNER_RATIFIED',authority_sha256=AUTHORITY_SHA256,
        protected_paths=33,baselines=10,special_artifacts=3,refused_mutations=refused),indent=2))

if __name__=='__main__': main()
