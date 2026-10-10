#!/usr/bin/env python3
"""Verify a manual handoff using Python stdlib only; never execute scientific code."""
import argparse,hashlib,json,pathlib

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def tree(p):
 if p.is_file():return sha(p)
 assert p.is_dir(),str(p)
 files=sorted(q for q in p.rglob('*') if q.is_file());assert files,str(p)
 h=hashlib.sha256()
 for f in files:
  n=f.relative_to(p).as_posix().encode();h.update(len(n).to_bytes(8,'big'));h.update(n);h.update(bytes.fromhex(sha(f)))
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=pathlib.Path,required=True);a=ap.parse_args();r=a.root.resolve();names=[]
 for line in (r/'SHA256SUMS').read_text().splitlines():
  h,n=line.split('  ',1);p=pathlib.PurePosixPath(n)
  assert not p.is_absolute() and '..' not in p.parts,n
  f=r/n;assert f.is_file() and not f.is_symlink(),n;assert sha(f)==h,n;names.append(n)
 actual={p.relative_to(r).as_posix() for p in r.rglob('*') if p.is_file() and p!=r/'SHA256SUMS'}
 assert len(names)==len(set(names)) and set(names)==actual,'file inventory mismatch'
 assert not any(p.is_symlink() for p in r.rglob('*')),'unexpected symlink'
 inputs=json.loads((r/'scientific/input-map.json').read_text())
 for rec in inputs:assert tree(r/'scientific'/rec['relative_path'])==rec['sha256'],rec['key']
 for name in ('CompactStar','Zaki','CONFIND'):
  m=json.loads((r/'bundles'/(name+'.json')).read_text());assert sha(r/'bundles'/(name+'.bundle'))==m['bundle_sha256']
 print('PASS:',len(names),'payload files;',len(inputs),'scientific identities; 3 bundle digests')
if __name__=='__main__':main()
