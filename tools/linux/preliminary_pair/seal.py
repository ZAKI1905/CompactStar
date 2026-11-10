#!/usr/bin/env python3
"""Hash and independently reread all external artifacts; preserve earlier attempts.

Run only after numerical jobs have exited. Destination must be outside the run.
The manifest and receipt remain external; commit the compact receipt to Git.
"""
import concurrent.futures,gzip,json,os,pathlib,sys,time
from run import digest,save,authenticate
out=pathlib.Path(sys.argv[1]).resolve();destination=pathlib.Path(sys.argv[2]).resolve()
assert out not in destination.parents and out!=destination
destination.mkdir(exist_ok=False)
post=out/'post-authentication';post.mkdir(exist_ok=False);authenticate(post)
start=time.monotonic()
paths=sorted(p for p in out.rglob('*') if p.is_file())
def record(p):return digest(p),str(p),p.stat().st_size
with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:entries=list(pool.map(record,paths))
manifest=destination/'external-SHA256SUMS.gz'
with manifest.open('wb') as raw:
 with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as f:
  for sha,path,size in entries:f.write(f'{sha}  {path}\n'.encode())
def verify(entry):
 sha,path,size=entry;p=pathlib.Path(path)
 return path if p.stat().st_size!=size or digest(p)!=sha else None
with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:failures=[x for x in pool.map(verify,entries) if x]
assert sorted(p for p in out.rglob('*') if p.is_file())==paths,'inventory changed during sealing'
receipt=dict(status='PASS' if not failures else 'FAIL',external_root=str(out),manifest=str(manifest),manifest_sha256=digest(manifest),files=len(entries),total_bytes=sum(x[2] for x in entries),verification='Every file reread and SHA-256 compared in a second complete pass; file set and sizes checked.',workers=24,elapsed_s=time.monotonic()-start,source_post_authentication=json.loads((post/'authentication.json').read_text()),failures=failures,script_sha256=digest(pathlib.Path(__file__)))
save(destination/'receipt.json',receipt);assert not failures
print(json.dumps({k:receipt[k] for k in ('status','files','total_bytes','manifest_sha256','elapsed_s')},indent=2))
