#!/usr/bin/env python3
"""Two-pass byte verification of completed external evidence, with compact roots.

Full manifests remain beside the heavy cluster artifacts. A compact receipt and
per-top-level-directory roots are suitable for Git; no active output is sealed.
"""
import concurrent.futures,gzip,hashlib,json,os,pathlib,sys,time
from run import save,digest
roots=[pathlib.Path(p).resolve() for p in sys.argv[2:]];destination=pathlib.Path(sys.argv[1]);destination.mkdir(exist_ok=False)
workers=int(os.environ.get('RF_SEAL_WORKERS','24'));start=time.monotonic()
def inventory():
 files=[];links=[]
 for root in roots:
  for parent,dirs,names in os.walk(root,followlinks=True):
   for name in dirs+names:
    p=pathlib.Path(parent)/name
    if p.is_symlink():links.append(dict(path=str(p),target=os.readlink(p)))
   for name in names:
    p=pathlib.Path(parent)/name
    if p.is_file():files.append(p)
 return sorted(files),sorted(links,key=lambda r:r['path'])
def one(p):
 before=p.stat();h=digest(p);after=p.stat()
 assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),f'active evidence: {p}'
 return dict(path=str(p),bytes=after.st_size,sha256=h)
files,links=inventory()
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:records=list(pool.map(one,files))
manifest=destination/'manifest.jsonl.gz'
with manifest.open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as f:
 for r in records:f.write((json.dumps(r,sort_keys=True,separators=(',',':'))+'\n').encode())
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:second=list(pool.map(one,files))
assert records==second,'byte verification changed between independent passes'
newfiles,newlinks=inventory();assert files==newfiles and links==newlinks,'evidence inventory changed'
groups={}
for r in records:
 p=pathlib.Path(r['path']);root=next(root for root in roots if p.is_relative_to(root));rel=p.relative_to(root)
 key=root.name+'/'+rel.parts[0];g=groups.setdefault(key,dict(files=0,bytes=0,hash=hashlib.sha256()))
 g['files']+=1;g['bytes']+=r['bytes'];g['hash'].update((str(rel)+'\0'+str(r['bytes'])+'\0'+r['sha256']+'\n').encode())
summary={k:dict(files=v['files'],bytes=v['bytes'],ordered_sha256=v['hash'].hexdigest()) for k,v in groups.items()}
receipt=dict(status='PASS_TWO_COMPLETE_BYTE_VERIFICATION_PASSES',roots=list(map(str,roots)),files=len(records),bytes=sum(r['bytes'] for r in records),manifest=str(manifest),manifest_sha256=digest(manifest),directory_roots=summary,symlinks=links,workers=workers,elapsed_s=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID'),script_sha256=digest(pathlib.Path(__file__)))
save(destination/'receipt.json',receipt)
print(json.dumps({k:receipt[k] for k in ('status','files','bytes','manifest_sha256','elapsed_s')}))
