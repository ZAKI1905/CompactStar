#!/usr/bin/env python3
"""Offline import of supplied filtered source bundles, with exact sparse checkout."""
import argparse, hashlib, json, os, pathlib, subprocess

def run(repo,*args,**kwargs):
 env={**os.environ,'GIT_NO_LAZY_FETCH':'1','GIT_TERMINAL_PROMPT':'0','GIT_ALLOW_PROTOCOL':'file'}
 return subprocess.run(['git','-C',str(repo),*args],env=env,check=True,**kwargs)
def escape(path):
 for c in ('\\','*','?','[',']','!','#',' '):path=path.replace(c,'\\'+c)
 return '/'+path

def main():
 p=argparse.ArgumentParser();p.add_argument('--bundles',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=False)
 for name in ['CompactStar','Zaki','CONFIND']:
  m=json.loads((a.bundles/(name+'.json')).read_text());bundle=(a.bundles/(name+'.bundle')).read_bytes()
  assert hashlib.sha256(bundle).hexdigest()==m['bundle_sha256'],name
  header,pack=bundle.split(b'\n\n',1)
  assert header.startswith(b'# v3 git bundle\n@object-format=sha1\n@filter=blob:none\n')
  assert dict((ref,oid) for oid,ref in (line.split() for line in header.decode().splitlines()[3:]))==m['refs']
  dest=a.output/name;run(a.output,'init','-q',str(dest))
  # Git requires both a promisor pack and a promisor remote to recognize omitted blobs.
  # This is deliberately a nonexistent LOCAL path; there is no network upstream.
  run(dest,'config','remote.omitted.promisor','true')
  run(dest,'config','remote.omitted.url','./OMITTED_OBJECTS_NOT_SUPPLIED')
  run(dest,'config','remote.omitted.partialclonefilter','blob:none')
  run(dest,'index-pack','--stdin','--fix-thin','--promisor',input=pack)
  for ref,oid in m['refs'].items():run(dest,'update-ref',ref,oid)
  run(dest,'config','core.sparseCheckout','true');run(dest,'config','core.sparseCheckoutCone','false')
  info=dest/'.git/info';info.mkdir(exist_ok=True)
  (info/'sparse-checkout').write_text('\n'.join(escape(f['path']) for f in m['included_tip_files'])+'\n')
  run(dest,'checkout','--detach',m['checkout_sha'])
  assert run(dest,'status','--porcelain',capture_output=True,text=True).stdout=='', name
  for f in m['included_tip_files']:
   q=dest/f['path'];assert hashlib.sha256(q.read_bytes()).hexdigest()==f['sha256'],str(q)
   assert not q.is_symlink(),str(q)
  for f in m['omitted_tip_files']:assert not (dest/f['path']).exists(),f['path']
  commits=run(dest,'rev-list',*m['refs'].values(),capture_output=True,text=True).stdout.splitlines()
  assert len(commits)==m['commit_count']
  assert hashlib.sha256(('\n'.join(sorted(commits))+'\n').encode()).hexdigest()==m['commits_sha256']
  run(dest,'fsck','--full')
  print(name,'PASS: refs, history, source bytes, sparse status, fsck',flush=True)
if __name__=='__main__':main()
