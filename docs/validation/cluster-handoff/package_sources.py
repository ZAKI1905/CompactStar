#!/usr/bin/env python3
"""Build source-only v3 bundles: exact commit/tree history plus selected source blobs.

Uses the documented filter=blob:none capability; source blobs are extra objects.
No repository is rewritten or mutated. Not a complete binary-history backup.
"""
import argparse, gzip, hashlib, json, os, pathlib, subprocess
ROOTS = {
 'CompactStar': ('CompactStar/','cmake/','tests/','docs/','main/','dependencies/include/'),
 'Zaki': ('Zaki/','cmake/','tests/','docs/','Examples/'),
 'CONFIND': ('include/','source/','cmake/','tests/','docs/','legal/'),
}
EXCLUDES=('build_xcode/','docs/doxygen_output/','main/Test/results/')
BINARY_SUFFIXES=('.a','.o','.so','.dylib','.dll','.exe','.png','.jpg','.jpeg','.gif','.ico','.pdf','.zip','.tar','.tar.gz','.tgz','.7z','.woff','.ttf','.DS_Store')
SECRET_MARKERS=(b'-----BEGIN PRIVATE KEY-----',b'-----BEGIN RSA PRIVATE KEY-----',b'-----BEGIN OPENSSH PRIVATE KEY-----')
def git(repo,*args,**kw):
 return subprocess.check_output(['git','-C',str(repo),*args],**kw)
def eligible(name,path):
 p=pathlib.PurePosixPath(path)
 return (('/' not in path or path.startswith(ROOTS[name])) and not path.startswith(EXCLUDES)
         and not path.endswith(BINARY_SUFFIXES) and not any(x in p.parts for x in ('.git','.ssh','.aws','.venv','__pycache__','node_modules'))
         and p.name not in ('.env','.DS_Store','credentials','id_rsa','id_ed25519'))
def permitted_content(path,data):
 if any(line in SECRET_MARKERS for line in data.splitlines()): return False
 if path.endswith('.gz'):
  try: data=gzip.decompress(data)
  except (OSError,EOFError): return False
 try: data.decode('utf-8')
 except UnicodeDecodeError: return False
 return b'\0' not in data

def build(name,repo,refs,checkout,out):
 refs={ref:git(repo,'rev-parse',sha+'^{commit}').decode().strip() for ref,sha in refs.items()}
 commits=git(repo,'rev-list',*refs.values()).decode().splitlines()
 lines=git(repo,'rev-list','--objects',*refs.values()).decode().splitlines()
 oid_paths={line.split(' ',1)[0]:line.split(' ',1)[1] if ' ' in line else '' for line in lines}
 process=subprocess.Popen(['git','-C',str(repo),'cat-file','--batch'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
 keep=[]; omitted=[]; blob_info={}
 for oid,path in oid_paths.items():
  process.stdin.write((oid+'\n').encode());process.stdin.flush()
  hdr=process.stdout.readline().decode().split(); size=int(hdr[2]);data=process.stdout.read(size);assert process.stdout.read(1)==b'\n'
  typ=hdr[1]
  if typ!='blob' or eligible(name,path) and permitted_content(path,data):
   keep.append(oid)
   if typ=='blob':blob_info[oid]={'sha256':hashlib.sha256(data).hexdigest(),'bytes':size}
  else:omitted.append({'oid':oid,'representative_path':path,'bytes':size})
 process.stdin.close();process.wait();assert process.returncode==0
 bundle=out/(name+'.bundle')
 with bundle.open('wb') as f:
  f.write(b'# v3 git bundle\n@object-format=sha1\n@filter=blob:none\n')
  for ref,oid in sorted(refs.items()): f.write((oid+' '+ref+'\n').encode())
  f.write(b'\n'); f.flush(); subprocess.run(['git','-C',str(repo),'pack-objects','--stdout','--threads=1','--no-reuse-delta'],input=('\n'.join(sorted(keep))+'\n').encode(),stdout=f,check=True)
 tip=refs[checkout]; included=[]; excluded=[]
 for entry in git(repo,'ls-tree','-rz',tip).split(b'\0'):
  if not entry:continue
  meta,path=entry.split(b'\t',1);mode,typ,oid=meta.decode().split();path=path.decode()
  if oid in blob_info:included.append({'path':path,'mode':mode,'oid':oid,**blob_info[oid]})
  else:excluded.append({'path':path,'mode':mode,'oid':oid})
 result={'name':name,'refs':refs,'checkout_ref':checkout,'checkout_sha':tip,'commit_count':len(commits),'commits_sha256':hashlib.sha256(('\n'.join(sorted(commits))+'\n').encode()).hexdigest(),'included_objects':len(keep),'omitted_blobs':omitted,'included_tip_files':included,'omitted_tip_files':excluded,'bundle_sha256':hashlib.sha256(bundle.read_bytes()).hexdigest(),'filter':'blob:none with extra source blobs; all reachable original commits and trees retained'}
 (out/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
 print(name,'commits',len(commits),'tip files',len(included),'omitted tip',len(excluded),'bundle bytes',bundle.stat().st_size,flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',type=pathlib.Path,required=True);p.add_argument('--compactstar',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=False)
 build('CompactStar',a.compactstar,{'refs/heads/master':'9a5c1eca40758e4c39f9f80db7324e12ad23e505','refs/heads/physics/phase6-breakpoint-aware-evolution':'0c3f9269c884161c9d89163660a47b141a93ace6','refs/heads/physics/phase6-bnv-campaign-readiness':'5433fc2da71ba380264f490b87926d850ff9f375','refs/heads/docs/phase6-cluster-preparation':'HEAD'},'refs/heads/docs/phase6-cluster-preparation',a.output)
 build('Zaki',a.workspace/'external/Zaki',{'refs/heads/master':'7aca56cfdf3875cf2a8c0acea4acdac44f69c9f5','refs/heads/modernize/zaki-2.0.1-fp-preservation':'7aca56cfdf3875cf2a8c0acea4acdac44f69c9f5','refs/handoff/accepted-production':'e263a6e180c5c417198e7778bd21fc9c0a32dc33'},'refs/handoff/accepted-production',a.output)
 build('CONFIND',a.workspace/'external/CONFIND',{'refs/heads/master':'b0cbd510fd3fd0c772fa50499cd749287cb39e7b'},'refs/heads/master',a.output)
if __name__=='__main__':main()
