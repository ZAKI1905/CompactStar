#!/usr/bin/env python3
"""Incremental source/evidence handoff for office push; never exports credentials."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
root=Path(__file__).resolve().parents[2];records={}
for name,repo,base in [('CompactStar',root,'f4e7c18f568045febc5bf65c6e95200752196dad'),('CONFIND',root.parent/'CONFIND','b0cbd510fd3fd0c772fa50499cd749287cb39e7b'),('Zaki',root.parent/'Zaki','e263a6e180c5c417198e7778bd21fc9c0a32dc33')]:
 def git(*args):return subprocess.check_output(['git',*args],cwd=repo,text=True).strip()
 if git('status','--porcelain'):raise SystemExit('refusing dirty handoff: '+name)
 branch=git('branch','--show-current');head=git('rev-parse','HEAD');bundle=a.output/(name+'.bundle')
 subprocess.run(['git','bundle','create',str(bundle),base+'..'+branch],cwd=repo,check=True)
 result=subprocess.run(['git','bundle','verify',str(bundle)],cwd=repo,capture_output=True,text=True,check=True)
 (a.output/(name+'-verification.txt')).write_text(result.stdout+result.stderr)
 records[name]=dict(base=base,head=head,branch=branch,bundle=bundle.name,sha256=hashlib.sha256(bundle.read_bytes()).hexdigest())
(a.output/'identities.json').write_text(json.dumps(records,indent=2)+'\n')
(a.output/'OFFICE_PUSH.txt').write_text('''Incremental bundles require the recorded base commit in each office repository.
Verify SHA256SUMS first. For each repository, run from its authenticated checkout:
  git bundle verify /path/to/NAME.bundle
  git fetch /path/to/NAME.bundle refs/heads/BRANCH:refs/heads/BRANCH
  git push origin BRANCH:refs/heads/BRANCH
Replace NAME and BRANCH with identities.json entries. Never force-push, merge,
reset the owner checkout, or promote a baseline. If the named branch already
exists, inspect ancestry before updating it. These are candidate corrections
and dependency test overlays, not accepted dependency master changes.
Use the office transfer workflow; no home-Mac-to-cluster scp is required.
''')
(a.output/'SHA256SUMS').write_text(''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.name+'\n' for f in sorted(a.output.iterdir()) if f.is_file()))
zpath=a.output.with_suffix('.zip')
with zipfile.ZipFile(zpath,'x',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(a.output.iterdir()):
  if f.is_file():z.write(f,a.output.name+'/'+f.name)
zpath.with_suffix('.zip.sha256').write_text(hashlib.sha256(zpath.read_bytes()).hexdigest()+'  '+zpath.name+'\n')
print(zpath)
