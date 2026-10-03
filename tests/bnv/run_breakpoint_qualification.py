#!/usr/bin/env python3
"""Authenticate immutable inputs and atomically publish a complete bounded run."""
import argparse,csv,hashlib,json,pathlib,shutil,subprocess,time

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def tree(p):
    if p.is_file():return digest(p)
    h=hashlib.sha256()
    for f in sorted(x for x in p.rglob('*') if x.is_file()):
        n=f.relative_to(p).as_posix().encode();h.update(len(n).to_bytes(8,'big'));h.update(n);h.update(bytes.fromhex(digest(f)))
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=pathlib.Path,required=True);ap.add_argument('--mode',default='all',choices=['all','source','control','bounded','crossings','tail']);a=ap.parse_args()
    repo=pathlib.Path(__file__).resolve().parents[2];auth=json.loads((repo/'docs/validation/evidence/obs32-forensics/authentication.json').read_text())
    assert not a.root.exists(),'fresh root required';assert shutil.disk_usage(repo).free>3*1024**3,'insufficient output/build headroom'
    p={k:pathlib.Path(v) for k,v in auth['input_paths'].items()}
    actual={k:tree(v) for k,v in p.items()};assert actual==auth['input_hashes'],actual
    a.root.mkdir(parents=True)
    executable=a.root/'qualification-executable'
    shutil.copy2(repo/'build/obs32/tests/phase6_breakpoint_qualification',executable)
    executable_hash=digest(executable)
    source_hashes={str(f.relative_to(repo)):digest(f) for base in ['CompactStar/Physics/BNV','tests/bnv'] for f in (repo/base).rglob('*') if f.is_file() and f.suffix in ['.hpp','.cpp','.py']}
    source_hashes.update({str(f.relative_to(repo)):digest(f) for f in [repo/'CompactStar/Physics/Evolution/StarContext.hpp',repo/'CompactStar/Physics/Evolution/src/StarContext.cpp']})
    bundle=pathlib.Path('/Users/keeper/Documents/CompactStar/data/phase6-breakpoint-candidate/obs32-review-evidence')
    for line in (bundle/'SHA256SUMS').read_text().splitlines():
        h,n=line.split(maxsplit=1);assert digest(bundle/n.lstrip('*'))==h,n
    cmd=[str(executable),str(p['profile_tree']),str(p['certificate']),str(p['thermal']),str(a.root/'work'),str(p['entry']),str(p['frozen']),str(p['coefficients']),str(bundle/'inputs'),str(a.root/'output.partial'),a.mode]
    started=time.time()
    with (a.root/'console.log').open('w') as log:r=subprocess.run(cmd,cwd=repo,stdout=log,stderr=subprocess.STDOUT)
    record={'command':cmd,'returncode':r.returncode,'seconds':time.time()-started,'inputs':actual,'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),'executable_sha256':executable_hash,'source_sha256':source_hashes,'published':False}
    if r.returncode==0:
        files=list((a.root/'output.partial').glob('*.tsv'));assert files
        for f in files:
            content=f.read_bytes();assert content.endswith(b'\n');rows=list(csv.reader(content.decode().splitlines(),delimiter='\t'));assert all(len(row)==len(rows[0]) for row in rows)
        (a.root/'output.partial').rename(a.root/'output');record['published']=True;record['outputs']={f.name:digest(f) for f in (a.root/'output').iterdir()}
    (a.root/'execution.json').write_text(json.dumps(record,indent=2)+'\n')
    executable.unlink()
    print(json.dumps(record,indent=2));raise SystemExit(r.returncode)
if __name__=='__main__':main()
