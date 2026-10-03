#!/usr/bin/env python3
"""Authenticate immutable inputs and atomically publish a complete bounded run."""
import argparse,csv,hashlib,json,pathlib,shutil,subprocess,time,math

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
def validate_output(root,mode):
    windows=['source','control','bounded','tail'] if mode=='all' else [mode] if mode in ['source','control','bounded','tail'] else []
    required=[]
    for name in windows:
        required += [name+'-'+kind+'.tsv' for kind in ['history','uncertainty','observations','statistics','events','summary','ledger','envelopes']]
        n={'source':65,'control':65,'bounded':241,'tail':33}[name]
        with (root/(name+'-observations.tsv')).open() as f: observations=list(csv.DictReader(f,delimiter='\t'))
        assert len(observations)==3*n, 'observation completeness'
        for tier in range(3):
            group=[r for r in observations if r['tier']==str(tier)]
            assert len(group)==n, 'tier completeness'
            times=[float(r['time']) for r in group]
            assert all(math.isfinite(t) for t in times) and all(a<b for a,b in zip(times,times[1:])), 'time corruption'
            assert all(math.isfinite(float(r[k])) for r in group for k in ['x','e','m']), 'nonfinite state'
    if mode in ['all','events']:required+=['historical-reproduction.tsv','reconstruction-events.tsv']
    if mode in ['all','events','crossings']:
        required+=['crossings.tsv']
        with (root/'crossings.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
        assert len(r)==333 and len({(x['run'],x['step']) for x in r})==37, 'crossing completeness'
    if mode=='supplement':
        required+=['source-run4-comparison.tsv','control-run4-comparison.tsv','tail-jacobian.tsv','tail-stability.tsv']
        for name,n in [('source',490),('control',1277)]:
            with (root/(name+'-run4-comparison.tsv')).open() as f:r=list(csv.DictReader(f,delimiter='\t'))
            assert len(r)==3*n and len({x['observation'] for x in r})==n, 'run4 completeness'
    assert required and all((root/n).is_file() for n in required), 'required artifact missing'
    for f in root.glob('*.tsv'):
        content=f.read_bytes();assert content.endswith(b'\n'), 'unterminated artifact'
        r=list(csv.reader(content.decode().splitlines(),delimiter='\t'));assert all(len(row)==len(r[0]) for row in r), 'schema corruption'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=pathlib.Path,required=True);ap.add_argument('--mode',default='all',choices=['all','source','control','bounded','crossings','tail','events','supplement']);a=ap.parse_args()
    repo=pathlib.Path(__file__).resolve().parents[2];auth=json.loads((repo/'docs/validation/evidence/obs32-forensics/authentication.json').read_text())
    assert not a.root.exists(),'fresh root required';assert shutil.disk_usage(repo).free>3*1024**3,'insufficient output/build headroom'
    p={k:pathlib.Path(v) for k,v in auth['input_paths'].items()}
    actual={k:tree(v) for k,v in p.items()};assert actual==auth['input_hashes'],actual
    a.root.mkdir(parents=True)
    executable=a.root/'qualification-executable'
    shutil.copy2(repo/'build/obs32/tests/phase6_breakpoint_qualification',executable)
    executable_hash=digest(executable)
    source_hashes={str(f.relative_to(repo)):digest(f) for base in ['CompactStar/Physics/BNV','tests/bnv'] for f in (repo/base).rglob('*') if f.is_file() and f.suffix in ['.hpp','.cpp','.py']}
    source_hashes.update({str(f.relative_to(repo)):digest(f) for f in [repo/'CompactStar/Physics/Evolution/StarContext.hpp',repo/'CompactStar/Physics/Evolution/src/StarContext.cpp',repo/'CompactStar/Analysis/src/ParticleNumberResponse.cpp',repo/'CompactStar/Analysis/src/ChemicalResponse.cpp']})
    bundle=pathlib.Path('/Users/keeper/Documents/CompactStar/data/phase6-breakpoint-candidate/obs32-review-evidence')
    for line in (bundle/'SHA256SUMS').read_text().splitlines():
        h,n=line.split(maxsplit=1);assert digest(bundle/n.lstrip('*'))==h,n
    cmd=[str(executable),str(p['profile_tree']),str(p['certificate']),str(p['thermal']),str(a.root/'work'),str(p['entry']),str(p['frozen']),str(p['coefficients']),str(bundle/'inputs'),str(a.root/'output.partial'),a.mode]
    source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
    started=time.time()
    with (a.root/'console.log').open('w') as log:r=subprocess.run(cmd,cwd=repo,stdout=log,stderr=subprocess.STDOUT)
    record={'command':cmd,'returncode':r.returncode,'seconds':time.time()-started,'inputs':actual,'source_head':source_head,'executable_sha256':executable_hash,'source_sha256':source_hashes,'published':False}
    if r.returncode==0:
        try:
            validate_output(a.root/'output.partial',a.mode)
            (a.root/'output.partial').rename(a.root/'output');record['published']=True;record['outputs']={f.name:digest(f) for f in (a.root/'output').iterdir()}
        except Exception as error:
            record['validation_error']=str(error);record['returncode']=2
    metadata=json.dumps(record,indent=2)+'\n'
    (a.root/'execution.json.partial').write_text(metadata)
    assert (a.root/'execution.json.partial').read_text()==metadata
    (a.root/'execution.json.partial').rename(a.root/'execution.json')
    executable.unlink()
    print(json.dumps(record,indent=2));raise SystemExit(record['returncode'])
if __name__=='__main__':main()
