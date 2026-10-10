#!/usr/bin/env python3
"""Read-only authentication of original Mac packages and frozen Phase-6 inputs."""
import hashlib, json, pathlib, subprocess
R = pathlib.Path(__file__).resolve().parents[3]
O = pathlib.Path(__file__).resolve().parent / 'evidence'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()
def files(p):
    assert p.is_dir(), str(p)
    result=sorted(x for x in p.rglob('*') if x.is_file())
    assert result, str(p)
    return result
def tree(p):
    if p.is_file(): return sha(p)
    h=hashlib.sha256()
    for f in files(p):
        n=f.relative_to(p).as_posix().encode()
        h.update(len(n).to_bytes(8,'big')); h.update(n); h.update(bytes.fromhex(sha(f)))
    return h.hexdigest()
def main():
    lockpath=R/'cmake/dependency-lock.json'
    lock=json.loads(lockpath.read_text())
    packages=[]
    for rec in lock:
        p=pathlib.Path(rec['prefix']).resolve(strict=True)
        def manifest(base):
            return hashlib.sha256(''.join(sha(f)+'  '+f.relative_to(p).as_posix()+'\n' for f in files(p/base)).encode()).hexdigest()
        actual={'archive_sha256':sha(p/'lib'/('lib'+rec['name']+'.a')), 'headers_sha256':manifest('include'),'configs_sha256':manifest('lib/cmake')}
        assert all(actual[k]==rec[k] for k in actual), (rec['name'],rec['mode'],actual)
        prov=json.loads((p/'provenance.json').read_text())
        assert prov==rec, 'complete package provenance differs from lock'
        archived=R/'docs/validation/evidence/external-deps-migration'/f"{rec['name']}-{rec['mode']}-provenance.json"
        assert (p/'provenance.json').read_bytes()==archived.read_bytes(), 'archived provenance bytes differ'
        packages.append({'name':rec['name'],'mode':rec['mode'],'original_prefix':rec['prefix'],'resolved_prefix':str(p),'actual':actual,'provenance_sha256':sha(p/'provenance.json'),'archived_provenance_exact':True,'status':'PASS'})
    for mode in ['Debug','Release']:
        pr={r['name']:r['prefix'] for r in lock if r['mode']==mode}
        subprocess.run(['python3',str(R/'cmake/verify_dependencies.py'),'--mode',mode,'--zaki',pr['Zaki'],'--confind',pr['CONFIND'],'--output',str(O/('original-verifier-'+mode+'.json'))],check=True)
    auth=json.loads((R/'docs/validation/evidence/obs32-forensics/authentication.json').read_text())
    inputs=[]
    for key,original in auth['input_paths'].items():
        p=pathlib.Path(original); actual=tree(p)
        assert actual==auth['input_hashes'][key], (key,actual)
        inputs.append({'key':key,'original':original,'resolved':str(p.resolve(strict=True)),'sha256':actual,'status':'PASS','files':len(files(p)) if p.is_dir() else 1,'bytes':sum(f.stat().st_size for f in files(p)) if p.is_dir() else p.stat().st_size})
    record={'dependency_lock_sha256':sha(lockpath),'package_artifacts':packages,'scientific_inputs':inputs,'authority':'Existing lock and obs32-forensics/authentication.json; original bytes only. No replacement packages built.'}
    (O/'mac-recovery.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'packages':len(packages),'inputs':inputs,'lock_sha256':record['dependency_lock_sha256']},indent=2))

if __name__=="__main__":main()
