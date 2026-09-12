#!/usr/bin/env python3
"""Execute qualified controlled tests. Qualification inputs are candidate evidence,
not governed baselines. A missing/mismatched qualification refuses; never skips.
The qualification recipe and original/revised fixed inputs are in docs/validation.
"""
import argparse,hashlib,json,subprocess,tempfile,sys,re
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare(qualification,output):
    profiles=list(qualification.glob('chemical-characterization/run-*/t8192-r80000'))
    if len(profiles)!=1:raise RuntimeError('requires exactly one authenticated radial80000 qualification')
    profile=profiles[0];certificate=qualification/'certificate-t8192-r80000.txt'
    for p,h in [(profile/'freegas.tsv','7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a'),(profile/'profile.tsv','e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890'),(profile/'model.txt','3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4'),(certificate,'7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453')]:
        if digest(p)!=h:raise RuntimeError('qualification input changed: '+str(p))
    thermal=output/'thermal';thermal.mkdir()
    d=np.loadtxt(profile/'profile.tsv',skiprows=1);meta=(profile/'model.txt').read_text().splitlines();m=np.array([float(x.split()[0]) for x in meta[:4]]);hc=float(meta[0].split()[1])
    d=d[d[:,3]>0];d=d[np.argsort(d[:,3])];nb,idx=np.unique(d[:,3],return_index=True);dens=d[idx,6:10]
    pf=hc*np.cbrt(3*np.pi**2*dens);mu=np.hypot(pf,m);slope=np.sum(pf*mu,axis=1)/(3*hc**3*nb)
    if not np.all(np.isfinite(slope)&(slope>0)):raise RuntimeError('invalid entropy adapter')
    temperatures=[0,1,2,4];charges=[0,1]
    for name,values in [('t',temperatures),('nb',nb),('yq',charges)]:
        (thermal/('eos.'+name)).write_text('1 '+str(len(values))+'\n'+'\n'.join(format(x,'.17g') for x in values)+'\n')
    with (thermal/'eos.thermo').open('w') as stream:
        stream.write(f'{m[0]:.17g} {m[1]:.17g} 0\n')
        for it,t in enumerate(temperatures):
            for ib,a in enumerate(slope):
                for iy in range(len(charges)):stream.write(f'{it+1} {ib+1} {iy+1} 0 {a*t:.17g} 0 0 0 0 0 0\n')
    return profile,certificate,thermal

def check_manifest(entry):
    e=json.loads(entry.read_text())
    expected=dict((p,h) for h,p in re.findall(r'^\| `([0-9a-f]{64})` \| `([^`]+)` \|$',(ROOT/'docs/validation/PHASE5D1_FROZEN_CONTEXT_PROVENANCE_PLAN.md').read_text(),re.M))
    observed={x['path']:x['actual'] for x in e['protected']}
    if len(expected)!=33 or len(e['protected'])!=33 or observed!=expected:raise RuntimeError('manifest differs from exact 33-path authority')
    if e['entry_sha']!='f7116c1408c06f976527f86d4397ad6d4540dedf':raise RuntimeError('entry SHA mismatch')
    for path,value in expected.items():
        if digest(ROOT/path)!=value:raise RuntimeError('protected hash mismatch: '+path)
    paths=subprocess.check_output(['git','ls-tree','-r','--name-only',e['entry_sha'],'tests/baselines'],cwd=ROOT,text=True).splitlines()
    if len(paths)!=10 or len(e['baselines'])!=10 or {x['path'] for x in e['baselines']}!=set(paths):raise RuntimeError('baseline manifest incomplete')
    for item in e['baselines']:
        original=subprocess.check_output(['git','show',e['entry_sha']+':'+item['path']],cwd=ROOT)
        if hashlib.sha256(original).hexdigest()!=item['sha256'] or digest(ROOT/item['path'])!=item['sha256']:raise RuntimeError('baseline changed')
    expected_special={'tests/baselines/phase5b_structural_response.json':'7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa','tests/baselines/phase5c_chemical_coefficients.json':'7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7','docs/validation/phase5c_chemical_coefficients_candidate.json':'a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b'}
    if len(e['special'])!=3 or {x['path']:x['expected'] for x in e['special']}!=expected_special:raise RuntimeError('special artifact manifest mismatch')
    for path,value in expected_special.items():
        if digest(ROOT/path)!=value:raise RuntimeError('protected artifact changed')
    return e

def pretrajectory_gate(out,entry):
    def protected():check_manifest(entry)
    protected();codes={}
    for name in ['phase5b_structural_response_regression','phase5c_chemical_coefficient_regression']:
        with (out/(name+'-pretrajectory.log')).open('w') as log:
            codes[name]=subprocess.call(['ctest','--test-dir',str(ROOT/'build'),'-R','^'+name+'$','--output-on-failure','--no-tests=error'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        (out/(name+'-pretrajectory.rc')).write_text(str(codes[name])+'\n')
        if codes[name]:raise RuntimeError('STOP BEFORE TRAJECTORY: '+name+' rc='+str(codes[name]))
    protected()
    (out/'pretrajectory-gate.json').write_text(json.dumps({'protected_count':33,'all_equal_entry':True,'regression_rc':codes,'entry_manifest_sha256':digest(entry)},indent=2)+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('executable',type=Path);p.add_argument('mode',choices=['oracles','trajectory']);p.add_argument('--qualification',type=Path,default=ROOT/'build/phase5d-audit/qualification-fresh-80000');p.add_argument('--output',type=Path);args=p.parse_args()
    out=args.output or Path(tempfile.mkdtemp(prefix='phase5d1c2-',dir=ROOT/'build'))
    out=out.resolve();out.mkdir(exist_ok=True)
    entry=ROOT/'build/phase5d1c2-audit/entry-hashes.json'
    e=check_manifest(entry)
    if args.mode=='trajectory':
        oracle=ROOT/'build/phase5d1c2-audit/oracles-final'
        if (oracle/'raw-rc.txt').read_text().strip()!='0':raise RuntimeError('coupled oracles have not passed')
        if json.loads((oracle/'run-provenance.json').read_text())['executed_sha256']!=digest(args.executable.resolve()):raise RuntimeError('oracle executable differs from trajectory executable')
    for name in ['phase5b_structural_response_regression','phase5c_chemical_coefficient_regression']:
        if (ROOT/('build/phase5d1c2-audit/'+name+'-entry.rc')).read_text().strip()!='0':raise RuntimeError('entry regression not passed')
    profile,certificate,thermal=prepare(args.qualification.resolve(),out)
    command=[str(args.executable.resolve()),str(profile),str(certificate),str(thermal),str(out),args.mode,str(entry)]
    provenance={'command':command,'executable_sha256':digest(args.executable.resolve()),'inputs':{str(p):digest(p) for p in [profile/'freegas.tsv',profile/'profile.tsv',profile/'model.txt',certificate,*thermal.glob('eos.*')]},'cmake_cache_sha256':digest(ROOT/'build/CMakeCache.txt'),'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}
    (out/'run-provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('QUALIFIED COMMAND',json.dumps(command),flush=True)
    # Keep the exact executed image even when later source/build batches change.
    import shutil
    executed=out/'phase5d_evolution.executed';shutil.copy2(args.executable.resolve(),executed);command[0]=str(executed)
    provenance['executed_image']=str(executed);provenance['executed_sha256']=digest(executed)
    (out/'run-provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    with (out/'console.log').open('w') as log:
        process=subprocess.Popen(command,cwd=ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        gate_error=None
        for line in process.stdout:
            log.write(line);log.flush()
            if line.strip()=='PRE_TRAJECTORY_READY':
                try:
                    pretrajectory_gate(out,entry)
                    process.stdin.write('PROTECTED_GOVERNED_GATES_PASS\n');process.stdin.flush()
                except Exception as error:
                    gate_error=str(error);process.stdin.close()
        rc=process.wait()
        if gate_error:
            (out/'gate-error.txt').write_text(gate_error+'\n');rc=rc or 1
    (out/'raw-rc.txt').write_text(str(rc)+'\n');print((out/'console.log').read_text(),flush=True)
    print('RAW EXECUTABLE RC',rc,'EVIDENCE',out,flush=True)
    return rc
if __name__=='__main__':raise SystemExit(main())
