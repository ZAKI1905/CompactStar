#!/usr/bin/env python3
"""Read-only provenance audit; missing artifacts never hash as empty trees."""
import argparse,gzip,hashlib,json,pathlib,subprocess,sys
sys.dont_write_bytecode=True
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tests/rotochemical'),str(ROOT/'tests/bnv')]
from manifest import protected_authority
from run_breakpoint_qualification import tree

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def command(argv):
    r=subprocess.run(argv,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    return {'command':argv,'returncode':r.returncode,'output':r.stdout}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=pathlib.Path,required=True);a=ap.parse_args()
    out=a.output;out.mkdir(parents=True,exist_ok=True)
    accepted='0c3f9269c884161c9d89163660a47b141a93ace6'
    d={'accepted':accepted,'head':command(['git','rev-parse','HEAD']), 'master':command(['git','rev-parse','master']), 'ancestry':command(['git','log','--format=%H %s','master..'+accepted]),'environment':{}}
    for name,cmd in {'compiler':['/usr/bin/clang++','--version'],'cmake':['cmake','--version'],'gsl':['gsl-config','--version'],'os':['sw_vers'],'kernel':['uname','-a'],'disks':['diskutil','list'],'python':[sys.executable,'--version']}.items(): d['environment'][name]=command(cmd)
    d['protected']=[]
    for name,old in sorted(protected_authority(ROOT).items()):
        actual=digest(ROOT/name)
        reviewed=hashlib.sha256(subprocess.check_output(['git','show',accepted+':'+name],cwd=ROOT)).hexdigest()
        assert actual==reviewed,'accepted protected bytes changed: '+name
        d['protected'].append({'path':name,'historical_phase5d_pin':old,'accepted_amendment_sha256':reviewed,'actual':actual,'status':'HISTORICAL_PIN_MATCH' if old==actual else 'ACCEPTED_AMENDMENT_DELTA'})
    d['all_baselines']={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'tests/baselines').rglob('*')) if p.is_file()}
    d['accepted_implementation_diff']=command(['git','diff',accepted,'--','CompactStar','cmake','CMakeLists.txt','dependencies','tests/baselines'])
    assert not d['accepted_implementation_diff']['output']
    auth=json.loads((ROOT/'docs/validation/evidence/obs32-forensics/authentication.json').read_text());d['inputs']={}
    recovered=json.loads((out/'recovery-input-trees.json').read_text());small=json.loads((out/'recovery-search.json').read_text())['files']
    for k,s in auth['input_paths'].items():
        p=pathlib.Path(s);expected=auth['input_hashes'][k];copies=[]
        for x in recovered:
            if x['kind']==k and x['match']:copies.append(x['path'])
        for x in small:
            if x['sha256']==expected:copies.append(x['path'])
        actual=tree(p) if p.exists() else None
        alternatives=[{'path':s,'actual':tree(pathlib.Path(s))} for s in copies]
        assert all(x['actual']==expected for x in alternatives)
        d['inputs'][k]={'expected':expected,'original_path':str(p),'original_exists':p.exists(),'actual':actual,'authenticated_alternatives':alternatives,'status':'PASS' if actual==expected or alternatives else 'UNQUALIFIED'}
    d['packages']=[]
    for rec in json.loads((ROOT/'cmake/dependency-lock.json').read_text()):
        prefix=pathlib.Path(rec['prefix']);p=prefix/'lib'/('lib'+rec['name']+'.a')
        d['packages'].append({**rec,'exists':p.is_file(),'actual':digest(p) if p.is_file() else None,'status':'UNQUALIFIED' if not p.is_file() else 'REQUIRES_FULL_PACKAGE_AUTHENTICATION'})
    arch=ROOT/'docs/validation/evidence/phase6-breakpoint-candidate';m=json.loads((arch/'manifest.json').read_text());checked=[]
    for rec in m['files']:
        p=arch/rec['path'];assert digest(p)==rec['sha256'],str(p)
        if 'uncompressed_sha256' in rec:assert hashlib.sha256(gzip.decompress(p.read_bytes())).hexdigest()==rec['uncompressed_sha256']
        checked.append(rec['path'])
    d['archived_evidence_verified']=checked
    d['local_gsl']={str(p):digest(p) for p in [pathlib.Path('/opt/local/lib/libgsl.27.dylib'),pathlib.Path('/opt/local/lib/libgslcblas.0.dylib'),pathlib.Path('/opt/local/include/gsl/gsl_odeiv2.h')]}
    d['standalone_compile_command']=['/usr/bin/clang++','-std=c++17','-g','-I.','-I/opt/local/include','tests/bnv/breakpoint_contract.cpp','-L/opt/local/lib','-lgsl','-lgslcblas','-o','build/readiness-standalone/breakpoint_contract']
    d['standalone_contract']=command(['build/readiness-standalone/breakpoint_contract'])
    assert d['standalone_contract']['returncode']==0
    d['standalone_executable_sha256']=digest(ROOT/'build/readiness-standalone/breakpoint_contract')
    d['standalone_source_sha256']=digest(ROOT/'tests/bnv/breakpoint_contract.cpp')
    d['standalone_linkage']=command(['otool','-L','build/readiness-standalone/breakpoint_contract'])
    d['qualification_gap']='13-test qualified-stack rerun unavailable: CMake refused absent Scriptorium package prefix. Standalone GSL test is a distinct build, not a substitute.'
    (out/'environment.json').write_text(json.dumps(d,indent=2)+'\n')
    print('protected',len(d['protected']),'baselines',len(d['all_baselines']),'archive',len(checked),'inputs',{k:v['status'] for k,v in d['inputs'].items()})
if __name__=='__main__':main()
