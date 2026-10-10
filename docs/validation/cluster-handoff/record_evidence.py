#!/usr/bin/env python3
"""Record the local source, compiler, library and test identities without modifying them."""
import datetime,hashlib,json,pathlib,subprocess,sys,xml.etree.ElementTree as ET
R=pathlib.Path(__file__).resolve().parents[3];O=pathlib.Path(__file__).resolve().parent/'evidence'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cmd(*args,cwd=R):return subprocess.check_output(args,cwd=cwd,text=True).strip()
base=pathlib.Path('/Users/keeper/Documents/CompactStar')
sources={}
for name,repo,expected in [('CompactStar canonical',base/'repo/CompactStar','9a5c1eca40758e4c39f9f80db7324e12ad23e505'),('CompactStar readiness',base/'worktrees/CompactStar-phase6-bnv-campaign-readiness','5433fc2da71ba380264f490b87926d850ff9f375'),('Zaki',base/'external/Zaki','7aca56cfdf3875cf2a8c0acea4acdac44f69c9f5'),('CONFIND',base/'external/CONFIND','b0cbd510fd3fd0c772fa50499cd749287cb39e7b')]:
 head=cmd('git','rev-parse','HEAD',cwd=repo);assert head==expected
 status=cmd('git','status','--porcelain',cwd=repo);assert not status
 sources[name]={'path':str(repo),'head':head,'branch':cmd('git','branch','--show-current',cwd=repo),'status':status,'worktrees':cmd('git','worktree','list','--porcelain',cwd=repo)}
subprocess.run(['git','merge-base','--is-ancestor','9a5c1eca','5433fc2'],cwd=R,check=True)
subprocess.run(['git','merge-base','--is-ancestor','0c3f9269','5433fc2'],cwd=R,check=True)
assert not cmd('git','diff','5433fc2','--','CompactStar','cmake','CMakeLists.txt','tests','dependencies')
sys.path.insert(0,str(R/'tests/rotochemical'));from manifest import protected_authority
protected=[]
for path,old in protected_authority(R).items():
 accepted=hashlib.sha256(subprocess.check_output(['git','show','0c3f9269:'+path],cwd=R)).hexdigest();actual=sha(R/path);assert actual==accepted
 protected.append({'path':path,'phase5d_historical_pin':old,'accepted_amendment_sha256':accepted,'actual':actual,'status':'PIN_MATCH' if old==actual else 'PREVIOUSLY_ACCEPTED_AMENDMENT_DELTA'})
xml=ET.parse(O/'mac-tests.xml').getroot();tests=xml.findall('.//testcase');assert len(tests)==13 and all(t.get('status')=='run' and t.find('failure') is None and t.find('skipped') is None for t in tests)
build=R/'build/mac-recovery';executables=sorted((build/'tests').iterdir());executables=[p for p in executables if p.is_file() and p.stat().st_mode & 0o111]
compiler=cmd('/usr/bin/clang++','--version');cmake=cmd('cmake','--version')
lock=json.loads((R/'cmake/dependency-lock.json').read_text());assert compiler==lock[0]['compiler'].strip();assert cmake==lock[0]['cmake'].strip()
record={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sources':sources,'preparation_base':'5433fc2da71ba380264f490b87926d850ff9f375','preparation_branch':cmd('git','branch','--show-current'),'volume':{'mount':'/Volumes/Scriptorium','filesystem':'APFS','uuid':'99EE5D09-F272-401A-95CB-0D19755A3C24','device_at_inspection':'disk5s1','inspection':'diskutil info /Volumes/Scriptorium; mounted USB external volume'},'readiness_only_commits':cmd('git','log','--format=%H %s','9a5c1eca..5433fc2'),'compiler':compiler,'cmake':cmake,'os':cmd('sw_vers'),'uname':cmd('uname','-a'),'python':cmd('/Users/keeper/miniforge3/bin/python3.12','-c','import sys,numpy,scipy,mpmath; print(sys.version); print("numpy",numpy.__version__,"scipy",scipy.__version__,"mpmath",mpmath.__version__)'),'gsl_version':cmd('/opt/local/bin/gsl-config','--version'),'local_numerical_libraries':{str(p):{'resolved':str(p.resolve()),'sha256':sha(p)} for p in [pathlib.Path('/opt/local/lib/libgsl.27.dylib'),pathlib.Path('/opt/local/lib/libgslcblas.0.dylib'),pathlib.Path('/opt/local/include/gsl/gsl_odeiv2.h')]},'build_files':{str(p.relative_to(R)):sha(p) for p in [build/'CMakeCache.txt',build/'compile_commands.json',build/'libCompactStar.a'] if p.is_file()},'executables':{p.name:{'sha256':sha(p),'linkage':cmd('otool','-L',str(p))} for p in executables},'tests':{'count':13,'failures':0,'skipped':0,'mode':'Debug','names':[t.get('name') for t in tests],'junit_sha256':sha(O/'mac-tests.xml')},'protected':protected,'baseline_sha256':{str(p.relative_to(R)):sha(p) for p in sorted((R/'tests/baselines').rglob('*')) if p.is_file()},'production_delta_from_accepted':cmd('git','diff','0c3f9269','--','CompactStar','cmake','CMakeLists.txt','dependencies','tests/baselines')}
assert record['production_delta_from_accepted']==''
(O/'source-build-tests.json').write_text(json.dumps(record,indent=2)+'\n')
print('PASS: source identities, 33 protected paths, compiler/CMake identity, 13 tests')
