#!/usr/bin/env python3
"""Authenticated exploratory evidence; retained production files never modified."""
import argparse, concurrent.futures, gzip, hashlib, json, os, pathlib, shlex, subprocess, sys, time
ROOT=pathlib.Path(__file__).resolve().parents[3]
BASE=ROOT.parents[1]
sys.path.insert(0,str(ROOT/'tools/linux/native_pair'))
from snapshot import load_pair, digest
from runtime import runtime_binding
SNAP=ROOT/'docs/validation/results-first-bnv-pilot/parallel-preparation/snapshots'
RETAINED=BASE/'parallel-preparation-20261110-r3/build/tests'

def save(path,x):
 path.write_text(json.dumps(x,indent=2,default=str)+'\n')

def build(out):
 out.mkdir(parents=True,exist_ok=False)
 pair=load_pair(SNAP)
 original=(ROOT/'tests/native/stage.cpp').read_text()
 replacements={
  'void stage_role(Ctx p,double end,const std::filesystem::path& out) {':'void stage_role(Ctx p,double end,const std::filesystem::path& out,int selected_reference,int selected_tier) {',
  'for(int i=0;i<5;++i) {':'for(int i=0;i<5;++i) { if(i!=selected_reference)continue;',
  'const auto& ref=refs.back();':'const auto& ref=refs.back();',
  'for(int tier=0;tier<3;++tier)for(int policy=0;policy<4;++policy) {':'for(int tier=0;tier<3;++tier)for(int policy=0;policy<4;++policy) { if(tier!=selected_tier)continue;',
 }
 generated=original
 for old,new in replacements.items():
  assert generated.count(old)==1,old
  generated=generated.replace(old,new)
 (out/'control_stage.inc').write_text(generated)
 flags_path=RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make'
 link_path=RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt'
 flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in flags_path.read_text().splitlines() if ' = ' in l}
 link=shlex.split(link_path.read_text());commands=[]
 for name in ('control',):
  obj=out/(name+'.o');binary=out/name
  command=[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(ROOT/'tests/native'),'-I'+str(out),*flags['CXX_FLAGS'],'-c',str(ROOT/'tools/linux/preliminary_pair'/f'{name}.cpp'),'-o',str(obj)]
  commands.append(command);subprocess.run(command,cwd=RETAINED,check=True)
  command=[str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]
  commands.append(command);subprocess.run(command,cwd=RETAINED,check=True)
  assert runtime_binding(binary)==pair[0]['runtime_bindings'][pair[0]['stage_binary']]
 paths=[ROOT/'tests/native/stage.cpp',ROOT/'tools/linux/preliminary_pair/control.cpp',pathlib.Path(__file__),out/'control_stage.inc',out/'control',flags_path,link_path]
 paths += [(RETAINED/x).resolve() for x in link if x.endswith(('.a','.so'))]
 save(out/'build.json',dict(status='PASS',scheduling_substitutions=replacements,commands=commands,cwd=RETAINED,files={str(p):digest(p) for p in paths},snapshot_identities=json.loads((SNAP/'identities.json').read_text()),production_changed=False))

def control(out,kind,index):
 s=load_pair(SNAP)[1]
 build=out/'build'
 manifest=json.loads((build/'build.json').read_text())
 for path,sha in manifest['files'].items():assert digest(path)==sha,path
 work=out/f'{kind}-{index}';work.mkdir(exist_ok=False)
 args=[str(build/'control'),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(work/'assembly'),str(out/'control'),str(index if kind=='reference' else -1),str(index if kind=='tier' else -1)]
 save(work/'invocation.json',dict(args=args,control_snapshot_sha256=digest(SNAP/'control.json'),build_sha256=digest(build/'build.json'),exploratory=True,qualified=False))
 start=time.monotonic()
 with (work/'run.log').open('w') as log:
  code=subprocess.run(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT).returncode
 save(work/'completion.json',dict(exit_code=code,wall_seconds=time.monotonic()-start,slurm_job=os.environ.get('SLURM_JOB_ID'),affinity=sorted(os.sched_getaffinity(0))))
 assert code==0,(kind,index,code)

def authenticate(out):
 load_pair(SNAP)
 manifest=ROOT/'docs/validation/results-first-bnv-pilot/parallel-source-assessment/evidence/external-SHA256SUMS.gz'
 with gzip.open(manifest,'rt') as f:entries=[line.rstrip('\n').split('  ',1) for line in f]
 def check(row):
  sha,path=row
  return path if digest(path)!=sha else None
 start=time.monotonic()
 with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:
  failures=[p for p in pool.map(check,entries) if p]
 save(out/'authentication.json',dict(status='FAIL' if failures else 'PASS',files=len(entries),failures=failures,manifest_sha256=digest(manifest),snapshots=json.loads((SNAP/'identities.json').read_text()),wall_seconds=time.monotonic()-start))
 assert not failures

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','authenticate','reference','tier']);p.add_argument('out',type=pathlib.Path);p.add_argument('index',nargs='?',type=int)
 a=p.parse_args()
 if a.mode=='build':build(a.out/'build')
 elif a.mode=='authenticate':authenticate(a.out)
 else:control(a.out,a.mode,a.index)
