#!/usr/bin/env python3
"""Stage B scheduling only. Refuses mutated bindings and failed preflight."""
import argparse,json,os,pathlib,shlex,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools/linux/preliminary_pair'))
from run import BASE,RETAINED,SNAP,load_pair,digest,save,runtime_binding

def build(out):
 pair=load_pair(SNAP);target=out/'build';target.mkdir(exist_ok=False)
 original=(ROOT/'tests/native/stage.cpp').read_text()
 substitutions={
 'void stage_role(Ctx p,double end,const std::filesystem::path& out) {':'void stage_role(Ctx p,double end,const std::filesystem::path& out,int selected_reference,int selected_tier) {',
 'for(int i=0;i<5;++i) {':'for(int i=0;i<5;++i) { if(i!=selected_reference)continue;',
 'for(int tier=0;tier<3;++tier)for(int policy=0;policy<4;++policy) {':'for(int tier=0;tier<3;++tier)for(int policy=0;policy<4;++policy) { if(tier!=selected_tier)continue;',
 }
 generated=original
 for old,new in substitutions.items():
  assert generated.count(old)==1,old
  generated=generated.replace(old,new)
 (target/'stage_b.inc').write_text(generated)
 flags_path=RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make';link_path=RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt'
 flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in flags_path.read_text().splitlines() if ' = ' in l}
 link=shlex.split(link_path.read_text());obj=target/'driver.o';binary=target/'driver'
 commands=[[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(ROOT/'tests/native'),'-I'+str(target),*flags['CXX_FLAGS'],'-c',str(ROOT/'tools/linux/stage_b/driver.cpp'),'-o',str(obj)],
 [str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]]
 for c in commands:subprocess.run(c,cwd=RETAINED,check=True)
 assert runtime_binding(binary)==pair[0]['runtime_bindings'][pair[0]['stage_binary']]
 files=[ROOT/'tests/native/stage.cpp',ROOT/'tools/linux/stage_b/driver.cpp',pathlib.Path(__file__),target/'stage_b.inc',binary,flags_path,link_path]
 files += [(RETAINED/x).resolve() for x in link if x.endswith(('.a','.so'))]
 save(target/'build.json',dict(status='PASS',commands=commands,substitutions=substitutions,files={str(p):digest(p) for p in files},production_changed=False))

def run(out,kind,index):
 pair=load_pair(SNAP)
 assert json.loads((out/'preflight.json').read_text())['status']=='PASS'
 role='source' if index<(5 if kind=='reference' else 3) else 'control';i=index%(5 if kind=='reference' else 3)
 s=pair[role=='control'];build=out/'build'
 for p,h in json.loads((build/'build.json').read_text())['files'].items():assert digest(p)==h,p
 work=out/f'{role}-{kind}-{i}';work.mkdir(exist_ok=False)
 args=[str(build/'driver'),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(work/'assembly'),str(out/'native'/role),role,str(i if kind=='reference' else -1),str(i if kind=='tier' else -1)]
 save(work/'invocation.json',dict(preflight_sha256=digest(out/'preflight.json'),args=args,snapshot=digest(SNAP/(role+'.json')),build=digest(build/'build.json'),commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),job=os.environ.get('SLURM_JOB_ID')))
 start=time.monotonic()
 with (work/'run.log').open('w') as f:code=subprocess.run(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT).returncode
 save(work/'completion.json',dict(exit_code=code,wall_seconds=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID'),affinity=sorted(os.sched_getaffinity(0))))
 if code:raise SystemExit(code)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','reference','tier']);p.add_argument('out',type=pathlib.Path);p.add_argument('index',nargs='?',type=int);a=p.parse_args()
 if a.mode=='build':build(a.out)
 else:run(a.out,a.mode,a.index)
