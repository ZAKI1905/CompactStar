#!/usr/bin/env python3
"""Extend accepted Stage A actual-observation tracing to all Stage B endpoints."""
import json,os,pathlib,shlex,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools/linux/preliminary_pair'))
from run import BASE,RETAINED,SNAP,load_pair,digest,save,runtime_binding
out=pathlib.Path(sys.argv[1]);mode=sys.argv[2];build=out/'observation-trace-build-v2';pair=load_pair(SNAP)
if mode=='build':
 build.mkdir(exist_ok=False)
 old=ROOT/'tools/linux/preliminary_pair/interpolation.cpp'
 src=old.read_text().replace('#include "../../../tests/native/stage.cpp"','#include "'+str(ROOT/'tests/native/stage.cpp')+'"').replace('"retained_helpers.inc"','"trace_context.inc"')
 extra=r'''
  // Append accepted main and restart times not already on the scientific grid.
  const auto supplement=out.parent_path()/(role+"-supplement");
  require(std::filesystem::exists(supplement),"accepted endpoint evidence absent");
  for(const auto& file:std::filesystem::directory_iterator(supplement)){
   if(file.path().filename().string().rfind("accepted-",0)!=0||file.path().extension()!=".tsv")continue;
   std::ifstream f(file.path());std::string line;std::getline(f,line);double target=0;std::vector<State> ys;
   while(std::getline(f,line)){std::istringstream r(line);int k;double t;State y;r>>k>>t;for(auto& v:y)r>>v;require(bool(r)&&k==int(ys.size()),"bad accepted reference row");if(k==0)target=t;else require(t==target,"accepted reference times differ");ys.push_back(y);}
   require(ys.size()==5,"missing accepted references");
   if(std::find(times.begin(),times.end(),target)!=times.end())continue;
   times.push_back(target);for(int k=0;k<5;++k)observations[k].push_back(ys[k]);
  }
  const auto count=times.size();
'''
 needle='  NativeQualification::ProcessPool(5*1921,24,out/"workers",'
 assert src.count(needle)==1
 src=src.replace(needle,extra+'  NativeQualification::ProcessPool(5*count,StageWorkers(),out/"workers",')
 src=src.replace('const auto k=task/1921,index=task%1921;','const auto k=task/count,index=task%count;')
 src=src.replace('PASS all 9605 reference observations reproduced exactly','PASS all scientific and accepted reference observations reproduced exactly')
 src=src.replace('#include "trace_context.inc"','#include "trace_context.inc"\n#include "AuditStepper.hpp"\nusing AuditStepper=BNV::AuditStepper;')
 src=src.replace('TraceContext context(p,1,rhs,diag);','std::ofstream applies(dir/"applies.tsv");TraceContext context(p,1,rhs,diag);BindAudit(context,applies);').replace('Stepper solver(','AuditStepper solver(')
 (build/'observations.cpp').write_text(src)
 from audit_trace import generate
 generate(ROOT,build)
 flags_path=RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make';link_path=RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt'
 flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in flags_path.read_text().splitlines() if ' = ' in l}
 link=shlex.split(link_path.read_text());obj=build/'observations.o';binary=build/'observations'
 commands=[[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(build),*flags['CXX_FLAGS'],'-c',str(build/'observations.cpp'),'-o',str(obj)],
 [str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]]
 for cmd in commands:subprocess.run(cmd,cwd=RETAINED,check=True)
 assert runtime_binding(binary)==pair[0]['runtime_bindings'][pair[0]['stage_binary']]
 files=[old,ROOT/'tests/native/stage.cpp',ROOT/'tests/diagnostics/checkpoint_replay.cpp',pathlib.Path(__file__),build/'observations.cpp',build/'trace_context.inc',build/'AuditStepper.hpp',ROOT/'tools/linux/stage_b/audit_trace.py',ROOT/'CompactStar/Physics/BNV/BreakpointStepper.hpp',binary,flags_path,link_path]
 files += [(RETAINED/x).resolve() for x in link if x.endswith(('.a','.so'))]
 save(build/'build.json',dict(commands=commands,files={str(p):digest(p) for p in files}))
else:
 assert mode in ('source','control');s=pair[mode=='control']
 for p,h in json.loads((build/'build.json').read_text())['files'].items():assert digest(p)==h,p
 args=[str(build/'observations'),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(out/(mode+'-observation-assembly')),str(out/'native'/mode),str(out/(mode+'-interpolation')),mode]
 save(out/(mode+'-observation-invocation.json'),dict(args=args,build=digest(build/'build.json'),job=os.environ.get('SLURM_JOB_ID')))
 start=time.monotonic();code=subprocess.run(args,cwd=ROOT).returncode
 save(out/(mode+'-observation-completion.json'),dict(exit_code=code,elapsed_s=time.monotonic()-start));assert code==0
