#!/usr/bin/env python3
"""Build and run isolated source reconstructions against retained numerical archives."""
import json,os,pathlib,shlex,subprocess,sys,time
from run import ROOT,BASE,RETAINED,SNAP,load_pair,digest,save,runtime_binding

def build(out):
 out.mkdir(exist_ok=False)
 load_pair(SNAP)
 trace=(ROOT/'tests/diagnostics/checkpoint_replay.cpp').read_text()
 reader=(ROOT/'tests/native/parallel_regression.cpp').read_text()
 (out/'retained_helpers.inc').write_text(trace[trace.index('struct TraceContext'):trace.index('std::string Bytes')]+reader[reader.index('BNV::MainTrajectoryResult ReadRetainedMain'):reader.index('void WriteLedger')])
 flags_path=RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make';link_path=RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt'
 flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in flags_path.read_text().splitlines() if ' = ' in l}
 link=shlex.split(link_path.read_text());obj=out/'recover.o';binary=out/'recover'
 commands=[[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(out),*flags['CXX_FLAGS'],'-c',str(ROOT/'tools/linux/preliminary_pair/recover.cpp'),'-o',str(obj)],
 [str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]]
 for c in commands:subprocess.run(c,cwd=RETAINED,check=True)
 paths=[ROOT/'tests/native/stage.cpp',ROOT/'tests/native/parallel_regression.cpp',ROOT/'tests/diagnostics/checkpoint_replay.cpp',ROOT/'tools/linux/preliminary_pair/recover.cpp',pathlib.Path(__file__),out/'retained_helpers.inc',binary,flags_path,link_path]
 paths += [(RETAINED/x).resolve() for x in link if x.endswith(('.a','.so'))]
 save(out/'build.json',dict(commands=commands,files={str(p):digest(p) for p in paths},runtime=runtime_binding(binary),production_changed=False))

def recover(out,mode):
 s=load_pair(SNAP)[0];build=out/'recovery-build'
 for p,h in json.loads((build/'build.json').read_text())['files'].items():assert digest(p)==h,p
 assert json.loads((out/'authentication.json').read_text())['status']=='PASS'
 samples=[0,1,2,4,8,16,32,64,128,256,512,960,1440,1918,1919,1920]
 indices=samples if mode=='sample' else [i for i in range(1921) if i not in samples]
 if mode=='full':
  assert json.loads((out/'sample-completion.json').read_text())['exit_code']==0
 index_file=out/(mode+'-indices.txt');index_file.write_text('\n'.join(map(str,indices))+'\n')
 retained=BASE/'stage-a-parallel-20261110-r1/source-stage/run/trajectory/source/tier-0-cap-0'
 args=[str(build/'recover'),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(out/(mode+'-assembly')),str(retained),str(out/mode),str(index_file)]
 save(out/(mode+'-invocation.json'),dict(args=args,snapshot=digest(SNAP/'source.json'),build=digest(build/'build.json')))
 start=time.monotonic()
 code=subprocess.run(args,cwd=ROOT).returncode
 save(out/(mode+'-completion.json'),dict(exit_code=code,wall_seconds=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID')))
 assert code==0

if __name__=='__main__':
 mode=sys.argv[1];out=pathlib.Path(sys.argv[2])
 if mode=='build':build(out/'recovery-build')
 else:recover(out,mode)
