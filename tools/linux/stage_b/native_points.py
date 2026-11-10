#!/usr/bin/env python3
"""Build/dispatch retained-reference evidence, no trajectory reruns."""
import json,os,pathlib,shlex,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools/linux/preliminary_pair'))
from run import BASE,RETAINED,SNAP,load_pair,digest,save,runtime_binding
out=pathlib.Path(sys.argv[1]);mode=sys.argv[2];build=out/'native-points-build';pair=load_pair(SNAP)
if mode=='build':
 build.mkdir(exist_ok=False)
 flags_path=RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make';link_path=RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt'
 flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in flags_path.read_text().splitlines() if ' = ' in l}
 link=shlex.split(link_path.read_text());obj=build/'native_points.o';binary=build/'native_points'
 commands=[[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(build),*flags['CXX_FLAGS'],'-c',str(ROOT/'tools/linux/stage_b/native_points.cpp'),'-o',str(obj)],
 [str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]]
 for cmd in commands:subprocess.run(cmd,cwd=RETAINED,check=True)
 assert runtime_binding(binary)==pair[0]['runtime_bindings'][pair[0]['stage_binary']]
 files=[ROOT/'tests/native/stage.cpp',ROOT/'tests/diagnostics/checkpoint_replay.cpp',ROOT/'tools/linux/stage_b/native_points.cpp',pathlib.Path(__file__),binary,flags_path,link_path]
 files += [(RETAINED/x).resolve() for x in link if x.endswith(('.a','.so'))]
 save(build/'build.json',dict(commands=commands,files={str(p):digest(p) for p in files}))
else:
 assert mode in ('source','control');s=pair[mode=='control']
 for p,h in json.loads((build/'build.json').read_text())['files'].items():assert digest(p)==h,p
 args=[str(build/'native_points'),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(out/(mode+'-native-points-assembly')),str(out/'native'/mode),str(out/(mode+'-native-points')),mode]
 save(out/(mode+'-native-points-invocation.json'),dict(args=args,build=digest(build/'build.json'),job=os.environ.get('SLURM_JOB_ID')))
 start=time.monotonic();code=subprocess.run(args,cwd=ROOT).returncode
 save(out/(mode+'-native-points-completion.json'),dict(exit_code=code,elapsed_s=time.monotonic()-start));assert code==0
