#!/usr/bin/env python3
import json,pathlib,shlex,subprocess,sys,time,os
from run import ROOT,BASE,RETAINED,SNAP,load_pair,digest,save
out=pathlib.Path(sys.argv[1]);mode=sys.argv[2];build=out/'history-build'
if mode=='build':
 build.mkdir(exist_ok=False)
 flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in (RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make').read_text().splitlines() if ' = ' in l}
 link=shlex.split((RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt').read_text());obj=build/'history_trace.o';binary=build/'history_trace'
 commands=[[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(out/'recovery-build'),*flags['CXX_FLAGS'],'-c',str(ROOT/'tools/linux/preliminary_pair/history_trace.cpp'),'-o',str(obj)],
 [str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]]
 for cmd in commands:subprocess.run(cmd,cwd=RETAINED,check=True)
 save(build/'build.json',dict(commands=commands,files={str(p):digest(p) for p in [ROOT/'tools/linux/preliminary_pair/history_trace.cpp',pathlib.Path(__file__),binary,out/'recovery-build/retained_helpers.inc']}))
else:
 s=load_pair(SNAP)[0 if mode=='source' else 1]
 for path,h in json.loads((build/'build.json').read_text())['files'].items():assert digest(path)==h
 retained=BASE/'stage-a-parallel-20261110-r1/source-stage/run/trajectory/source' if mode=='source' else out/'control'
 args=[str(build/'history_trace'),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(out/f'{mode}-history-assembly'),str(retained),str(out/f'{mode}-history-traces'),mode]
 save(out/f'{mode}-history-invocation.json',dict(args=args,snapshot=digest(SNAP/(mode+'.json'))))
 start=time.monotonic();rc=subprocess.run(args,cwd=ROOT).returncode
 save(out/f'{mode}-history-completion.json',dict(exit_code=rc,elapsed_s=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID')));assert rc==0
