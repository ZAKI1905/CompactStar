#!/usr/bin/env python3
"""Compile retained-reference accepted-endpoint observations; no source rerun."""
import pathlib,shlex,subprocess,sys,time,os
from run import ROOT,BASE,RETAINED,SNAP,load_pair,digest,save
from high_precision import SOURCE
out=pathlib.Path(sys.argv[1]);build=out/'accepted-build';build.mkdir(exist_ok=False)
s=load_pair(SNAP)[0]
flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in (RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make').read_text().splitlines() if ' = ' in l}
link=shlex.split((RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt').read_text());obj=build/'accepted.o';binary=build/'accepted'
commands=[[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(out/'recovery-build'),*flags['CXX_FLAGS'],'-c',str(ROOT/'tools/linux/preliminary_pair/accepted.cpp'),'-o',str(obj)],
[str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]]
for cmd in commands:subprocess.run(cmd,cwd=RETAINED,check=True)
args=[str(binary),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(out/'accepted-assembly'),str(SOURCE),str(out/'control'),str(out/'accepted-observations')]
start=time.monotonic();rc=subprocess.run(args,cwd=ROOT).returncode
save(out/'accepted-completion.json',dict(exit_code=rc,elapsed_s=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID'),commands=commands,args=args,files={str(p):digest(p) for p in [ROOT/'tools/linux/preliminary_pair/accepted.cpp',pathlib.Path(__file__),binary]}));assert rc==0
