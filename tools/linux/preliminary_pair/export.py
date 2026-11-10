#!/usr/bin/env python3
import json,pathlib,shlex,subprocess,sys
from run import ROOT,BASE,RETAINED,SNAP,load_pair,digest,save
out=pathlib.Path(sys.argv[1]);build=out/'export-build';build.mkdir(exist_ok=False)
s=load_pair(SNAP)[0]
flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in (RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make').read_text().splitlines() if ' = ' in l}
link=shlex.split((RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt').read_text());obj=build/'export.o';binary=build/'export'
commands=[[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],'-c',str(ROOT/'tools/linux/preliminary_pair/export.cpp'),'-o',str(obj)],
[str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]]
for cmd in commands:subprocess.run(cmd,cwd=RETAINED,check=True)
args=[str(binary),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(out/'export-assembly'),str(out/'thermal-model')]
subprocess.run(args,cwd=ROOT,check=True)
paths=[ROOT/'tools/linux/preliminary_pair/export.cpp',pathlib.Path(__file__),binary]+list((out/'thermal-model').glob('*'))+[pathlib.Path(s['coefficients']),pathlib.Path(s['coefficients']).parent/'tangent.tsv']
save(out/'thermal-model/manifest.json',dict(files={str(p):digest(p) for p in paths},commands=commands,args=args,scope='Frozen binary64 cache cells and evaluated photon/direct-power anchors; static physical uncertainty remains conditional.'))
