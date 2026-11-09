#!/usr/bin/env python3
"""Build reporting-only collectors against the unchanged retained numerical archives."""
import hashlib,json,pathlib,shlex,subprocess,sys
root=pathlib.Path(__file__).resolve().parents[3]
retained=pathlib.Path(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
sys.path.insert(0,str(root/'tools/linux/native_pair'))
from snapshot import load_pair,digest
from runtime import runtime_binding
snapshots=root/'docs/validation/results-first-bnv-pilot/native-runtime-retry/snapshots'
pair=load_pair(snapshots)
stage=(root/'tests/native/stage.cpp').read_text()
native_main=stage[stage.index('BNV::MainTrajectoryResult native_main('):stage.index('void stage_role(')]
inc=out/'retained_native_main.inc';inc.write_text(native_main)
build=retained/'build/tests';cmake=build/'CMakeFiles/native_p2_stage.dir'
flags={line.split(' = ',1)[0]:shlex.split(line.split(' = ',1)[1]) for line in (cmake/'flags.make').read_text().splitlines() if ' = ' in line}
link=shlex.split((cmake/'link.txt').read_text());compiler=link[0]
commands=[]
for name in ['checkpoint_replay','checkpoint_evidence_test']:
 obj=out/(name+'.o');binary=out/name
 command=[compiler,*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(out),*flags['CXX_FLAGS'],'-c',str(root/'tests/diagnostics'/f'{name}.cpp'),'-o',str(obj)]
 commands.append(command);subprocess.run(command,cwd=build,check=True)
 command=[str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]
 commands.append(command);subprocess.run(command,cwd=build,check=True)
commands.append([str(out/'checkpoint_evidence_test'),str(out/'regression-evidence')]);subprocess.run(commands[-1],check=True)
files=[root/'tests/native/stage.cpp',inc,cmake/'flags.make',cmake/'link.txt',pathlib.Path(compiler)]
files+=list((root/'tests/diagnostics').glob('*'))
files+=list((root/'tools/linux/native_qualification').glob('*diagnostic*'))
files += [p.resolve() for x in link if x.endswith('.a') or x.endswith('.so') for p in [build/x if not pathlib.Path(x).is_absolute() else pathlib.Path(x)]]
files += [out/'checkpoint_replay',out/'checkpoint_evidence_test']
binding=runtime_binding(out/'checkpoint_replay')
assert binding==pair[0]['runtime_bindings'][pair[0]['stage_binary']], 'instrumentation changed runtime linkage'
files += [pathlib.Path(d['resolved']) for d in binding['dependencies']]
manifest=dict(schema='native-checkpoint-observability-build-v1',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),snapshots={r:digest(snapshots/(r+'.json')) for r in ['source','control']},commands=commands,working_directory=str(build),runtime_binding=binding,files=[dict(path=str(p.resolve()),sha256=digest(p)) for p in sorted(set(files))],native_main_verbatim_sha256=digest(inc),numerical_libraries_reused=True)
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('PASS diagnostic collector built; original snapshots and numerical archives unchanged')
