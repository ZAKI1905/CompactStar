#!/usr/bin/env python3
"""Build local reconstruction study against unchanged authenticated archives."""
import json,pathlib,shlex,subprocess,sys
root=pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(root/'tools/linux/native_pair'))
from snapshot import load_pair,digest
from runtime import runtime_binding
retained=pathlib.Path(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
snapshots=root/'docs/validation/results-first-bnv-pilot/native-runtime-retry/snapshots'
pair=load_pair(snapshots)
replay=(root/'tests/diagnostics/checkpoint_replay.cpp').read_text()
inc=out/'retained_trace_helpers.inc'
inc.write_text(replay[replay.index('struct TraceContext'):replay.index('int main(int argc')])
build=retained/'build/tests';cmake=build/'CMakeFiles/native_p2_stage.dir'
flags={line.split(' = ',1)[0]:shlex.split(line.split(' = ',1)[1]) for line in (cmake/'flags.make').read_text().splitlines() if ' = ' in line}
link=shlex.split((cmake/'link.txt').read_text());compiler=link[0]
name='checkpoint_refinement';obj=out/(name+'.o');binary=out/name
commands=[[compiler,*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],'-I'+str(out),*flags['CXX_FLAGS'],'-c',str(root/'tests/diagnostics'/f'{name}.cpp'),'-o',str(obj)]]
commands.append([str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link])
for command in commands:subprocess.run(command,cwd=build,check=True)
files=[root/'tests/native/stage.cpp',inc,cmake/'flags.make',cmake/'link.txt',pathlib.Path(compiler),binary]
files+=list((root/'tests/diagnostics').glob('*'))
files+=[root/'tools/linux/native_qualification'/p for p in ['build_refinement.py','run_refinement.py','refinement.sbatch']]
files+=[root/'docs/validation/results-first-bnv-pilot/echem-refinement/PREDECLARATION.md']
files += [p.resolve() for x in link if x.endswith('.a') or x.endswith('.so') for p in [build/x if not pathlib.Path(x).is_absolute() else pathlib.Path(x)]]
binding=runtime_binding(binary)
assert binding==pair[0]['runtime_bindings'][pair[0]['stage_binary']], 'study changed runtime linkage'
files += [pathlib.Path(d['resolved']) for d in binding['dependencies']]
manifest=dict(schema='echem-local-feasibility-build-v1',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),snapshots={r:digest(snapshots/(r+'.json')) for r in ['source','control']},commands=commands,working_directory=str(build),runtime_binding=binding,files=[dict(path=str(p.resolve()),sha256=digest(p)) for p in sorted(set(files))],numerical_libraries_reused=True)
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('PASS reconstruction study built; no production changes')
