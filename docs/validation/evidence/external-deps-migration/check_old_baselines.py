import pathlib,subprocess,hashlib,json,os,time
root=pathlib.Path(__file__).resolve().parent
out=root/'old-Debug-emitted';out.mkdir(exist_ok=False)
src=pathlib.Path('/Users/keeper/Documents/CompactStar/repo/CompactStar')
data='/Users/keeper/Documents/CompactStar/data/compose'
cases=[('tov_reference_cmf',['--emit',str(out/'tov_dscmf1_reference.tsv')],['tov_dscmf1_reference.tsv']),('hartle_moment_inertia_cmf',['--emit',str(out/'hartle_I_dscmf1_debug.tsv')],['hartle_I_dscmf1_debug.tsv']),('tov_path_equivalence_cmf',['--emit',str(out/'tov_path_equivalence_dscmf1.tsv')],['tov_path_equivalence_dscmf1.tsv']),('grid_convergence_cmf',['--emit-dir',str(out)],['grid_convergence_cmf_1p6_debug.tsv','grid_convergence_cmf_1p6_trajectory.tsv'])]
records=[]
for target,args,files in cases:
 cmd=[str(root/'old-Debug/tests'/target),data]+args
 start=time.monotonic()
 with (root/'evidence'/f'old-baseline-{target}.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 record=dict(command=cmd,returncode=r.returncode,wall_seconds=time.monotonic()-start,files=[])
 for name in files:
  expected=hashlib.sha256((src/'tests/baselines'/name).read_bytes()).hexdigest()
  actual=hashlib.sha256((out/name).read_bytes()).hexdigest() if (out/name).exists() else None
  record['files'].append(dict(name=name,expected=expected,actual=actual,equal=expected==actual))
 records.append(record);(root/'evidence/old-baseline-reproduction.json').write_text(json.dumps(records,indent=2)+'\n')
 print(json.dumps(record),flush=True)
 if r.returncode or not all(f['equal'] for f in record['files']):raise SystemExit(1)
