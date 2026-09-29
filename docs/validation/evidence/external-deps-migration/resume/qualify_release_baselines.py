from pathlib import Path
import subprocess,hashlib,json,time
q=Path(__file__).resolve().parent;data='/Users/keeper/Documents/CompactStar/data/compose';records=[]
cases=[('tov_reference_cmf','tov_dscmf1_reference.tsv'),('hartle_moment_inertia_cmf','hartle_I_dscmf1_debug.tsv'),('tov_path_equivalence_cmf','tov_path_equivalence_dscmf1.tsv'),('grid_convergence_cmf',None)]
for stack in ['old','new']:
 out=q/f'{stack}-Release-emitted';out.mkdir()
 for target,name in cases:
  cmd=[str(q/f'{stack}-Release/tests'/target),data]+(['--emit',str(out/name)] if name else ['--emit-dir',str(out)])
  start=time.monotonic()
  with (q/'evidence'/f'{stack}-Release-baseline-{target}.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=f)
  records.append(dict(stack=stack,target=target,command=cmd,returncode=p.returncode,wall_seconds=time.monotonic()-start))
  assert p.returncode==0
old={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (q/'old-Release-emitted').glob('*.tsv')}
new={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (q/'new-Release-emitted').glob('*.tsv')}
result=dict(runs=records,old=old,new=new,differences=[p for p in old.keys()|new.keys() if old.get(p)!=new.get(p)])
(q/'evidence/Release-baseline-comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(len(old),len(result['differences']))
assert old==new and len(old)==5
