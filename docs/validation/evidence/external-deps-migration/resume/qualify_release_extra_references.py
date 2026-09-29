from pathlib import Path
import subprocess,json,hashlib,time
q=Path(__file__).resolve().parent; records=[]
eos='/Users/keeper/Documents/CompactStar/data/compose'
for target in ['hartle_monopole_regression','baryon_number_cmf']:
 for stack in ['old','new']:
  path=q/f'{stack}-Release-emitted'/f'{target}.tsv'
  assert not path.exists()
  command=[str(q/f'{stack}-Release/tests'/target),eos,'--emit',str(path)]
  start=time.monotonic()
  with (q/'evidence'/f'{stack}-Release-{target}-emit.log').open('w') as f:p=subprocess.run(command,stdout=f,stderr=f)
  assert p.returncode==0
  records.append({'target':target,'stack':stack,'source':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'wall_seconds':time.monotonic()-start,'command':command})
 old,new=records[-2:];assert old['sha256']==new['sha256']
(q/'evidence/Release-extra-reference-comparison.json').write_text(json.dumps(records,indent=2)+'\n')
print('Both Release references are byte-identical OLD versus NEW')
