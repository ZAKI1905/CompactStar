from pathlib import Path
import subprocess,json,shlex,os
q=Path(__file__).resolve().parent
src=Path('/Users/keeper/Documents/CompactStar/repo/CompactStar'); test=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration/tests/dependency_migration')
for mode,short in [('Debug','D'),('Release','R')]:
 out=q/f'old-plot-free-{mode}';out.mkdir(exist_ok=True);objs=[]
 commands=json.loads((q/f'old-{mode}/compile_commands.json').read_text())
 for rel in ['CompactStar/Core/src/TaskManager.cpp','CompactStar/Extensions/MixedStar/src/DarkCore_Analysis.cpp']:
  c=next(c for c in commands if c['file']==str(src/rel));cmd=shlex.split(c['command']);obj=out/(Path(rel).name+'.o');objs.append(str(obj));cmd[cmd.index('-o')+1]=str(obj);cmd[cmd.index('-c')+1]=str(q/'provider-probe-source'/rel)
  with (out/(obj.name+'.log')).open('w') as f:subprocess.run(cmd,cwd=c['directory'],stdout=f,stderr=f,check=True)
 cmd=json.loads((q/f'provider-{mode}/old-link.json').read_text());cmd[cmd.index('-o')+1]=str(out/'taskmanager');cmd=[x for x in cmd if not x.startswith('-Wl,-map,')];cmd[1:1]=objs
 with (out/'link.log').open('w') as f:subprocess.run(cmd,cwd=q/f'old-{mode}/tests',stdout=f,stderr=f,check=True)
 for fixture in ['T1','T2']:
  subprocess.run(['python3',test/'check_fixture.py','--exe',out/'taskmanager','--mode',fixture,'--build-mode',mode,'--eos','/Users/keeper/Documents/CompactStar/data/compose/DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos','--report',out/f'{fixture}-result.json'],check=True)
