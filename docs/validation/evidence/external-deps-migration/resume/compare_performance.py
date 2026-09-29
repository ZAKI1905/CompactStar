from pathlib import Path
import subprocess, json, sys
q=Path(__file__).resolve().parent
src=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration')
reports=[]
for mode in ['Debug','Release']:
    for case in ['T1','T2']:
        for stack in ['old','new']:
            exe=q/f'provider-{mode}/old-taskmanager' if stack=='old' else q/f'new-{mode}/tests/dependency_migration/taskmanager_migration'
            output=q/'evidence'/f'performance-{stack}-{mode}-{case}.json'
            command=[sys.executable,str(src/'tests/dependency_migration/check_fixture.py'),'--exe',str(exe),'--build-mode',mode,'--mode',case,'--eos','/Users/keeper/Documents/CompactStar/data/compose/DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos','--report',str(output)]
            subprocess.run(command,check=True)
            reports.append({'stack':stack,**json.loads(output.read_text()),'command':command})
(q/'evidence/paired-performance.json').write_text(json.dumps({'load_note':'Back-to-back OLD/NEW measurements; existing long governed tests continue on this Mac. Not an idle-machine benchmark.','runs':reports},indent=2)+'\n')
