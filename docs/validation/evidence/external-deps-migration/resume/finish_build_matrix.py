from pathlib import Path
import subprocess,time,json,hashlib,sys
q=Path(__file__).resolve().parent
src=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration')
logs=[q/'evidence'/f'{stack}-{mode}-full-suite.log' for stack in ['old','new'] for mode in ['Debug','Release']]
while not all('Total Test time (real)' in p.read_text() for p in logs) or not (q/'adr0017-new-Debug-attempt2/result.json').exists():
 time.sleep(20)
assert json.loads((q/'adr0017-new-Debug-attempt2/result.json').read_text())['pass']
records=[]
def run(cmd,log):
 start=time.monotonic()
 with (q/'evidence'/log).open('w') as f:r=subprocess.run([str(x) for x in cmd],stdout=f,stderr=f,cwd=src)
 records.append({'command':[str(x) for x in cmd],'returncode':r.returncode,'wall_seconds':time.monotonic()-start,'log':log})
 (q/'evidence/final-build-checks.json').write_text(json.dumps(records,indent=2)+'\n')
 print(log,r.returncode,flush=True)
 if r.returncode:raise SystemExit(r.returncode)
for mode in ['Debug','Release']:
 build=q/f'new-{mode}'
 run(['/opt/homebrew/bin/cmake','-S',src,'-B',build],f'final-{mode}-configure.log')
 run(['/opt/homebrew/bin/cmake','--build',build,'--parallel','4'],f'final-{mode}-build.log')
 run(['/opt/homebrew/bin/cmake','--install',build,'--prefix',q/f'install-{mode}'],f'final-{mode}-install.log')
 run(['/opt/homebrew/bin/ctest','--test-dir',build,'--output-on-failure','-j','2','-R','taskmanager_|unavailable_spin|smoke|tov_|hartle_|grid_convergence|rotochemical_local'],f'final-{mode}-affected-tests.log')
 (q/'evidence'/f'final-{mode}-test-inventory.json').write_bytes(subprocess.check_output(['/opt/homebrew/bin/ctest','--test-dir',build,'--show-only=json-v1']))
run([sys.executable,q/'audit_images.py','final'],'final-all-images-console.log')
run([sys.executable,q/'audit_final_links.py'],'final-links-console.log')
run([sys.executable,q/'audit_final_providers.py'],'final-providers-console.log')
old=json.loads((q/'evidence/tested-all-images.json').read_text());new=json.loads((q/'evidence/final-all-images.json').read_text())
changed={mode:[{'file':a['file'],'old':a['sha256'],'new':b['sha256']} for a,b in zip(old[mode],new[mode]) if a['sha256']!=b['sha256']] for mode in old}
(q/'evidence/final-image-changes.json').write_text(json.dumps(changed,indent=2)+'\n')
print('FINAL MATRIX BUILD CHECKS COMPLETE',flush=True)
