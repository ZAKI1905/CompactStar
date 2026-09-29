from pathlib import Path
import subprocess,time,json,concurrent.futures
q=Path(__file__).resolve().parent
while not (q/'evidence/final-image-changes.json').exists():
 time.sleep(20)
records=[]
def heat(mode,index):
 cmd=[str(q/f'new-{mode}/tests/heat_capacity_v1')]
 start=time.monotonic();log=q/'evidence'/f'final-{mode}-heat-concurrent-{index}.log'
 with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=f)
 return dict(command=cmd,returncode=r.returncode,wall_seconds=time.monotonic()-start,log=log.name)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 records=list(pool.map(lambda args:heat(*args),[(m,i) for m in ['Debug','Release'] for i in range(2)]))
(q/'evidence/final-concurrent-heat-capacity.json').write_text(json.dumps(records,indent=2)+'\n')
assert all(r['returncode']==0 for r in records)
cmd=['/opt/homebrew/bin/ctest','--test-dir',str(q/'new-Release'),'--output-on-failure','-j','4','-E','^phase5d1_controlled_evolution_regression$']
start=time.monotonic()
with (q/'evidence/final-Release-applicable-suite.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=f)
(q/'evidence/final-Release-applicable-suite.json').write_text(json.dumps(dict(command=cmd,returncode=r.returncode,wall_seconds=time.monotonic()-start,exclusion_reason='Phase-5D1 independently rerun in full after restoring protected install CMake files; retain new-Release-phase5d1-repaired.log.'),indent=2)+'\n')
assert r.returncode==0
print('Final concurrent fixture controls and Release applicable suite PASS',flush=True)
