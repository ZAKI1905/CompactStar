from pathlib import Path
import subprocess,json,re,hashlib,gzip
q=Path(__file__).resolve().parent
reports=[]
for mode in ['Debug','Release']:
 b=q/f'new-{mode}';bad=[];links=[];task=[]
 for p in sorted(b.rglob('link.txt')):
  text=p.read_text();links.append(dict(file=str(p.relative_to(b)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),line=text.strip()))
  if re.search(r'libpython|dynamic_lookup|libomp|dependencies/(?:lib|include)|Python3::',text):bad.append(str(p))
 commands=json.loads((b/'compile_commands.json').read_text())
 for c in commands:
  if 'dependencies/include' in c['command'] or '-fopenmp' in c['command']:bad.append(c['file'])
 for f in sorted(b.rglob('*.o.d')):
  if 'dependencies/include' in f.read_text():bad.append(str(f))
 executable=b/'tests/dependency_migration/taskmanager_migration'
 for flag,name in [(['nm','-m'],'nm'),(['otool','-L'],'dylibs'),(['otool','-tvV'],'disassembly')]:
  data=subprocess.check_output([*flag,str(executable)])
  (q/'evidence'/f'final-{mode}-{name}.txt.gz').write_bytes(gzip.compress(data,mtime=0))
  if name=='nm':
   suspicious=[l for l in data.decode().splitlines() if re.search(r'\b_?(?:Py_|_Py|numpy|matplotlib|omp_|GOMP_|__kmpc)',l)]
   if suspicious:bad+=suspicious
 mapdata=(b/'tests/dependency_migration/taskmanager.map').read_bytes()
 (q/'evidence'/f'final-{mode}-taskmanager.map.gz').write_bytes(gzip.compress(mapdata,mtime=0))
 reports.append(dict(mode=mode,failures=bad,link_count=len(links),compile_count=len(commands),taskmanager_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),links=links))
(q/'evidence/final-link-audit.json').write_text(json.dumps(reports,indent=2)+'\n')
for r in reports:print(r['mode'],r['link_count'],'links',r['compile_count'],'TUs',len(r['failures']),'failures')
assert not any(r['failures'] for r in reports)
