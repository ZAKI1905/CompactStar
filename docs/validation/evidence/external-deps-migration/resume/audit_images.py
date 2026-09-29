from pathlib import Path
import subprocess,shlex,hashlib,json,re,sys
q=Path(__file__).resolve().parent
label=sys.argv[1] if len(sys.argv)>1 else 'final'
report={}
for mode in ['Debug','Release']:
 b=q/f'new-{mode}';images=[]
 for f in sorted(b.rglob('link.txt')):
  if 'phase5d1-regression-evidence' in f.parts:continue
  words=shlex.split(f.read_text())
  if '-o' not in words:continue
  p=(f.parents[2]/words[words.index('-o')+1]).resolve()
  if p not in images:images.append(p)
 items=[]
 for p in images:
  output=subprocess.check_output(['/usr/bin/nm','-m',str(p)],text=True)
  forbidden=[l for l in output.splitlines() if re.search(r'\b_?(?:Py_|_Py|numpy|matplotlib|omp_|GOMP_|__kmpc|gROOT)|_ZN[0-9]+(?:ROOT|TROOT|TCanvas|TFile|TTree)',l)]
  dylibs=subprocess.check_output(['/usr/bin/otool','-L',str(p)],text=True)
  if re.search(r'libpython|libomp|libgomp|libiomp|lib(?:ROOT|Core|RIO|Hist|Graf|Gpad|Tree)\.',dylibs):forbidden.append(dylibs)
  items.append({'file':str(p.relative_to(b)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'forbidden_symbols_or_dylibs':forbidden,'dylibs':dylibs.splitlines()[1:]})
 report[mode]=items
 assert not any(x['forbidden_symbols_or_dylibs'] for x in items)
 print(mode,len(items),'images clean',flush=True)
(q/'evidence'/f'{label}-all-images.json').write_text(json.dumps(report,indent=2)+'\n')
