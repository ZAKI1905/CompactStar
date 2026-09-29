from pathlib import Path
import subprocess,re,json,gzip,runpy
q=Path(__file__).resolve().parent
# Reuse the same map/body parser as the historical diagnostic, without its run loop.
ns={};exec((q/'audit_provider_probe.py').read_text().split('result={}')[0].replace("r=pathlib.Path(__file__).resolve().parent","r=pathlib.Path.cwd()"),ns)
symbols,bodies=ns['symbols'],ns['bodies'];reports={}
for mode in ['Debug','Release']:
 old=symbols(q/f'provider-{mode}/old.map');new=symbols(q/f'new-{mode}/tests/dependency_migration/taskmanager.map')
 oldasm=bodies(q/f'provider-{mode}/old-disassembly.txt')
 newpath=q/'evidence'/f'final-{mode}-disassembly.txt';newpath.write_bytes(gzip.decompress(newpath.with_suffix('.txt.gz').read_bytes()));newasm=bodies(newpath)
 changes=[]
 for sym in sorted(old.keys() & new.keys()):
  if old[sym]['class']==new[sym]['class']:continue
  item=dict(symbol=sym,old=old[sym],new=new[sym])
  for name,asm in [('old',oldasm),('new',newasm)]:item[name+'_fp']=[line for line in asm.get(sym,[]) if re.search(r'\t(f(?:madd|msub|nmadd|nmsub|add|sub|mul|div|sqrt|cmp|csel)|[su]cvtf|fcvt)',line)]
  changes.append(item)
 coord={name:{sym:{**provider,'arithmetic':[l for l in asm.get(sym,[]) if re.search(r'\tf\w+',l)]} for sym,provider in table.items() if 'Coord3D' in sym and any(x in sym for x in ['XYDist','ltE','eqE'])} for name,table,asm in [('old',old,oldasm),('new',new,newasm)]}
 reports[mode]=dict(common=len(old.keys()&new.keys()),provider_changes=changes,coord3d=coord)
 print(mode,'common',reports[mode]['common'],'changes',len(changes),'FP-bearing',sum(bool(c['old_fp'] or c['new_fp']) for c in changes))
(q/'evidence/final-provider-audit.json').write_text(json.dumps(reports,indent=2)+'\n')
