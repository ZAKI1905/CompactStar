from pathlib import Path
import subprocess,json,gzip,hashlib
q=Path(__file__).resolve().parent
src=Path('/Users/keeper/Documents/CompactStar/worktrees/ZakiLib-2.0.1-fp-preservation/tests/fp_preservation')
records=[]
for mode in ['Debug','Release']:
 exe=q/f'characterize-{mode}'
 cmd=['/usr/bin/clang++','-std=c++17','-O0','-ffp-contract=off','-arch','arm64','-I',str(q/mode/'Zaki/include'),'-isystem','/opt/local/include',str(src/'characterize.cpp'),str(q/f'new-{mode}/libCompactStar.a'),str(q/mode/'Zaki/lib/libZaki.a'),'/opt/local/lib/libgsl.dylib','/opt/local/lib/libgslcblas.dylib','-lz','-o',str(exe)]
 subprocess.run(cmd,check=True)
 for stress in [False,True]:
  dest=q/'evidence'/f'characterize-{mode}{"-stress" if stress else ""}.tsv'
  args=[str(exe),str(dest),'stress' if stress else 'fixed','/Users/keeper/Documents/CompactStar/data/compose/DS-CMF-1-with-crust','DS(CMF)-1_with_crust.eos']
  subprocess.run(args,check=True)
  oracle=gzip.decompress((src/'fixtures'/('oracle-stress.tsv.gz' if stress else 'oracle.tsv.gz')).read_bytes())
  lines=dest.read_bytes().splitlines();old=oracle.splitlines();diff=sum(a!=b for a,b in zip(lines,old))+abs(len(lines)-len(old))
  records.append(dict(mode=mode,stress=stress,records=len(lines),differences=diff,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),command=cmd))
  (q/'evidence/integrated-zaki-characterization.json').write_text(json.dumps(records,indent=2)+'\n')
  print(mode,stress,len(lines),diff,flush=True)
  assert diff==0
