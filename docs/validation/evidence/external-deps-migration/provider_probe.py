import pathlib, subprocess, json, shlex, re, concurrent.futures
r=pathlib.Path(__file__).resolve().parent
src=pathlib.Path('/Users/keeper/Documents/CompactStar/repo/CompactStar')
def run(cmd,log,cwd=None):
 with log.open('w') as f:res=subprocess.run(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT)
 if res.returncode:raise RuntimeError(str(log))
for mode in ('Debug','Release'):
 build=r/f'old-{mode}';probe=r/f'provider-{mode}';probe.mkdir()
 flags=['-g'] if mode=='Debug' else ['-O3','-DNDEBUG']
 common=['/usr/bin/clang++','-std=c++17','-arch','arm64','-pthread',*flags]
 driver=probe/'driver.o'
 run([*common,'-I',str(src),'-I',str(build/'generated/include'),'-isystem',str(src/'dependencies/include'),'-isystem','/opt/local/include','-c',str(r/'taskmanager_driver.cpp'),'-o',str(driver)],probe/'driver.log')
 oldlink=shlex.split((build/'tests/CMakeFiles/compactstar_library_smoke.dir/link.txt').read_text())
 oldlink=[str(driver) if x.endswith('compactstar_library_smoke.cpp.o') else x for x in oldlink]
 oldlink[oldlink.index('-o')+1]=str(probe/'old-taskmanager')
 oldlink+=['-Wl,-map,'+str(probe/'old.map')]
 run(oldlink,probe/'old-link.log',build/'tests')
 (probe/'old-link.json').write_text(json.dumps(oldlink,indent=2)+'\n')
 members=set(re.findall(r'libCompactStar\.a\(([^)]+)\)',(probe/'old.map').read_text()))
 commands=json.loads((build/'compile_commands.json').read_text())
 records=[]
 for c in commands:
  if pathlib.Path(c['file']).name+'.o' not in members:continue
  tokens=shlex.split(c['command']);rel=pathlib.Path(c['file']).relative_to(src);dest=probe/(pathlib.Path(c['file']).name+'.o'); args=[tokens[0]];i=1
  while i<len(tokens):
   t=tokens[i]
   if t in ['-isystem','-I','-o','-c']:i+=2;continue
   if t.startswith('-I'):i+=1;continue
   args.append(t);i+=1
  args+=['-I',str(r/'provider-probe-source'),'-I',str(build/'generated/include'),'-isystem',str(r/mode/'Zaki/include'),'-isystem',str(r/mode/'CONFIND/include'),'-isystem','/opt/local/include','-c',str(r/'provider-probe-source'/rel),'-o',str(dest)]
  records.append({'source':str(rel),'object':str(dest),'command':args})
 assert len(records)==len(members),(len(records),len(members))
 (probe/'new-compile.json').write_text(json.dumps(records,indent=2)+'\n')
 def compileone(c):run(c['command'],pathlib.Path(c['object']+'.log'))
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:list(pool.map(compileone,records))
 # Preserve original archive member ordering for selected members.
 order=subprocess.check_output(['ar','-t',str(build/'libCompactStar.a')],text=True).splitlines()
 objects={pathlib.Path(c['object']).name:c['object'] for c in records}
 run(['ar','rcs',str(probe/'libCompactStar.a'),*[objects[x] for x in order if x in objects]],probe/'ar.log')
 newdriver=probe/'new-driver.o'
 run([*common,'-I',str(r/'provider-probe-source'),'-I',str(build/'generated/include'),'-isystem',str(r/mode/'Zaki/include'),'-isystem',str(r/mode/'CONFIND/include'),'-isystem','/opt/local/include','-c',str(r/'taskmanager_driver.cpp'),'-o',str(newdriver)],probe/'new-driver.log')
 newlink=[*common,str(newdriver),str(probe/'libCompactStar.a'),str(r/mode/'Zaki/lib/libZaki.a'),str(r/mode/'CONFIND/lib/libCONFIND.a'),'/opt/local/lib/libgsl.dylib','/opt/local/lib/libgslcblas.dylib','/opt/local/lib/libomp/libomp.dylib','-lz','-Wl,-map,'+str(probe/'new.map'),'-o',str(probe/'new-taskmanager')]
 run(newlink,probe/'new-link.log')
 (probe/'new-link.json').write_text(json.dumps(newlink,indent=2)+'\n')
 for arm in ['old','new']:
  run(['nm','-m',str(probe/f'{arm}-taskmanager')],probe/f'{arm}-nm.txt')
  run(['otool','-tvV',str(probe/f'{arm}-taskmanager')],probe/f'{arm}-disassembly.txt')
 print(mode,len(members),'linked and disassembled',flush=True)
