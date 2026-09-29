import pathlib,re,json
r=pathlib.Path(__file__).resolve().parent

def symbols(p):
 objects={}; out={}
 for line in p.read_text().splitlines():
  m=re.match(r'\[\s*(\d+)\]\s+(.+)',line)
  if m:objects[m[1]]=m[2]
  m=re.match(r'0x[0-9A-Fa-f]+\s+0x[0-9A-Fa-f]+\s+\[\s*(\d+)\]\s+(.+)',line)
  if m and ('Zaki' in m[2] or 'CONFIND' in m[2]):
   provider=objects[m[1]]
   cls='Zaki' if 'libZaki.a(' in provider else ('CONFIND' if re.search(r'lib(?:CONFIND|Confind)\.a\(',provider) else 'consumer')
   out[m[2]]={'class':cls,'provider':provider}
 return out

def bodies(p):
 out={};name=None
 for line in p.read_text().splitlines():
  if line.endswith(':') and not line.startswith((' ','\t')):name=line[:-1];out[name]=[]
  elif name:out[name].append(line)
 return out
result={}
for mode in ('Debug','Release'):
 p=r/f'provider-{mode}'
 old,new=[symbols(p/f'{a}.map') for a in ['old','new']]
 asm={a:bodies(p/f'{a}-disassembly.txt') for a in ['old','new']}
 changes=[]
 for symbol in sorted(old.keys() & new.keys()):
  if old[symbol]['class']==new[symbol]['class']:continue
  entry={'symbol':symbol,'old':old[symbol],'new':new[symbol]}
  for a in ['old','new']:
   b=asm[a].get(symbol,[])
   entry[a+'_fp']=[l for l in b if re.search(r'\t(f(?:madd|msub|nmadd|nmsub|add|sub|mul|div|sqrt|cmp|csel)|[su]cvtf|fcvt)',l)]
   entry[a+'_instructions']=len(b)
  changes.append(entry)
 result[mode]={'common':len(old.keys() & new.keys()),'changes':changes}
 (r/'evidence/provider-probe-audit.json').write_text(json.dumps(result,indent=2)+'\n')
 print(mode,json.dumps(result[mode],indent=2))
