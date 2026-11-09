#!/usr/bin/env python3
"""Read-only checks plus private temporary corruption; never launches a numerical executable."""
import copy,json,pathlib,sys,tempfile
from snapshot import load_pair,digest,encoded
original=pathlib.Path(sys.argv[1]);load_pair(original)
def refused(directory,label):
 try:load_pair(directory)
 except (ValueError,FileNotFoundError):print('PASS refusal:',label)
 else:raise AssertionError(label)
with tempfile.TemporaryDirectory(prefix='native-snapshot-') as temp:
 d=pathlib.Path(temp)
 for p in original.glob('*.json'):(d/p.name).write_bytes(p.read_bytes())
 s=json.loads((d/'source.json').read_text());s['initial_state']=[1,0,0]
 (d/'source.json').write_bytes(encoded(s));refused(d,'tampered snapshot identity')
 for name,value in [('missing dependency',{'path':str(d/'absent-certificate'),'sha256':'0'*64}),('mismatched dependency',{'path':str(d/'identities.json'),'sha256':'0'*64})]:
  for role in ('source','control'):
   obj=json.loads((original/(role+'.json')).read_text());obj['files'].append(value);(d/(role+'.json')).write_bytes(encoded(obj))
  (d/'identities.json').write_bytes(encoded({r:digest(d/(r+'.json')) for r in ('source','control')}));refused(d,name)
 for role in ('source','control'):
  obj=json.loads((original/(role+'.json')).read_text());obj['initial_state']=[0,1,0];(d/(role+'.json')).write_bytes(encoded(obj))
 (d/'identities.json').write_bytes(encoded({r:digest(d/(r+'.json')) for r in ('source','control')}));refused(d,'invalid serialized initial state')
print('PASS snapshot authentication and refusal tests; no numerical execution')
