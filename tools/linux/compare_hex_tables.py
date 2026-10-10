#!/usr/bin/env python3
"""Diagnostic comparison of historical hex-float tables; no equivalence verdict."""
import argparse
import json
import math
import re
import struct
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('reference',type=Path);p.add_argument('candidate',type=Path);p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
left=a.reference.read_text().splitlines();right=a.candidate.read_text().splitlines()
result=dict(status='DIAGNOSTIC_ONLY_NO_EQUIVALENCE_AUTHORITY',reference=str(a.reference),candidate=str(a.candidate),
            reference_rows=len(left),candidate_rows=len(right),changed_rows=0,changed_float_fields=0,
            changed_bit_encodings=0,other_changed_fields=0,max_absolute=0.,max_relative=0.,max_ulp=0,samples=[])
def order(x):
    i=struct.unpack('>Q',struct.pack('>d',x))[0]
    return (~i & ((1<<64)-1)) if i>>63 else i|(1<<63)
for i,(x,y) in enumerate(zip(left,right),1):
    if x==y: continue
    result['changed_rows']+=1
    xx=x.split('\t'); yy=y.split('\t')
    if len(xx)!=len(yy): result['other_changed_fields']+=1; continue
    for j,(u,v) in enumerate(zip(xx,yy)):
        if u==v: continue
        if '0x' in u and 'p' in u and '0x' in v and 'p' in v:
            f,g=float.fromhex(u),float.fromhex(v)
            if not math.isfinite(f) or not math.isfinite(g):
                result['other_changed_fields']+=1; continue
            if f==g: continue
            result['changed_float_fields']+=1
            absolute=abs(f-g);relative=absolute/max(abs(f),abs(g)) if f or g else 0
            ulp=abs(order(f)-order(g))
            result['max_absolute']=max(result['max_absolute'],absolute)
            result['max_relative']=max(result['max_relative'],relative)
            result['max_ulp']=max(result['max_ulp'],ulp)
            if len(result['samples'])<8:result['samples'].append(dict(row=i,column=j,reference=u,candidate=v,ulp=ulp))
        elif re.fullmatch('[0-9a-f]{16}',u) and re.fullmatch('[0-9a-f]{16}',v):
            result['changed_bit_encodings']+=1
        else: result['other_changed_fields']+=1
result['unpaired_rows']=abs(len(left)-len(right))
a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='samples'},indent=2))
