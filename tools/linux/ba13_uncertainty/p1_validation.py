#!/usr/bin/env python3
import json,sys
from pathlib import Path
import numpy as np
import mpmath as mp
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')));records=[]
for i in (0,10,20):
 d=np.loadtxt(R/'runs'/f'nominal-{i}'/'profile/profile.tsv',skiprows=1)
 for j in np.linspace(0,len(d)-1,9,dtype=int):
  pf=197.3269804*np.cbrt(3*np.pi**2*d[j,6]);answers=[]
  for prec in (80,112):
   with mp.workdps(prec):answers.append(float(mp.quad(lambda x:3*x*x*mp.sqrt(mp.mpf(939.56542052)**2+(mp.mpf(float(pf))*x)**2),[0,1])))
  assert answers[0]==answers[1]
  vals=[]
  for n in (12,24):
   x,w=np.polynomial.legendre.leggauss(n);t=(x+1)/2;vals.append(float(np.sum(w/2*3*t*t*np.hypot(939.56542052,pf*t))))
  error=max(abs(v-answers[0])/answers[0] for v in vals);assert error<128*sys.float_info.epsilon
  records.append(dict(index=i,node=int(j),pf=float(pf),relative_error=error,status='PASS_FIXED_INPUT_P1_ARITHMETIC'))
(R/'p1-validation.json').write_text(json.dumps(records,indent=2)+'\n');print('checks',len(records),'maximum relative error',max(v['relative_error'] for v in records))
