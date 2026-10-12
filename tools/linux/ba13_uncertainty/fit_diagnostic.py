#!/usr/bin/env python3
"""Unit-rescaling diagnostic using unchanged historical regression on copies."""
import json,math,sys
from pathlib import Path
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')));OLD=S.parents[1]/'ba13-recovery'
sys.path.insert(0,str(S/'tools/linux'));from ba13_refit import load
sys.path.insert(0,str(S/'tests/bnv'));from compare_candidate import regression
x,series=load(OLD/'run1/raw');results={}
for name,(v,e,scales,kind) in series.items():
 items=[]
 for j,scale in enumerate(scales):
  y=[q[j] for q in v[10:]];u=[q[j] for q in e[10:]];base=regression(x[10:],y,u)
  for power in (-40,0,40):
   factor=math.ldexp(1/scale,power);res=regression(x[10:],[z*factor for z in y],[z*factor for z in u])
   cov=res[2][0][0]/factor**2
   items.append(dict(component=j,power=power,intercept_change_scaled=abs(res[0]/factor-base[0])/scale,
    slope_change_scaled=abs(res[1]/factor-base[1])*1e-6/scale,
    intercept_covariance_ratio=cov/base[2][0][0] if base[2][0][0] else None,
    original_errors_below_floor=sum(ee<64*sys.float_info.epsilon*max(1,abs(yy)) for yy,ee in zip(y,u))))
 results[name]=dict(components=len(scales),measurements=items,max_slope_effect_over_window=max(z['slope_change_scaled'] for z in items),
  min_covariance_ratio=min(z['intercept_covariance_ratio'] for z in items if z['intercept_covariance_ratio'] is not None),
  max_covariance_ratio=max(z['intercept_covariance_ratio'] for z in items if z['intercept_covariance_ratio'] is not None))
(R/'fit-units.json').write_text(json.dumps(dict(status='DIAGNOSTIC_NOT_PHYSICAL_UNCERTAINTY',quantities=results),indent=2)+'\n')
for n,v in results.items():print(n,v['min_covariance_ratio'],v['max_covariance_ratio'],v['max_slope_effect_over_window'])
