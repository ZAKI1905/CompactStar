#!/usr/bin/env python3
"""Freeze holdout predictions before any withheld X2 result is generated."""
import csv
import json
import os
import numpy as np
from reference import *

if not os.environ.get('SLURM_JOB_ID'):
    raise RuntimeError('Slurm required')
if any(O.glob('holdout-*-root.json')):
    raise RuntimeError('withheld result already exists')
comp,m1,m3,m4,diag,x,refs,K = load_reference()
knots = x[[0,10,20]]
targets = list(csv.DictReader((RAW/'target_stars.tsv').open(),delimiter='\t'))
B0 = LD(targets[0]['B_solved'])
groups = sorted(n for n in refs if not n.startswith('t_'))
out = dict(schema='compactstar.ba13.frozen-holdouts.v1',interpolator='piecewise-linear-achieved-B',
           predeclaration='93c28c3bd85afea8a893a34a8a6dd92ce06c52be',
           reference_kind='X2',target_indices=[5,15],groups=groups,targets={},
           qualification='UNQUALIFIED',criterion_scope='prospective interpolator consistency only',
           conditions=['curvature continuation','transfer of retained numerical bands',
                       'positive B slope','unseen reference bias and higher-order remainders excluded'],
           source_hashes={str(p.relative_to(R)):sha(p) for p in (
               F/'m1.json',F/'m3.json',F/'m4.json',F/'component-data.json.gz',
               R/'ba13-closure-independent-review-20261020/outputs/reference-diagnostics.json')})
for i in (5,15):
    profile=RAW/f'chemical-profiles/target-{i}/profile.tsv'
    last=profile.read_text().splitlines()[-1].split('\t')
    # Same governed geometric-to-CGS constants used for retained X1/X2 targets.
    ps=float(last[4])*(2.99792458e10**4/(6.673e-8*1e10))
    a,b=(0,10) if i==5 else (10,20)
    hc0,hc1=[LD(json.loads((R/f'x2-fixed-step/t{j}-n4096.json').read_text(),
                          parse_float=str)['values'][8]) for j in (a,b)]
    weight=(x[i]-x[a])/(x[b]-x[a])
    record=dict(B=targets[i]['B_solved'],surfaceP_cgs=format(ps,'.17g'),
                surfacePgeom=last[4],surface_profile_sha256=sha(profile),
                pressure_conversion='Pgeom * c^4/(G*1e10); binary64 physical constants retained',
                H_guess=str(hc0+weight*(hc1-hc0)),achieved_x=float(x[i]),bracket=[a,b],groups={})
    for n in groups:
        prediction,_,_=interpolation([x[i]],knots,refs[n])
        rem=K[n]*float(abs((x[i]-x[a])*(x[i]-x[b])))/2
        bands={t:sum(m1['reference'][str(t)][n][k] for k in ('integration','eos_equilibrium','center'))
               for t in (0,10,20)}
        endpoint=(1-float(weight))*bands[a]+float(weight)*bands[b]
        holdout=max(bands.values())
        slope=float(np.max(abs(np.asarray(m4['common_B'][n]['sensitivities']))))
        btail=max(m1['root_alignment'][str(t)]['functional_GL4_tail_estimate'] for t in (0,10,20))
        root_endpoint=2*slope*max(m1['root_alignment'][str(t)]['computed_B_residual']+btail
                                for t in (a,b))/float(B0)/.999
        root_holdout=slope*(2e-16*float(targets[i]['B_solved'])+btail)/float(B0)/.999
        record['groups'][n]=dict(prediction=list(map(float,prediction[0])),K=K[n],
            interpolation=rem,endpoint_band=endpoint,holdout_band=holdout,
            root_endpoint=root_endpoint,root_holdout=root_holdout,
            criterion=rem+endpoint+holdout+root_endpoint+root_holdout)
    out['targets'][str(i)]=record
print(json.dumps(out,indent=2,allow_nan=False))
