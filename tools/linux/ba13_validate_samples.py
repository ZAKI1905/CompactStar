#!/usr/bin/env python3
"""Audit complete static BA13 samples and exact native repeat, no trajectory."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True


def rows(path,count):
    with path.open() as f:r=list(csv.DictReader(f,delimiter='\t'))
    if len(r)!=count:raise ValueError(f'{path}: incomplete rows')
    return [{k:float(v) for k,v in item.items()} for item in r]


def main():
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--repeat',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    files={'target_stars.tsv':21,'tangent.tsv':21,'coefficients.tsv':21,'thermal_surface.tsv':273,'profile_metrics.tsv':43029}
    data={n:rows(a.raw/n,count) for n,count in files.items()}
    errors=[];identities={}
    for n in files:
        left=(a.raw/n).read_bytes();right=(a.repeat/n).read_bytes()
        identities[n]=dict(sha256=hashlib.sha256(left).hexdigest(),repeat_sha256=hashlib.sha256(right).hexdigest(),exact=left==right)
        if left!=right:errors.append('repeat mismatch '+n)
        if not all(math.isfinite(v) for row in data[n] for v in row.values()):errors.append('nonfinite '+n)
    targets=data['target_stars.tsv'];tangents=data['tangent.tsv'];coefs=data['coefficients.tsv'];thermal=data['thermal_surface.tsv'];profile=data['profile_metrics.tsv'];B0=targets[0]['B_solved']
    for i,row in enumerate(targets):
        if row['index']!=i or row['target_fraction'] != -5e-8*i:errors.append('nominal target '+str(i))
        if abs(row['residual'])>5e-11*B0 or not 0<=row['bracket_width']<=5e-11*B0:errors.append('target tolerance '+str(i))
        if abs((row['B_solved']-row['B_target'])-row['residual'])>64*sys.float_info.epsilon*B0:errors.append('target consistency '+str(i))
        if tangents[i]['B_solved']!=row['B_solved'] or coefs[i]['B_solved']!=row['B_solved']:errors.append('owner B mismatch '+str(i))
        if abs(sum(tangents[i][n] for n in ('t_n','t_e','t_mu'))-1)>2*sys.float_info.epsilon:errors.append('closed tangent '+str(i))
        if abs(tangents[i]['raw_closure'])>tangents[i]['closure_budget']:errors.append('raw tangent closure '+str(i))
        for n in ('e','mu'):
            if tangents[i]['t_'+n]<=tangents[i]['u_t_'+n]:errors.append('lepton lower bound '+str(i))
        for k,t in enumerate(thermal[i*13:(i+1)*13]):
            expected=10**(6+.25*k)
            if t['index']!=i or abs(t['Tinf_K']-expected)>4*sys.float_info.epsilon*expected:errors.append('thermal grid '+str(i))
            if t['Cstar_erg_K']<=0 or t['Tsurface_inf_K']<=0:errors.append('thermal positivity '+str(i))
        if coefs[i]['Ltilde_Me']<=0 or coefs[i]['Ltilde_Mmu']<=0:errors.append('enabled channels '+str(i))
        if coefs[i]['P2_inf']!=0 or coefs[i]['u_P2']!=0:errors.append('exact P2 structural zero '+str(i))
        for j,rowp in enumerate(profile[i*2049:(i+1)*2049]):
            ref=profile[j]
            if rowp['index']!=i or rowp['q']!=ref['q'] or rowp['active_dimension']!=ref['active_dimension']:errors.append('sampled source support '+str(i));break
    profile_files={}
    for pth in sorted((a.raw/'chemical-profiles').rglob('*')):
        if pth.is_file():
            rel=pth.relative_to(a.raw);other=a.repeat/rel
            same=other.is_file() and pth.read_bytes()==other.read_bytes()
            profile_files[str(rel)]=same
            if not same:errors.append('profile/adapter repeat '+str(rel))
    result=dict(status='PASS_STATIC_SAMPLE_INTEGRITY' if not errors else 'FAIL',errors=errors,
        sample_count=21,quantities_expected=22,thermal_knots=13,profile_nodes=2049,
        max_target_residual_over_B0=max(abs(r['residual'])/B0 for r in targets),
        max_bracket_width_over_B0=max(r['bracket_width']/B0 for r in targets),
        target_bound=5e-11,fit_indices=list(range(10,21)),
        achieved_B_strictly_decreasing=all(a['B_solved']>b['B_solved'] for a,b in zip(targets,targets[1:])),
        files=identities,profile_adapter_files=len(profile_files),exact_profiles=all(profile_files.values()),
        frozen_certificate_qualified=False,
        limitation='Input integrity and sampled support only; final historical fit and complete uncertainty audit are separate')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(result['status'])
    if errors:raise SystemExit(1)

if __name__=='__main__':main()
