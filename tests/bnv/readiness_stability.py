#!/usr/bin/env python3
"""Archived-matrix stability derivation, independent of the nonlinear BNV solve."""
import argparse,csv,gzip,hashlib,io,json,pathlib,subprocess,sys
from fractions import Fraction as F
import numpy as np
import scipy
from scipy.optimize import brentq
import mpmath as mp

REPO=pathlib.Path(__file__).resolve().parents[2]
ARCH=REPO/'docs/validation/evidence/phase6-breakpoint-candidate'
def rows(p): return list(csv.DictReader(io.StringIO(gzip.decompress(p.read_bytes()).decode()),delimiter='\t'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def polynomial(b):
    A=[[F(0) for _ in range(6)] for _ in range(6)]
    for i,r in enumerate([[F(1,4)],[F(3,32),F(9,32)],[F(1932,2197),F(-7200,2197),F(7296,2197)],
                         [F(439,216),F(-8),F(3680,513),F(-845,4104)],
                         [F(-8,27),F(2),F(-3544,2565),F(1859,4104),F(-11,40)]],1): A[i][:len(r)]=r
    v=[F(1)]*6;c=[F(1)]
    for _ in range(6):
        c.append(sum(x*y for x,y in zip(b,v)))
        v=[sum(x*y for x,y in zip(row,v)) for row in A]
    return c
P5=polynomial([F(16,135),F(0),F(6656,12825),F(28561,56430),F(-9,50),F(2,55)])
P4=polynomial([F(25,216),F(0),F(1408,2565),F(2197,4104),F(-1,5),F(0)])
def matrix_poly(c,J,h):
    R=np.zeros_like(J);I=np.eye(3)
    for x in reversed(c):R=R@(h*J)+float(x)*I
    return R

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=pathlib.Path,required=True);ap.add_argument('--probe',type=pathlib.Path,required=True);args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    # Verify every consumed archived file against the accepted manifest.
    manifest=json.loads((ARCH/'manifest.json').read_text())
    names=['runs/supplement01/output/tail-jacobian.tsv.gz','runs/supplement01/output/tail-stability.tsv.gz','runs/final07/output/tail-observations.tsv.gz']
    for n in names:
        rec=next(x for x in manifest['files'] if x['path']==n);assert sha(ARCH/n)==rec['sha256']
    matrices={}
    for r in rows(ARCH/names[0]):
        key=(int(r['endpoint']),int(r['refinement']));matrices.setdefault(key,np.zeros((3,3)))[int(r['row']),int(r['column'])]=float(r['derivative_per_s'])
    direct={}
    for r in rows(ARCH/names[1]):
        key=(int(r['endpoint']),int(r['refinement']),float(r['h_s']));direct.setdefault(key,np.zeros((3,3)))[int(r['row']),int(r['column'])]=float(r['RKF45_amplification'])
    discrepancy=max(np.max(np.abs(matrix_poly(P5,matrices[e,r],h)-M)) for (e,r,h),M in direct.items())
    assert discrepancy<1e-10,'tableau does not match retained direct GSL'
    def radius(c):
        def f(x):return abs(np.polynomial.polynomial.polyval(-x,list(map(float,c))))-1
        grid=np.linspace(.01,10,10000)
        k=next(i for i,x in enumerate(grid) if f(x)>0)
        return brentq(f,grid[k-1],grid[k],xtol=1e-14)
    r5,r4=radius(P5),radius(P4)
    rates=[max(abs(np.linalg.eigvals(J[1:,1:]))) for J in matrices.values()]
    variation=max(np.linalg.norm(matrices[e,1][1:,1:]-matrices[e,0][1:,1:],2) for e in [0,1])
    upper=1.01*(max(rates)+2*variation)
    cap=float(.5*min(r5,r4)/upper)
    eqrow=rows(ARCH/names[2])[0];eq=np.array([float(eqrow[k]) for k in ['x','e','m']])
    scale=np.diag([1.,1e-7,1e-7]);records=[]
    for key,J in matrices.items():
        S=np.linalg.solve(scale,J@scale);eig,V=np.linalg.eig(S)
        for label,h in [('review_scale',4.8e10),('candidate',cap),('half',cap/2),('quarter',cap/4),('old_max',146614328829.19922)]:
            R=matrix_poly(P5,S,h)
            records.append({'endpoint':key[0],'refinement':key[1],'label':label,'h_s':h,'spectral_radius':float(max(abs(np.linalg.eigvals(R)))),'scaled_norm_2':float(np.linalg.norm(R,2)),'eigenvector_condition_2':float(np.linalg.cond(V))})
    J=matrices[0,1];inp=args.output/'linear-input.txt';inp.write_text(' '.join(map(repr,[*map(float,J.flat),*map(float,eq),cap]))+'\n')
    cmd=[str(args.probe.resolve()),str(inp.resolve())]
    with (args.output/'linear.tsv').open('w') as f:subprocess.run(cmd,stdout=f,check=True)
    table=list(csv.DictReader((args.output/'linear.tsv').open(),delimiter='\t'));groups={}
    for r in table:groups.setdefault(tuple(int(r[k]) for k in ['tier','axis','sign','policy']),[]).append(r)
    assert len(groups)==72
    mp.mp.dps=70;M=mp.matrix([[mp.mpf(float(x)) for x in row] for row in J]);E=mp.matrix([mp.mpf(float(x)) for x in eq]);summaries=[]
    # Compare against 70-digit exponential of the exact binary64 J and initial state.
    # Re-evaluate all reference endpoints at 100 digits to bound arithmetic error.
    max_refinement=0.
    for key,group in groups.items():
        tier,axis,sign,policy=key;initial=mp.matrix([mp.mpf(float(group[0][k])) for k in ['x','e','m']])-E
        ratios=[];times=[];ref_errors=[]
        for row in group:
            t=float(row['time']);times.append(t);mp.mp.dps=70;ref=E+mp.expm(M*mp.mpf(t))*initial
            mp.mp.dps=100;fine=E+mp.expm(M*mp.mpf(t))*initial
            y=np.array([float(row[k]) for k in ['x','e','m']]);rtol=10.**(-7-2*tier);atol=np.array([1e-12,1e-18,1e-18])*.01**tier
            for i in range(3):
                # Binary64 rounding of the displayed reference is not used for scoring.
                u=2*abs(fine[i]-ref[i])+mp.mpf('1e-65');D=mp.mpf(float(atol[i]))+mp.mpf(rtol)*max(abs(fine[i]),abs(mp.mpf(float(y[i]))))
                ratios.append(float((abs(mp.mpf(float(y[i]))-fine[i])+u)/D));ref_errors.append(float(u/D))
        assert times[0]==0 and times[-1]==5e11 and all(a<b for a,b in zip(times,times[1:]))
        end=group[-1];acc=len(group)-1;rej=int(end['rejected']);max_refinement=max(max_refinement,max(ref_errors))
        summaries.append({'tier':tier,'axis':axis,'sign':sign,'policy':policy,'accepted':acc,'rejected':rej,'rejection_fraction':rej/(acc+rej),'rhs':int(end['rhs']),'wall_s':float(end['wall_s']),'max_error_over_D':max(ratios),'max_step_s':max(float(r['h_s']) for r in group)})
    result={'status':'DIAGNOSTIC ONLY: archived fixed Jacobian, not nonlinear BNV qualification','polynomial_5':[str(x) for x in P5],'polynomial_4':[str(x) for x in P4],'negative_real_radius_5':r5,'negative_real_radius_4':r4,'chemical_rate_upper_screen_per_s':upper,'jacobian_refinement_norm_per_s':variation,'candidate_cap_s':cap,'direct_GSL_matrix_max_abs_discrepancy':float(discrepancy),'amplifications':records,'linear_perturbations':summaries,'linear_reference_max_uncertainty_over_D':max_refinement,'inputs':{n:sha(ARCH/n) for n in names},'probe_sha256':sha(args.probe),'script_sha256':sha(pathlib.Path(__file__)),'command':cmd,'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__,'python':sys.version,'reference_limit':'high precision bounds exponential arithmetic of rounded frozen matrix only; no nonlinear/Jacobian modeling bound'}
    (args.output/'assessment.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['candidate_cap_s','negative_real_radius_5','negative_real_radius_4','direct_GSL_matrix_max_abs_discrepancy','linear_reference_max_uncertainty_over_D']},indent=2))
if __name__=='__main__':main()
