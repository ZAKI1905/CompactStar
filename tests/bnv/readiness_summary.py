#!/usr/bin/env python3
"""Summarize immutable bounded evidence without promoting it to a fresh run."""
import csv,gzip,hashlib,io,json,pathlib
import numpy as np
from scipy.linalg import expm
ROOT=pathlib.Path(__file__).resolve().parents[2]
ARCH=ROOT/'docs/validation/evidence/phase6-breakpoint-candidate'
OUT=ROOT/'docs/validation/evidence/phase6-campaign-readiness'
def rows(p):return list(csv.DictReader(io.StringIO(gzip.decompress(p.read_bytes()).decode()),delimiter='\t'))
def allocation(row):
    U=float(row['U_O'])+float(row['U_event']);F=float(row['F_i'])
    if F==0:
        assert U==0, 'nonzero error with zero allocation'
        return 0. # stored exact endpoint, no reconstruction
    return U/(.2*F)
def main():
    manifest=json.loads((ARCH/'manifest.json').read_text())
    def read(name):
        rec=next(r for r in manifest['files'] if r['path']==name);p=ARCH/name
        assert hashlib.sha256(p.read_bytes()).hexdigest()==rec['sha256']
        return rows(p)
    result={'status':'AUTHENTICATED RETAINED EVIDENCE, not a new BNV run','windows':{}}
    for window in ['source','control','bounded','tail']:
        base='runs/final07/output/'+window+'-'
        obs=read(base+'observations.tsv.gz');stats=read(base+'statistics.tsv.gz');unc=read(base+'uncertainty.tsv.gz');r20=read(base+'r20.tsv.gz');ledger=read(base+'ledger.tsv.gz')
        tiers=[]
        for k,name in enumerate(['BASELINE','REFINED','ULTRA']):
            o=[r for r in obs if r['tier']==str(k)];u=[r for r in unc if r['tier']==str(k)];s=stats[k]
            tiers.append({'tier':name,'component_max_accuracy_ratios':{c:max(float(r['accuracy_'+c]) for r in o) for c in ['x','e','m']},'all_reconstructions_qualified':all(r['reconstruction_qualified']=='1' for r in o),'max_reconstruction_allocation_ratio':max(allocation(r) for r in u),'r20':r20[k],'statistics':s})
        result['windows'][window]={'tiers':tiers,'ledger_combined_refined_ultra':ledger,'limitation':'ledger max is aggregated over REFINED/ULTRA; BASELINE independent ledger accuracy was not scored by this harness'}
    from readiness_stability import P5,matrix_poly
    J=np.zeros((3,3))
    for r in read('runs/supplement01/output/tail-jacobian.tsv.gz'):
        if r['endpoint']=='0' and r['refinement']=='1':J[int(r['row']),int(r['column'])]=float(r['derivative_per_s'])
    cap=json.loads((OUT/'linear-stability/assessment.json').read_text())['candidate_cap_s']
    S=np.diag([1.,1e-7,1e-7]);A=np.linalg.solve(S,J@S);prop=[]
    for factor in [1.,.5,.25]:
        h=cap*factor;exact=expm(h*A);R=matrix_poly(P5,A,h)
        prop.append({'cap_factor':factor,'h_s':h,'numerical_norm2':float(np.linalg.norm(R,2)),'exact_flow_norm2':float(np.linalg.norm(exact,2)),'operator_error_norm2':float(np.linalg.norm(R-exact,2))})
    result['nonnormality']={'scales':[1.,1e-7,1e-7],'linear_flow_comparison':prop,'interpretation':'spectral radius below one does not imply componentwise contraction; much of scaled norm growth is present in exact frozen flow. Must transport coupled uncertainty.'}
    seconds=float(json.loads(gzip.decompress((ARCH/'runs/final07/execution.json.gz').read_bytes()))['seconds'])
    observations=3*(65+65+241+33);campaign_observations=6*8193
    result['cost_scenario']={'historical_final07_seconds':seconds,'bounded_scored_observations':observations,'campaign_scored_observations':campaign_observations,'all_in_work_proxy_s_per_observation':seconds/observations,'six_trajectory_proxy_hours':seconds/observations*campaign_observations/3600,'steps_per_500kyr_at_cap_minimum':int(np.ceil(500000*31557600/cap)),'limitation':'Planning proxy, not measured campaign runtime or upper bound. final07 also included assembly, references, envelope tests and crossings; warm windows differ from full paths. Global reference and nested R20 work may dominate; no qualified runtime ceiling exists.'}
    result['depletion']={'fixture_rate_per_year':-1e-12,'duration_year':500000,'absolute_fraction':5e-7,'certificate_ceiling':1e-6,'utilization':.5,'status':'ARITHMETIC_ONLY; unavailable certificate.tsv prevents fresh 21-star drift/certificate verification'}
    (OUT/'retained-summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'cost':result['cost_scenario'],'nonnormality':result['nonnormality']},indent=2))
if __name__=='__main__':main()
