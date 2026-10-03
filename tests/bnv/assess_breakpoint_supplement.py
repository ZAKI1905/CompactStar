#!/usr/bin/env python3
"""Common-time run4 comparison and bounded tail stability evidence; no ratification."""
import argparse,csv,json,math,os,pathlib
import numpy as np

def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def assess(root):
    result={'status':'PROPOSED successor evidence, not campaign qualification','run4':{},'tail_stability':[]}
    for name,n in [('source',490),('control',1277)]:
        r=rows(root/(name+'-run4-comparison.tsv'))
        assert len(r)==3*n
        for component in range(3):
            g=[x for x in r if x['component']==str(component)]
            assert len(g)==n
            times=[float(x['time']) for x in g]
            assert all(math.isfinite(t) for t in times) and all(a<b for a,b in zip(times,times[1:]))
        values={k:max(float(x[k]) for x in r) for k in ['archive_difference_over_DU','new_reference_uncertainty_over_DU','BASELINE_accuracy','REFINED_accuracy','ULTRA_accuracy']}
        values['observations']=n
        values['pass']=values['new_reference_uncertainty_over_DU']<=.05 and all(values[k]<=1 for k in ['BASELINE_accuracy','REFINED_accuracy','ULTRA_accuracy'])
        # Archived reference discrepancy is reported, not assigned an invented uncertainty.
        result['run4'][name]=values
    jac=rows(root/'tail-jacobian.tsv');amp=rows(root/'tail-stability.tsv')
    assert len(jac)==36 and len(amp)==108
    for endpoint in range(2):
        matrices=[]
        for refinement in range(2):
            J=np.zeros((3,3))
            for x in jac:
                if x['endpoint']==str(endpoint) and x['refinement']==str(refinement):J[int(x['row']),int(x['column'])]=float(x['derivative_per_s'])
            matrices.append(J)
            modes=np.linalg.eigvals(J)
            record={'endpoint':endpoint,'refinement':refinement,'eigenvalues_per_s':[[v.real,v.imag] for v in modes],'amplification':[]}
            for factor in [.5,1.,2.]:
                group=[x for x in amp if x['endpoint']==str(endpoint) and x['refinement']==str(refinement) and float(x['h_factor'])==factor]
                assert len(group)==9
                M=np.zeros((3,3))
                for x in group:M[int(x['row']),int(x['column'])]=float(x['RKF45_amplification'])
                record['amplification'].append({'h_factor':factor,'h_s':float(group[0]['h_s']),'spectral_radius':float(max(abs(np.linalg.eigvals(M))))})
            result['tail_stability'].append(record)
        # Rescale chemical coordinates to 1e-7 MeV before taking a matrix norm.
        scale=np.diag([1.,1e-7,1e-7]);inv=np.diag([1.,1e7,1e7])
        change=float(np.linalg.norm(inv@(matrices[0]-matrices[1])@scale,np.inf)/np.linalg.norm(inv@matrices[1]@scale,np.inf))
        for record in result['tail_stability']:
            if record['endpoint']==endpoint:record['jacobian_refinement_relative_change']=change
    result['units']='x dimensionless; eta in MeV; raw Jacobian entry J_ij has units (state_i/state_j)/s; derivative_per_s is only literally per second for equal-unit components. Eigenvalues are per second. Amplification matrix entries use the corresponding component ratio.'
    result['pass']=all(v['pass'] for v in result['run4'].values())
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=pathlib.Path);a=p.parse_args();r=assess(a.root)
    out=a.root/'assessment.json';temp=out.with_suffix('.json.partial');assert not out.exists() and not temp.exists()
    content=json.dumps(r,indent=2,allow_nan=False)+'\n'
    with temp.open('w') as f:f.write(content);f.flush();os.fsync(f.fileno())
    assert temp.read_text()==content and json.loads(temp.read_text())==r
    temp.rename(out);print(content)
