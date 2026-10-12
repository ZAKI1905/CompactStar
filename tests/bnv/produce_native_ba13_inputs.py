#!/usr/bin/env python3
"""Linux BA13 adapter: recompute the Phase-5C tail witness at each actual cut.

Original producer bytes remain in produce_frozen_certificates.py. The inherited
B0 background/center characterization is retained and explicitly identified; it
is not a new per-star refinement proof. No Phase-5 formula/assertion is changed.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import numpy as np
from produce_frozen_certificates import thermal_source


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,required=True)
    p.add_argument('--profiles',type=Path,required=True)
    p.add_argument('--phase5-characterization',type=Path,required=True)
    a=p.parse_args();source=a.source_root.resolve()
    sys.path.insert(0,str(source/'tests/analysis'))
    from chemical_production_evidence import certificates, transport
    from chemical_trackr_budget import Model, tail
    fixed=json.loads((source/'docs/validation/phase5c2_preproduction_evidence.json').read_text())
    parent=json.loads(a.phase5_characterization.read_text())
    dirs=sorted(a.profiles.glob('target-*'),key=lambda x:int(x.name.split('-')[-1]))
    if [x.name for x in dirs] != [f'target-{i}' for i in range(21)]:
        raise ValueError('exact 21-star input set required')
    for i,directory in enumerate(dirs):
        d=np.loadtxt(directory/'profile.tsv',skiprows=1);model=Model(directory)
        characterized=copy.deepcopy(parent)
        # g is used only to report the tail's propagated Z diagnostic. The tail
        # witness and analytic enclosure are determined from this actual profile.
        measured,_,_=tail(d,model,np.array(parent['G']))
        characterized['tail']=measured
        certificate,_=certificates(directory,characterized,fixed['goals'])
        transport(certificate,directory/'certificate.txt',fixed)
        receipt=dict(schema='compactstar.ba13.current-cut-input.v1',
            profile_sha256=hashlib.sha256((directory/'profile.tsv').read_bytes()).hexdigest(),
            inherited_characterization_sha256=hashlib.sha256(a.phase5_characterization.read_bytes()).hexdigest(),
            inherited_components=['background_characterization','provider_nodal_comparison','equilibrium_anchor','centre'],
            recomputed_component='tail: unchanged chemical_trackr_budget.tail on actual target profile',
            limitation='Inherited B0 errors are not independent per-target structural refinement evidence',
            measured_tail=measured)
        (directory/'tail-witness.json').write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+'\n')
        (directory/'certificate.json').write_text(json.dumps(certificate,indent=2,sort_keys=True,allow_nan=False)+'\n')
        thermal_source(directory)
        print(f'NATIVE_FROZEN_CHEMICAL_CERTIFICATE {i} PASS_CURRENT_CUT',flush=True)

if __name__=='__main__':main()
