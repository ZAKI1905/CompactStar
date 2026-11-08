#!/usr/bin/env python3
"""Exact binary64 and 80-digit witness; no stars, rates or equations are changed."""
import argparse
import csv
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path

ROOT=Path('/mnt/sdd/zaki/CompactStar/linux-transition-20261010')
REPO=Path(__file__).resolve().parents[3]
HISTORICAL=ROOT.parent/'incoming/CompactStar-EKU-handoff-20261010/scientific/phase6/coefficients.tsv'
NATIVE=ROOT/'ba13-recovery/run1/raw/coefficients.tsv'

def witness():
    pins=json.loads((REPO/'docs/validation/results-first-bnv-pilot/evidence/preflight.json').read_text())['input_sha256']
    if hashlib.sha256(NATIVE.read_bytes()).hexdigest()!=pins[str(NATIVE)]:
        raise ValueError('native coefficient identity changed')
    result={}
    with localcontext() as ctx:
        ctx.prec=80
        for label,p,literal in [('historical',HISTORICAL,-2.4136520263641375e37),('native',NATIVE,-2.413652026364133e37)]:
            with p.open() as f:row=next(csv.DictReader(f,delimiter='\t'))
            if int(row['index'])!=0:raise ValueError('not the B0 row')
            B=float(row['B_solved']);fraction=-1e-12;year=31557600.0
            derived=fraction*B/year
            if derived!=literal:raise ValueError(label+' derived binary64 identity changed')
            exact=Decimal.from_float(fraction)*Decimal.from_float(B)/Decimal.from_float(year)
            physical=Decimal('-1e-12')*Decimal.from_float(B)/Decimal('31557600')
            rounding=(Decimal.from_float(derived)-exact)/Decimal.from_float(math.ulp(derived))
            if abs(rounding)>1:raise ValueError('unexpected two-operation rounding')
            result[label]=dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                B0_count=B,B0_hex=B.hex(),fractional_rate_hex=fraction.hex(),seconds_per_year=year,
                derived_Bdot_count_s=derived,derived_hex=derived.hex(),expected_literal_hex=literal.hex(),
                exact_binary_inputs_witness=str(exact),decimal_fraction_witness=str(physical),
                binary64_rounding_error_ulp=str(rounding))
        Bm=Decimal.from_float(result['historical']['B0_count']);Bn=Decimal.from_float(result['native']['B0_count'])
        dm=Decimal.from_float(result['historical']['derived_Bdot_count_s']);dn=Decimal.from_float(result['native']['derived_Bdot_count_s'])
        result['comparison']=dict(B0_relative_difference=str((Bn-Bm)/Bm),
            B0_difference_ulp=str((Bn-Bm)/Decimal.from_float(math.ulp(float(Bm)))),
            Bdot_relative_difference=str((dn-dm)/dm),
            Bdot_difference_ulp=str((dn-dm)/Decimal.from_float(math.ulp(float(dm)))),
            duration_s=462269531250,expected_fractional_depletion=str(Decimal('1e-12')*462269531250/31557600),
            conclusion='Both exact literals are reproduced by the unchanged fractional-drive formula with their own authenticated B0. No new microscopic rate or tolerance.',
            historical_assertion_modified=False,trajectories_run=0)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    with a.output.open('x') as f:json.dump(witness(),f,indent=2);f.write('\n')
