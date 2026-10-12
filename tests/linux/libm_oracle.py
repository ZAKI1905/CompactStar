#!/usr/bin/env python3
"""100-digit exact-input scalar accuracy, near-tie and held-out cases; no Mac emulation."""
import argparse
import ctypes
from decimal import Decimal as D, localcontext
import json
import math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--libm',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
lib=ctypes.CDLL(a.libm)
lib.pow.restype=ctypes.c_double;lib.pow.argtypes=[ctypes.c_double,ctypes.c_double]
lib.log10.restype=ctypes.c_double;lib.log10.argtypes=[ctypes.c_double]
rows=[]
with localcontext() as ctx:
 ctx.prec=120
 L=float.fromhex('0x1.22fad6cb53501p-10')
 for x in [L,1-(1-L),-2.875,-0.3125,0.8125,3.0625,7.375,14.125,16.75]:
  exact=(D(10).ln()*D.from_float(x)).exp(); actual=lib.pow(10.,x)
  error=abs(D.from_float(actual)-exact);bound=D.from_float(math.ulp(float(exact)))
  assert error<=bound
  rows.append(dict(function='pow10',input=x.hex(),actual=actual.hex(),correctly_rounded=float(exact).hex(),exact=str(exact),absolute_error=str(error),bound=str(bound)))
 for x in [0.03125,0.8125,1.0625,2.375,13.25,139.75,1e14,3e15]:
  exact=D.from_float(x).ln()/D(10).ln();actual=lib.log10(x)
  error=abs(D.from_float(actual)-exact);bound=D.from_float(math.ulp(float(exact)))
  assert error<=bound
  rows.append(dict(function='log10',input=x.hex(),actual=actual.hex(),correctly_rounded=float(exact).hex(),exact=str(exact),absolute_error=str(error),bound=str(bound)))
a.output.write_text(json.dumps(dict(status='PASS_CANDIDATE_SCALAR_ONLY',oracle='Decimal exp/ln at 120 digits; exact binary inputs',rows=rows),indent=2)+'\n')
print('PASS near-tie and 15 additional independent scalar inputs; contour topology remains separately assessed')
