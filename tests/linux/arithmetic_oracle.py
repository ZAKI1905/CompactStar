#!/usr/bin/env python3
"""Independent arbitrary precision oracle, exact binary inputs; stdlib-only on Mac.

Machin pi = 16 atan(1/5)-4 atan(1/239); alternating-series remainder <
next term. Iteration stops below 10^(-precision-8), with 16 guard digits.
Onset uses the two-body Kallen momentum formula, independently of production mu.
80 and 120 decimal digit runs must give the same rounded binary64 result.
"""
import argparse
from decimal import Decimal as D, localcontext
import json
import math
from pathlib import Path
import subprocess


def pi(precision):
    def atan_inverse(n):
        x = D(1) / n
        term, total, i = x, x, 1
        while True:
            term *= -x*x
            add = term / (2*i+1)
            total += add
            if abs(add) < D(10) ** (-precision-8):
                return total
            i += 1
    return 16*atan_inverse(5)-4*atan_inverse(239)


def calculate(constants, precision):
    with localcontext() as ctx:
        ctx.prec = precision + 16
        mn, mp, me, k = map(D.from_float, constants)
        momentum2 = (mn-mp-me)*(mn+mp+me)*(mn-mp+me)*(mn+mp-me)/(4*mn*mn)
        truth = momentum2*momentum2.sqrt()*k**3/(3*pi(precision)**2)
        return +truth


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--executable', type=Path)
    p.add_argument('--input', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    raw = a.input.read_text() if a.input else subprocess.check_output([str(a.executable.resolve())], text=True)
    rows = [r.split('\t') for r in raw.splitlines() if r.startswith(('constants\t','onset\t','inverse\t','muon_onset\t','types\t','local\t'))]
    constants = [float.fromhex(x) for x in next(r[1:] for r in rows if r[0]=='constants')]
    t80, truth = calculate(constants,80), calculate(constants,120)
    assert float(t80)==float(truth) and abs(t80-truth)<D('1e-85')
    onset = float.fromhex(next(r[1] for r in rows if r[0]=='onset'))
    checks=[]
    with localcontext() as ctx:
        ctx.prec=136
        rel=abs(D.from_float(onset)/truth-1)
        assert rel<D('3e-16'), ('PE-V6',str(rel))
        for row in rows:
            if row[0]!='inverse':
                continue
            mu, actual = map(float.fromhex,row[1:])
            x=D.from_float(mu); me=D.from_float(constants[2]); k=D.from_float(constants[3])
            p2=x*x-me*me
            expected=p2*p2.sqrt()*k**3/(3*pi(120)**2)
            error=abs(D.from_float(actual)-expected)
            bound=16*D(2)**-52*abs(expected)
            assert error<=bound, (row,str(error),str(bound))
            checks.append(dict(input=mu.hex(),actual=actual.hex(),expected=str(expected),absolute_error=str(error),bound=str(bound)))
    local_checks=[]
    with localcontext() as ctx:
        ctx.prec=136
        # Evaluate uses the binary64 reciprocal stored by the production object.
        hc=D.from_float(1/constants[3]); high_pi=pi(120); eps=D(2)**-52
        for row in rows:
            if row[0]!='local': continue
            n,m,k_a,mu_a,en_a,der_a=map(float.fromhex,row[1:])
            n,m=map(D.from_float,(n,m))
            k=hc*((3*high_pi**2*n).ln()/3).exp(); mu=(m*m+k*k).sqrt();x=k/m
            aa=x*(1+x*x).sqrt()*(1+2*x*x);bb=(x+(1+x*x).sqrt()).ln()
            pref=m**4/(8*high_pi**2*hc**3);energy=pref*(aa-bb)
            derivative=k*k/(3*mu*n)
            eb=32*eps*pref*(abs(aa)+abs(bb)) if x>=D('.01') else 32*eps*abs(energy)+pref*7*x**13/416
            for name,actual,expected,bound in [('momentum',k_a,k,32*eps*abs(k)),('mu',mu_a,mu,32*eps*abs(mu)),('energy',en_a,energy,eb),('derivative',der_a,derivative,32*eps*abs(derivative))]:
                error=abs(D.from_float(actual)-expected)
                assert error<=bound, (name,row,str(error),str(bound))
                local_checks.append(dict(quantity=name,n=str(n),mass=str(m),actual=actual.hex(),absolute_error=str(error),bound=str(bound)))
    result=dict(status='PASS_CANDIDATE_ONLY',oracle='Decimal Machin pi and Kallen momentum, 80/120 digits',
        constants=[v.hex() for v in constants],exact_onset=str(truth),correctly_rounded=float(truth).hex(),
        actual=onset.hex(),relative_error=str(rel),onset_limit='3e-16 (unchanged)',inverse_checks=checks,
        local_checks=local_checks,types=next(r[1:] for r in rows if r[0]=='types'),raw=raw)
    text=json.dumps(result,indent=2)+'\n'
    if a.output: a.output.write_text(text)
    print(text)

if __name__=='__main__': main()
