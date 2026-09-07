#!/usr/bin/env python3
"""Durable GC9 analytic oracle + explicit wrong routes; future API adapter goes here."""
import json
import subprocess
import sys
sys.dont_write_bytecode = True
import numpy as np
from chemical_curved_reference import references


def radial_integral(order, power=-1, volume=-0.5, variable=False):
    x,w = np.polynomial.legendre.leggauss(order)
    r = 6*(x+1)
    mass = 1.8*(r/12)**3
    lapse = (3*np.sqrt(0.7)-np.sqrt(1-2*mass/r))/2
    return np.dot(6*w,4*np.pi*r*r*lapse**power*(1-2*mass/r)**volume*
                  (1+(r/12)**2 if variable else 1))


def run():
    output = subprocess.check_output([sys.argv[1], 'gc9'], text=True)
    production = {}
    for line in output.splitlines():
        if line.startswith('GC9\t'):
            _, shape, order, mutant, value = line.split('\t')
            production[int(shape), int(order), int(mutant)] = float(value)
    assert len(production) == 42
    ref70,ref100 = references(70),references(100)
    # Predeclared arithmetic envelope: 4096 elementary operations at binary64 u.
    # This fixture is smooth and well conditioned; it is not a general G tolerance.
    u = np.finfo(float).eps/2
    gamma = 4096*u/(1-4096*u)
    report = {'classification':'INDEPENDENT ANALYTIC ORACLE',
              'reference_method':'closed primitive + 70/100-digit theta integration',
              'production_path_exercised':True, 'reference':ref100,
              'roundoff_relative_envelope':gamma, 'shapes':[]}
    for shape in (False,True):
        expected = float(ref100[int(shape)])
        assert float(ref70[int(shape)]) == expected
        values = [production[int(shape), n, 0] for n in (16,32,64)]
        error = abs(values[-1]-expected)/expected
        truncation = abs(values[-1]-values[-2])/expected
        assert error <= gamma+truncation
        mutants = {'M6_omit_lapse':(0,-.5),'M7a_extra_inverse_lapse':(-2,-.5),
                   'M7b_nu_as_2Phi':(-.5,-.5),'M8_omit_volume':(-1,0),
                   'M19_sign_flipped_lapse':(1,-.5),'M20_inverse_volume':(-1,.5)}
        separations = {name:(production[int(shape), 64, index]-expected)/expected
                       for index,name in enumerate(mutants, start=1)}
        # Detector margin compared only to a predeclared numerical error, never
        # set a tolerance from any observed mutant percentage.
        assert all(abs(v)>2*(gamma+truncation) for v in separations.values())
        report['shapes'].append({'variable':shape,'correct_relative_error':error,
            'order32_64_relative_difference':truncation,'mutant_signed_relative_separation':separations})
    return report


if __name__ == '__main__':
    print(json.dumps(run(),indent=2))
