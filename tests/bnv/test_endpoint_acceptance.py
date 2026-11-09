#!/usr/bin/env python3
"""Scope and contradiction refusal tests for the owner's endpoint condition."""
import copy
import importlib.util
from pathlib import Path
import unittest
from decimal import localcontext
ROOT=Path(__file__).resolve().parents[2]
s=importlib.util.spec_from_file_location('endpoint',ROOT/'tools/linux/results_first/endpoint_acceptance.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class Authority(unittest.TestCase):
    def setUp(self):
        self.registry=m.read(m.DOC/'acceptance-registry.json')
        self.bounds=m.read(m.WITNESS/'bounds.json')
    def test_accepted_scope(self):m.authenticate(self.registry,self.bounds)
    def test_broadened_scope_refused(self):
        for k,v in [('scope','another-star'),('maximum_extension_x','3.023e-26'),('physical_interval_x',['-2e-6','0']),('rigorous_endpoint_bound',True),('paired_trajectories_authorized',True)]:
            r=copy.deepcopy(self.registry);r[k]=v
            with self.subTest(k=k),self.assertRaises(ValueError):m.authenticate(r,self.bounds)
    def test_over_cap_witness_refused(self):
        self.bounds['eta']='3.0220001e-26'
        with self.assertRaises(ValueError):m.authenticate(self.registry,self.bounds)
    def test_known_transition_refused(self):
        data=m.read(m.BASE/'cut-assessment/retained-inputs.json.gz')
        data['centers']['0']['2048']['Hs']='0.2'
        with localcontext() as ctx:
            ctx.prec=110
            with self.assertRaisesRegex(ValueError,'topology contradiction'):
                m.transition_check(self.bounds,data,m.read(ROOT/'docs/validation/ba13-evidence/evidence/continuation.json'))
    def test_nonpositive_baryon_derivative_refused(self):
        data=m.read(m.BASE/'cut-assessment/retained-inputs.json.gz')
        pair=next(iter(data['widths']['20'].values()))
        pair['1']['primary'][0]=pair['-1']['primary'][0]
        with localcontext() as ctx:
            ctx.prec=110
            with self.assertRaisesRegex(ValueError,'baryon derivative'):
                m.transition_check(self.bounds,data,m.read(ROOT/'docs/validation/ba13-evidence/evidence/continuation.json'))
if __name__=='__main__':unittest.main()
