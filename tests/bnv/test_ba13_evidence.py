"""Assessor fault tests; the mutation runner executes these against altered code."""
import copy
import os
import pathlib
import sys
import unittest
import numpy as np

SOURCE=pathlib.Path(__file__).resolve().parents[2]/'tools/linux/ba13_evidence'
sys.path.insert(0,os.environ.get('BA13_EVIDENCE_MODULES',str(SOURCE)))
from evidence import contribution,qualify,m1_decision,FRACTIONS
from reference import interpolation,paired_comparison


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.provenance=[dict(path='test-fixture',sha256='0'*64)]
        self.h=dict(drift_bound=.1,threshold=1.,utilization=.1)
        self.registry=dict(accepted_conditions={},accepted_standards={})
        self.terms=[contribution(c,c,0.,evidence_type='analytic_bound',
                   provenance=self.provenance,derivation='fixture exact zero') for c in FRACTIONS]

    def check(self,terms=None):
        return qualify(self.h,terms or self.terms,self.registry,frozen_historical=self.h)

    def test_no_default_evidence_class(self):
        with self.assertRaises(TypeError):
            contribution('bad','integration',1e-8,provenance=self.provenance)

    def test_nonzero_provenance_required(self):
        with self.assertRaises(ValueError):
            contribution('bad','integration',1e-8,evidence_type='measured_difference',provenance=[])

    def test_analytic_needs_derivation(self):
        with self.assertRaises(ValueError):
            contribution('bad','integration',1e-8,evidence_type='analytic_bound',provenance=self.provenance)

    def test_evidence_class_preserved(self):
        t=contribution('x','integration',1e-8,evidence_type='conditional_estimate',
                       provenance=self.provenance,condition_id='test-hypothesis')
        self.assertEqual(t['evidence_type'],'conditional_estimate')

    def test_conditional_acceptance_itemized(self):
        self.terms[3]=contribution('test','integration',1e-8,evidence_type='conditional_estimate',
                     provenance=self.provenance,condition_id='test-hypothesis',authority='fixture-owner')
        self.assertEqual(self.check()['status'],'UNQUALIFIED')
        self.registry['accepted_conditions']['test-hypothesis']='different-decision'
        self.assertEqual(self.check()['status'],'UNQUALIFIED')
        self.registry['accepted_conditions']['test-hypothesis']='fixture-owner'
        self.assertEqual(self.check()['status'],'PASS')

    def test_numerical_standard_not_automatic(self):
        self.terms[3]=contribution('test','integration',1e-8,evidence_type='validated_numerical_estimate',
                     provenance=self.provenance,standard_id='test-standard')
        self.assertEqual(self.check()['status'],'UNQUALIFIED')

    def test_missing_term_and_category_refuse(self):
        t=copy.deepcopy(self.terms);t[3]['value']=None;t[3]['missing']='no reference'
        self.assertEqual(self.check(t)['status'],'UNQUALIFIED')
        with self.assertRaises(ValueError):self.check(t[1:])

    def test_analytic_cannot_upgrade_inputs(self):
        self.terms.append(contribution('propagation','normalization',1e-8,evidence_type='analytic_bound',
                          provenance=self.provenance,derivation='valid conditional theorem',dependencies=['absent']))
        self.assertEqual(self.check()['status'],'UNQUALIFIED')

    def test_wrong_reference_rejected_by_comparison(self):
        p=np.array([[1.],[1.2]]);r=p.copy()
        with self.assertRaises(ValueError):paired_comparison(p,r,.01,reference_kind='R')
        with self.assertRaises(ValueError):interpolation([.5],[0,1],p,reference_kind='R')

    def test_achieved_coordinate_and_linear_remainder(self):
        y,_,_=interpolation([.3],[0,.2,1],np.array([[0.],[.04],[1.]]))
        self.assertAlmostEqual(float(y[0,0]),.16)
        self.assertLessEqual(abs(float(y[0,0])-.3**2),2*(.3-.2)*(1-.3)/2+1e-16)
        with self.assertRaises(ValueError):interpolation([1.01],[0,1],[[0],[1]])

    def test_actual_paired_reference_and_root_charge(self):
        p=np.array([[1.],[1.2],[1.4]]);x=np.array([[1.],[1.1],[1.2]])
        e,d,root=paired_comparison(p,x,.025)
        self.assertAlmostEqual(float(d[-1]),.2)
        self.assertEqual(root,.025)
        with self.assertRaises(ValueError):paired_comparison(p,x,None)

    def test_M1_budget_not_resolution(self):
        self.assertEqual(m1_decision(1e-14,1e-8,1e-5,complete=True,accepted=True),'PASS')
        self.assertEqual(m1_decision(1e-4,None,1e-5,complete=False,accepted=False),'FAIL')
        self.assertEqual(m1_decision(1e-14,1e-8,1e-5,complete=False,accepted=True),'UNQUALIFIED')

    def test_frozen_history_and_invalid_values(self):
        with self.assertRaises(ValueError):
            qualify(dict(self.h,threshold=2),self.terms,self.registry,frozen_historical=self.h)
        for bad in (-1,float('nan'),float('inf')):
            with self.assertRaises(ValueError):contribution('x','integration',bad,evidence_type='measured_difference',provenance=self.provenance)


if __name__=='__main__':unittest.main(verbosity=2)
