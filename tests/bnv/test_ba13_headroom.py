"""Faults exercise quantitative counterexamples, not only required-field flags."""
import pathlib,sys,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]/'tools/linux/ba13_headroom'))
from enclosure import *

class EnclosureTests(unittest.TestCase):
 def test_budget_units_and_completeness(self):
  for h in [2.106701e-7,9.8540992e-5]:self.assertAlmostEqual(sum(budget(h).values())/h,1.)
 def test_support_edge_omission(self):
  # Constant integrand on a moving support has exact integral equal to its edge.
  edge0,edgei=1.,1.00003;nom0,nomi=1.00009,1.00011
  truth=abs(edgei/edge0-1);observed=abs(nomi-nom0)
  self.assertGreater(truth,observed)
  correction=paired_affine([edgei-nomi],[edge0-nom0])
  self.assertGreaterEqual(normalized_drift(observed,correction,nom0,abs(edge0-nom0)),truth)
 def test_incorrect_shared_cancellation(self):
  self.assertEqual(paired_affine([2e-5],[1e-5]),1e-5)
  self.assertGreater(abs((1.00002/1.00001)-1),0.)
  # Multiplicative thermal factor really does cancel in a same-knot ratio.
  self.assertAlmostEqual((1.000655*1.000003)/(1.000655*1.)-1,3e-6)
 def test_missing_B_functional(self):
  # Computed B=2H+0.01; target=2. Root residual zero, actual root H=1.
  self.assertEqual(root_error(0,.01,2),.005)
  self.assertGreater(abs(.995-1),root_error(0,0,2))
 def test_inverse_remainder(self):
  g=np.diag([1.,2.,3.]);e=np.diag([.1,.2,.3]);a,total,rem,rho=inverse_majorant(g,e)
  actual=abs(np.linalg.inv(g-e)-a)
  self.assertTrue(np.all(actual<=total*(1+1e-14)))
  self.assertGreater(actual[0,0],(total-rem)[0,0])
 def test_changing_inverse_shared_error(self):
  g0=np.diag([1.,2.,3.]);gi=np.diag([1.2,2.1,2.8]);e=np.diag([.1,.1,.1])
  base=L.T@(np.linalg.inv(gi)-np.linalg.inv(g0))@L
  pert=L.T@(np.linalg.inv(gi-e)-np.linalg.inv(g0-e))@L
  b=paired_inverse(g0,gi,e)
  self.assertGreater(np.max(abs(pert-base)),0.)
  self.assertTrue(np.all(abs(pert-base)<=b['total']))
 def test_underestimated_center_and_integration(self):
  # y'=y, y(0)=1: missing center value and constant defect have known transport.
  center,defect=.001,.002
  exact=center*math.e+defect*(math.e-1)
  self.assertGreater(exact,center+defect)
  self.assertAlmostEqual(exact,math.e*(center+defect*(1-math.exp(-1))))
 def test_covariance_not_physical(self):
  terms={k:dict(value=0.,evidence_type='analytic_bound') for k in FRACTIONS}
  terms['integration']['source']='fit_covariance'
  with self.assertRaises(ValueError):qualify(dict(drift_bound=0.,threshold=1e-4),terms)
 def test_all_missing_refuses(self):
  terms={k:dict(value=None,evidence_type='unresolved') for k in FRACTIONS}
  self.assertEqual(qualify(dict(drift_bound=0.,threshold=1e-4),terms)['status'],'UNQUALIFIED')
 def test_frozen_threshold_failure(self):
  terms={k:dict(value=1e-5,evidence_type='analytic_bound') for k in FRACTIONS}
  h=dict(drift_bound=0.,threshold=1e-4)
  self.assertEqual(qualify(h,terms)['status'],'FAIL');self.assertEqual(h['threshold'],1e-4)
 def test_unseen_bump_and_curvature(self):
  # Bump on [0,h] vanishes at both sampled ends. Sampled maxima cannot bound it.
  h=5e-8;amp=2e-4
  K=2*amp*math.pi**2/h**2
  self.assertGreater(amp,0.)
  self.assertGreaterEqual(between_samples(K,h),amp)
 def test_source_residual_correlation(self):
  x=np.linspace(-5e-7,-1e-6,11);y=1+2*x
  h=projection(x,np.full(11,1e-8),y)
  self.assertLess(np.max(residual_error(h,np.ones((11,1))*1e-5,np.zeros(11))),1e-19)
  self.assertLess(np.max(abs(h@y-y)),1e-14)
 def test_normalization_and_singular_refusal(self):
  self.assertEqual(normalized_drift(.1,.01,1,.1),.11/.9)
  with self.assertRaises(ValueError):normalized_drift(0,0,1,1)
  with self.assertRaises(ValueError):inverse_majorant(np.eye(3),np.eye(3))

if __name__=='__main__':unittest.main(verbosity=2)
