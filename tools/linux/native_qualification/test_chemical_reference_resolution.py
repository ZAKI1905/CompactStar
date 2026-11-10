#!/usr/bin/env python3
"""Analytic tests, never retained-data reference integrations."""
import unittest
from decimal import Decimal as D, localcontext

from resolve_chemical_reference import aligned_rk4, segment_value, sensitivity_bound, KS
from verify_chemical_reference import Chemistry, PRECISION, PI


class AlignedTests(unittest.TestCase):
    def setUp(self):
        self.context = localcontext()
        self.context.__enter__().prec = PRECISION

    def tearDown(self):
        self.context.__exit__(None, None, None)

    def test_tiny_segment_exact_cubic_forcing_integral(self):
        # A C1 piecewise cubic with an extremely steep final piece. Uniform
        # grids cannot resolve it. Integral on a segment is h*(x0+x1)/2
        # + h^2*(dx0-dx1)/12. RK4/Simpson integrates this forcing exactly.
        nodes = [(D(7), D(1), D(2)), (D(9), D(3), D(-1)),
                 (D(9)+D('1e-30'), D(5), D(0))]
        exact = sum((b[0]-a[0])*(a[1]+b[1])/2 + (b[0]-a[0])**2*(a[2]-b[2])/12
                    for a, b in zip(nodes, nodes[1:]))
        for k in KS:
            actual = aligned_rk4(lambda x, y: (x, -2*x), nodes, (D(11), D(13)), k)
            self.assertLess(abs(actual[0]-(11+exact)), D('1e-105'))
            self.assertLess(abs(actual[1]-(13-2*exact)), D('1e-105'))

    def test_normalized_coordinate_and_offset(self):
        nodes = [(D(1024), D(3), D(0)), (D(1024)+D('1e-20'), D(3), D(0))]
        expected = (D('3.5e-20'), D('-3.5e-20'))
        actual = aligned_rk4(lambda x, y: (x, -x), nodes, (D(0), D(0)), 8, D('.5'))
        self.assertEqual(actual, expected)
        self.assertEqual(segment_value(nodes[0], nodes[1], D(0)), nodes[0][1])
        self.assertEqual(segment_value(nodes[0], nodes[1], D(1)), nodes[1][1])

    def test_fourth_order_across_unequal_segments(self):
        nodes = [(D(0), D(0), D(0)), (D('.2'), D(0), D(0)), (D(1), D(0), D(0))]
        values = [aligned_rk4(lambda x, y: y, nodes, (D(1), D(2)), k) for k in KS]
        errors = [abs(v[0]-D(1).exp()) for v in values]
        for old, new in zip(errors, errors[1:]):
            self.assertTrue(D('.06') < new/old < D('.07'))
        last = abs(values[-1][0]-values[-2][0])
        self.assertLess(errors[-1], 2*last)

    def test_positive_polynomial_matches_reaction(self):
        model = Chemistry(dict(Z00=2, Z01=.1, Z10=.1, Z11=3,
                               Ltilde_Me=1e-60, Ltilde_Mmu=2e-60), dict(t_e=0, t_mu=0))
        T, eta = D('1.2e8'), D('.02')
        for j in range(2):
            polynomial = sum(model.luminosity[j]*c*eta**(2*k+1)*T**(6-2*k)
                /(model.kberg*11513*model.kb**(2*k+1)*PI**(2*k+2))
                for k, c in enumerate((14680,7560,840,24)))
            direct = model.luminosity[j]/model.kberg*T**7*model.reaction(T, eta)
            self.assertLess(abs(polynomial/direct-1), D('1e-105'))

    def test_sensitivity_zero_reaction_and_domain(self):
        model = Chemistry(dict(Z00=2, Z01=0, Z10=0, Z11=3,
                               Ltilde_Me=0, Ltilde_Mmu=0), dict(t_e=0, t_mu=0))
        bound = sensitivity_bound(model, [(D(0),D(0),D(0)), (D(1),D(0),D(0))], (D(1),D(2)))
        self.assertEqual(bound['status'], 'PASS')
        self.assertEqual(bound['uniform_thermal_bound_MeV'], 0)
        self.assertEqual(bound['lipschitz_per_s'], 0)

    def test_sensitivity_dominates_known_scalar_linear_case(self):
        # Tiny eta suppresses higher powers. Independently evaluate analytic
        # leading-term solution of y'=-a*T^6*y, with a thermal perturbation.
        model = Chemistry(dict(Z00=1e-10, Z01=0, Z10=0, Z11=1e-10,
                               Ltilde_Me=1e-70, Ltilde_Mmu=1e-70), dict(t_e=0, t_mu=0))
        eta, H = D('1e-30'), D(10)
        bound = sensitivity_bound(model, [(D(0),D(0),D(0)), (H,D(0),D(0))], (eta, eta))
        self.assertEqual(bound['status'], 'PASS')
        rate = model.z[0][0]*model.luminosity[0]*14680*D('1e48')/(model.kberg*11513*model.kb*PI**2)
        difference = eta*abs((-rate*H*D('6e-12').exp()).exp()-(-rate*H).exp())
        self.assertGreater(difference, 0)
        self.assertGreater(bound['uniform_thermal_bound_MeV'], difference)


if __name__ == '__main__':
    unittest.main(verbosity=2)
