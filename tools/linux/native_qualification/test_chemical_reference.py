#!/usr/bin/env python3
"""Analytic tests for the independent D1 reference; no retained data required."""
import unittest
from decimal import Decimal as D, localcontext

from verify_chemical_reference import (
    Chemistry, Hermite, NS, PRECISION, convergence, corner_uncertainty, energy,
    fidelity, lift, rk4, ulp,
)


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.context = localcontext()
        self.context.__enter__().prec = PRECISION

    def tearDown(self):
        self.context.__exit__(None, None, None)

    def test_binary64_lift_and_ulp(self):
        self.assertEqual(lift('0.1'), D(3602879701896397) / 36028797018963968)
        self.assertNotEqual(lift('0.1'), D('0.1'))
        self.assertEqual(ulp(1), D(2)**-52)
        self.assertEqual(ulp(0), D.from_float(float.fromhex('0x0.0000000000001p-1022')))

    def test_constant_rhs_and_nonzero_time_origin(self):
        self.assertEqual(rk4(lambda t, y: (D(2), D(-3)), D(7), D(11), (D(1), D(2)), 64),
                         (D(9), D(-10)))

    def test_cubic_forcing_exact_integral(self):
        y = rk4(lambda t, y: (4*t**3, -8*t**3), D(0), D(2), (D(3), D(5)), 64)
        self.assertEqual(y, (D(19), D(-27)))

    def test_exponential_solution_and_fourth_order(self):
        values = {n: rk4(lambda t, y: y, D(0), D(1), (D(1), D(2)), n) for n in NS}
        exact = D(1).exp()
        errors = [abs(values[n][0]-exact) for n in NS]
        self.assertLess(errors[-1], D('4e-13'))
        for old, new in zip(errors, errors[1:]):
            self.assertGreater(new/old, D('0.06'))
            self.assertLess(new/old, D('0.065'))
        self.assertTrue(all(r['status']=='PASS' for r in convergence(values)))

    def test_energy_with_cross_term(self):
        # Z^-1 = [[3,-1],[-1,2]]/5; E = eta^T Z^-1 eta / 2.
        z = [[D(2), D(1)], [D(1), D(3)]]
        self.assertEqual(energy(z, (D(1), D(2))), D('0.7'))
        self.assertEqual(energy(z, (D(-1), D(-2))), D('0.7'))
        self.assertEqual(energy(z, (D(0), D(0))), D(0))

    def test_four_corner_energy(self):
        z = [[D(2), D(0)], [D(0), D(4)]]
        uncertainty, corners = corner_uncertainty(lambda y: energy(z, y), (D(2), D(4)), (D('0.1'), D('0.2')))
        self.assertEqual(len(corners), 4)
        self.assertEqual(uncertainty, D('0.3075'))

    def test_hermite_reproduces_cubic(self):
        def f(t):
            return 2*t**3-3*t**2+5*t+7
        def df(t):
            return 6*t**2-6*t+5
        path = Hermite([(t, f(t), df(t)) for t in map(D, (0, 1, 3))])
        for t in (D(0), D('0.25'), D(1), D('1.5'), D(3)):
            self.assertEqual(path(t), f(t))

    def test_zero_convergence_differences(self):
        same = {n: (D(0), D(2)) for n in NS}
        for record in convergence(same):
            self.assertEqual(record['status'], 'PASS')
            self.assertTrue(all(c['zero_case']=='both_zero' and c['ratio'] is None for c in record['checks']))
        changed = {64: (D(1), D(1)), 128: (D(1), D(1)), 256: (D(2), D(2)), 512: (D(3), D(3))}
        self.assertTrue(all(r['status']=='FAIL' for r in convergence(changed)))
        self.assertEqual(convergence(changed)[0]['checks'][0]['zero_case'], 'zero_denominator_nonzero_numerator')

    def test_convergence_floor(self):
        values = {64: (D(1), D(1)), 128: (D(1), D(1)),
                  256: (D(1)+D('1e-65'), D(1)+D('1e-65')),
                  512: (D(1)+D('4e-65'), D(1)+D('4e-65'))}
        self.assertTrue(all(r['status']=='PASS' for r in convergence(values)))

    def test_reaction_zero_and_odd_symmetry(self):
        model = Chemistry(dict(Z00=2, Z01=0, Z10=0, Z11=4, Ltilde_Me=1, Ltilde_Mmu=1), dict(t_e=0, t_mu=0))
        self.assertEqual(model.rhs(D(0), (D(0), D(0))), (D(0), D(0)))
        eta = D('0.01')
        self.assertEqual(model.reaction(D('1e8'), -eta), -model.reaction(D('1e8'), eta))

    def test_fidelity_exact_two_ulp_boundary(self):
        class Stub:
            def rhs(self, x, eta):
                return D(1)+2*ulp(1), D(1)+3*ulp(1)
        result = fidelity(Stub(), [dict(context='7', t_s='0', x='0', e='0', mu='0', de='1', dmu='1')])
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(len(result['failures']), 1)
        self.assertEqual(result['failures'][0]['line'], 2)
        self.assertEqual(result['failures'][0]['component'], 'dmu')
        self.assertEqual(result['failures'][0]['difference_ulp'], 3)


if __name__ == '__main__':
    unittest.main(verbosity=2)
