#!/usr/bin/env python3
"""Known-answer crossing/propagation tests, independent of stellar trajectories."""
import unittest
from decimal import Decimal as D,localcontext
from piecewise import step,Coupled,MODEL,domain,SCALE

class PiecewiseTests(unittest.TestCase):
 def test_transverse_crossing_enclosure(self):
  with localcontext() as c:
   c.prec=110
   # x'=1+b*x on either side, continuous at 0 but derivative jumps 2 -> 5.
   def rhs(y):return (1+(D(2) if y[0]<0 else D(5))*y[0],D(0),D(0))
   for h in (D('.01'),D('.001'),D('.0001')):
    crossing=h/3;x0=((D(-2)*crossing).exp()-1)/2
    exact=((5*(h-crossing)).exp()-1)/5
    result=step(rhs,(x0,D(0),D(0)),h)[0]
    # Entire segment in [-.1,.1]: M=1.5, L=5.
    bound=D(5)*D('1.5')*h*h*(10*h).exp()
    self.assertLess(abs(result-exact),bound)
 def test_prefix_propagation_across_knot(self):
  with localcontext() as c:
   c.prec=110;h=D('.001');L=D(5);M=D('1.5');y=(D('-.02'),D(0),D(0));u=D('1e-12')
   initial=y[0]+u;crossing=-(1+2*initial).ln()/2
   def rhs(y):return (1+(D(2) if y[0]<0 else D(5))*y[0],D(0),D(0))
   for i in range(50):
    y=step(rhs,y,h);u=(L*h).exp()*u+L*M*h*h*(2*L*h).exp()
   exact=((5*(50*h-crossing)).exp()-1)/5
   self.assertLess(abs(y[0]-exact),u)
 def test_actual_cell_crossing_local_remainders(self):
  with localcontext() as c:
   c.prec=110
   for role in ('source','control'):
    m=Coupled(MODEL,role);box=domain(m);h=D('.01')
    for knot in (93,94,95,96):
     boundary=m.knots[knot];eta=(D('-2e-7'),D('-5e-7')) if role=='source' else (D(0),D(0))
     v=m.coupled_rhs((boundary,*eta))[0]
     initial=(boundary-v*h/3,*eta)
     coarse=step(m.coupled_rhs,initial,h);fine=initial
     for _ in range(32):fine=step(m.coupled_rhs,fine,h/32)
     bound=box['L']*box['M']*h*h*(2*box['L']*h).exp()
     self.assertLess(max(abs(a-b)/s for a,b,s in zip(coarse,fine,SCALE)),bound)
     self.assertLess((initial[0]-boundary)*(coarse[0]-boundary),0)

 def test_rk4_order_known_exponential(self):
  with localcontext() as c:
   c.prec=110;errors=[]
   for n in (8,16,32):
    y=(D(1),D(0),D(0));h=D('.1')/n
    for _ in range(n):y=step(lambda z:(z[0],D(0),D(0)),y,h)
    errors.append(abs(y[0]-D('.1').exp()))
   self.assertTrue(all(D(15)<a/b<D(17) for a,b in zip(errors,errors[1:])))

if __name__=='__main__':unittest.main()
