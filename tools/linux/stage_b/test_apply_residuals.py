#!/usr/bin/env python3
"""Known-answer arithmetic propagation through discarded trials and RK retries."""
import csv,pathlib,tempfile,unittest
from decimal import Decimal as D,localcontext
from unittest.mock import patch
import apply_residuals as a
class Linear:
 def __init__(self,*args):pass
 def coupled_rhs(self,y):return (y[0],D(0),D(0))
def box(model):return [[D(1),D(0),D(0)],[D(0)]*3,[D(0)]*3],dict(x_min=D(-10),x_max=D(10),eta_abs=[D(10)]*2)
class ArithmeticTraceTests(unittest.TestCase):
 def check(self,discard=False,reject=False,cached=False):
  with localcontext() as c,tempfile.TemporaryDirectory() as tmp:
   c.prec=110;p=pathlib.Path(tmp);A,W=a.old.tableau('rkf45');n=len(W);trace=[];calls=[]
   def row(t,x,cell=0):return dict(context=1,t_s=float(t),cell=cell,x=float(x),e=0,mu=0,dx=float(x),de=0,dmu=0)
   trace.append(row(0,1,-1))
   def attempt(h):
    fs=[];ys=[]
    for j in range(n):
     y=D(1)+h*sum(A[j][k]*fs[k] for k in range(j));ys.append(y);fs.append(y)
    endpoint=D(1)+h*sum(W[j]*fs[j] for j in range(n))
    return [row(h*sum(A[j]),ys[j]) for j in range(1,n)]+[row(h,endpoint)],endpoint
   def call(h,rejected=0):
    begin=len(trace)
    if not cached:trace.append(row(0,1))
    if rejected:
     rejected_rows,_=attempt(D(1));trace.extend(rejected_rows)
    rr,y=attempt(h);trace.extend(rr)
    calls.append(dict(context=1,apply=len(calls)+1,left=0,right=float(h),row_begin=begin,row_end=len(trace),rejected=rejected,x0=1,e0=0,mu0=0,x1=float(y),e1=0,mu1=0))
    return y
   if cached:
    trace.append(row(0,1,0))
    probe=row(0,1,-1);probe['dx']=100;trace.append(probe)
   if discard:call(D(1))
   exact=call(D('.1'),int(reject));native=float(exact)
   mapping=[dict(context=1,left=0,right=.1,h=.1,x=native,e=0,mu=0,replay_x=native,replay_e=0,replay_mu=0,steps=1,rejected=int(reject))]
   for name,data in [('local-rhs.tsv',trace),('applies.tsv',calls),('maps.tsv',mapping)]:
    with (p/name).open('w') as f:w=csv.DictWriter(f,fieldnames=data[0],delimiter='\t');w.writeheader();w.writerows(data)
   with patch.object(a,'Coupled',Linear),patch.object(a,'jacobian',box):r=a.audit('source','manufactured',p,[D('1e-10'),D(0),D(0)],'rkf45')
   # Linear RK stability polynomial multiplies an incoming perturbation.
   true_discrete_error=abs(D.from_float(native)-exact*(1+D('1e-10')))
   self.assertLessEqual(true_discrete_error,r['final_enclosure'][0])
   self.assertLess(r['final_enclosure'][0],D('1.5e-10'))
   self.assertEqual(r['rejected_attempts'],int(reject))
   self.assertEqual(r['event_trial_resets'],int(discard))
 def test_simple_map(self):self.check()
 def test_cached_derivative_ignores_selection_probe(self):self.check(cached=True)
 def test_discarded_event_trial_not_serial(self):self.check(discard=True)
 def test_rejected_attempt_not_serial(self):self.check(reject=True)
 def test_combined_discard_and_rejection(self):self.check(discard=True,reject=True)
if __name__=='__main__':unittest.main()
