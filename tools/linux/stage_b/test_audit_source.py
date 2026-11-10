#!/usr/bin/env python3
"""Verify isolated instrumentation retains every production Stepper operation."""
import pathlib,tempfile,unittest
from audit_trace import generate
ROOT=pathlib.Path(__file__).resolve().parents[3]
class AuditSourceTests(unittest.TestCase):
 def test_only_class_name_and_passive_boundary_hooks(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=pathlib.Path(tmp);generate(ROOT,p)
   instrumented=(p/'AuditStepper.hpp').read_text()
   self.assertEqual(instrumented.count('::AuditApplyBegin(t,y);'),1)
   self.assertEqual(instrumented.count('::AuditApplyEnd(t,y,evolve_->failed_steps-rejected);'),1)
   restored=instrumented.replace('AuditStepper','BreakpointStepper').replace('::AuditApplyBegin(t,y);\n        ','').replace('\n        ::AuditApplyEnd(t,y,evolve_->failed_steps-rejected);','')
   self.assertEqual(restored,(ROOT/'CompactStar/Physics/BNV/BreakpointStepper.hpp').read_text())
if __name__=='__main__':unittest.main()
