#!/usr/bin/env python3
"""Recovery integrity and falsifiers, independent of any trajectory."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/linux'))
from ba13_stage import rows,verify


class RecoveryTests(unittest.TestCase):
    def test_recovered_bytes(self):
        records=json.loads((ROOT/'docs/validation/ba13-recovery/recovered-identities.json').read_text())
        for item in records:
            self.assertEqual(hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest(),item['sha256'])
            archived=ROOT/'docs/validation/ba13-recovery/historical'/Path(item['path']).name
            self.assertEqual(archived.read_bytes(),(ROOT/item['path']).read_bytes())

    def test_truncated_and_nonfinite_tables_refuse(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'sample.tsv'
            p.write_text('index\tvalue\n0\t1\n')
            with self.assertRaises(RuntimeError):rows(p,21)
            for value in ('nan','inf','-inf'):
                p.write_text('index\tvalue\n0\t'+value+'\n')
                with self.assertRaises(RuntimeError):rows(p,1)

    def test_incomplete_and_mutated_stage_refuse(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp);f=d/'input';f.write_text('original')
            with self.assertRaises(FileNotFoundError):verify(d)
            (d/'COMPLETE.json').write_text(json.dumps({'files':{'input':hashlib.sha256(f.read_bytes()).hexdigest()}}))
            verify(d)
            f.write_text('changed')
            with self.assertRaises(RuntimeError):verify(d)

if __name__=='__main__':unittest.main()
