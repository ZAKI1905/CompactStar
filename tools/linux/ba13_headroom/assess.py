#!/usr/bin/env python3
"""Superseded assessment entry: execute the corrected, isolated closure layer.

The old assessor's production-profile reference assignment was invalid (D1).
The historical outputs stay read-only under ba13-headroom. Regeneration now
targets the separate ba13-closure-20261019 artifact directory. See the committed
predeclaration; the historical fitter is never changed or called by this entry.
"""
import pathlib
import runpy

if __name__=='__main__':
    runpy.run_path(str(pathlib.Path(__file__).resolve().parents[1]/'ba13_closure/assess.py'),run_name='__main__')
