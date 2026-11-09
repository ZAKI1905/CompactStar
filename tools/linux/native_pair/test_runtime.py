#!/usr/bin/env python3
"""Focused dependency-discovery and refusal tests; no numerical executable runs."""
import os
import unittest
from unittest.mock import patch
from runtime import interpreter_from_headers, resolved_paths, runtime_binding, verify_runtime


class RuntimeTests(unittest.TestCase):
    def test_actual_interpreter_not_host_ldd_alias(self):
        self.assertEqual(interpreter_from_headers('[Requesting program interpreter: /shared/ld.so]'), '/shared/ld.so')
        self.assertEqual(resolved_paths('linux-vdso.so.1 (0x12)\nlibc.so.6 => /shared/libc.so.6 (0x34)\n/shared/ld.so (0x56)'),
                         ['/shared/ld.so', '/shared/libc.so.6'])
        with self.assertRaises(ValueError):
            resolved_paths('/shared/ld.so => /lib64/ld.so (0x56)')

    def test_malformed_or_missing_linkage_refuses(self):
        for headers in ('', '[Requesting program interpreter: relative]',
                        '[Requesting program interpreter: /a]\n[Requesting program interpreter: /b]'):
            with self.assertRaises(ValueError):
                interpreter_from_headers(headers)
        for listing in ('', 'libc.so.6 => not found', 'unexpected line'):
            with self.assertRaises(ValueError):
                resolved_paths(listing)

    def test_injection_refuses_before_discovery(self):
        for key in ('LD_PRELOAD', 'LD_AUDIT'):
            with patch.dict(os.environ, {key: '/unbound/library.so'}):
                with self.assertRaises(ValueError):
                    runtime_binding('/unused')

    def test_changed_resolution_and_unbound_loader_refuse(self):
        expected = {'interpreter': '/shared/ld.so', 'dependencies': [
            {'path': '/shared/ld.so', 'resolved': '/shared/ld.so'}]}
        with patch('runtime.runtime_binding', return_value=expected):
            verify_runtime({'/binary': expected}, {'/binary': 'hash', '/shared/ld.so': 'hash'})
            with self.assertRaises(ValueError):
                verify_runtime({'/binary': expected}, {'/binary': 'hash'})
            with self.assertRaises(ValueError):
                verify_runtime({'/binary': {}}, {'/binary': 'hash', '/shared/ld.so': 'hash'})
            with self.assertRaises(ValueError):
                verify_runtime({'/binary': expected}, {})


if __name__ == '__main__':
    unittest.main()
