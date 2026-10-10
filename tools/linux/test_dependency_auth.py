#!/usr/bin/env python3
"""Exercise refusal of substituted source, platform, closure and package bytes."""
import copy
import json
import platform
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'cmake'))
from verify_dependencies import digest, manifest
from linux_dependencies import check_linux, SOURCES

class Authentication(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.args = SimpleNamespace(linux_lock=root/'lock.json', mode='Debug', system='Linux',
            processor=platform.machine(), compiler=str(root/'compiler'), compiler_version='14.1.0',
            cxx_flags='-ffp-contract=off', mode_flags='-g', gsl_include='gsl/include',
            gsl_library='gsl/a', gsl_cblas='gsl/b', zlib_include='z/include', zlib_library='z/a')
        Path(self.args.compiler).write_text('compiler')
        closure = {}
        for key in ('gsl', 'sysroot', 'compiler', 'tools'):
            path = root/('closure-'+key); path.mkdir(); (path/'file').write_text(key)
            closure[key] = dict(root=str(path), sha256=manifest(path,path.rglob('*')))
        packages = {}
        for name, sha in SOURCES.items():
            prefix=root/name
            for sub in ('lib/cmake','include'): (prefix/sub).mkdir(parents=True)
            (prefix/'lib'/('lib'+name+'.a')).write_text('archive')
            (prefix/'include/header').write_text('header')
            (prefix/'lib/cmake/config').write_text('config')
            rec=dict(prefix=str(prefix),source_sha=sha,mode='Debug',dirty=False,
                archive_sha256=digest(prefix/'lib'/('lib'+name+'.a')),
                headers_sha256=manifest(prefix,(prefix/'include').rglob('*')),
                configs_sha256=manifest(prefix,(prefix/'lib/cmake').rglob('*')))
            (prefix/'provenance.json').write_text(json.dumps(rec)); packages[name]=rec
            setattr(self.args,name.lower(),prefix)
        self.lock=dict(schema='compactstar-linux-dependencies-v1',status='CANDIDATE_UNQUALIFIED',
            platform=dict(system='Linux',machine=platform.machine()),
            toolchain=dict(compiler=self.args.compiler,compiler_sha256=digest(Path(self.args.compiler)),
                executables={n:dict(path=self.args.compiler,sha256=digest(Path(self.args.compiler))) for n in ('cmake','as','ld','ar','ranlib','git')},compiler_version='14.1.0',cxx_flags=self.args.cxx_flags,mode_flags={'Debug':'-g','Release':'-O3 -DNDEBUG'}),
            gsl_source_sha256='dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b',
            closure=closure,packages={'Debug':packages,'Release':{}},
            link_inputs={k:getattr(self.args,k) for k in ('gsl_include','gsl_library','gsl_cblas','zlib_include','zlib_library')})
    def run_check(self):
        if self.args.linux_lock: self.args.linux_lock.write_text(json.dumps(self.lock))
        return check_linux(self.args,digest,manifest)
    def test_valid(self): self.assertEqual(self.run_check()['status'],'CANDIDATE_UNQUALIFIED')
    def test_selection_refusals(self):
        for key,value in [('linux_lock',None),('system','Darwin'),('processor','arm64'),
                          ('compiler_version','11.5.0'),('cxx_flags','-ffast-math'),
                          ('mode_flags','-O3'),('gsl_library','other'),('mode','Release')]:
            old=getattr(self.args,key)
            with self.subTest(key=key),self.assertRaises(SystemExit):
                setattr(self.args,key,value); self.run_check()
            setattr(self.args,key,old)
    def test_changed_artifacts(self):
        paths=[self.args.zaki/'include/header',self.args.confind/'lib/cmake/config',
               self.args.zaki/'lib/libZaki.a',self.args.zaki/'provenance.json',
               Path(self.lock['closure']['gsl']['root'])/'file',Path(self.args.compiler)]
        for p in paths:
            old=p.read_bytes(); p.write_bytes(old.replace(b'e263a6e', b'0000000') if p.name == 'provenance.json' else old+b' ')
            with self.subTest(path=p),self.assertRaises(SystemExit): self.run_check()
            p.write_bytes(old)
    def test_wrong_source_and_inventory(self):
        for mutate in (lambda x:x['packages']['Debug']['Zaki'].update(source_sha='0'*40),
                       lambda x:x['packages']['Debug'].update(unknown={}),
                       lambda x:x['closure'].pop('sysroot'),
                       lambda x:x['platform'].update(system='Darwin')):
            old=copy.deepcopy(self.lock);mutate(self.lock)
            with self.assertRaises(SystemExit): self.run_check()
            self.lock=old
    def test_nonlinux(self):
        with patch('linux_dependencies.platform.system',return_value='Darwin'),self.assertRaises(SystemExit): self.run_check()
    def test_mac_lock_unchanged(self):
        self.assertEqual(digest(ROOT/'cmake/dependency-lock.json'),
            'e1307a042180034ad028f2ff07dbc7fd80752494ff4b2ff61cc8ca0eec2089c5')

if __name__ == '__main__': unittest.main()
