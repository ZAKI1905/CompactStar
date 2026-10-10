#!/usr/bin/env python3
"""Record candidate artifacts after native builds; never declare equivalence."""
import argparse
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'cmake'))
from verify_dependencies import digest, manifest
from linux_dependencies import SOURCES

p = argparse.ArgumentParser()
p.add_argument('--toolchain', type=Path, required=True)
p.add_argument('--sources', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
tc = a.toolchain.resolve()

def cache(path):
    return dict(line.split('=', 1) for line in path.read_text().splitlines()
                if '=' in line and not line.startswith(('#', '//')))

packages = {}
for mode in ('Debug', 'Release'):
    packages[mode] = {}
    for name, source_sha in SOURCES.items():
        source = a.sources/name
        head = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
        dirty = subprocess.check_output(['git', '-C', str(source), 'status', '--porcelain'], text=True)
        assert head == source_sha and not dirty, (name, head, dirty)
        prefix = tc/'install'/mode/name
        build = tc/'build'/(name+'-'+mode)
        c = cache(build/'CMakeCache.txt')
        record = dict(name=name, mode=mode, source_sha=head, dirty=False, prefix=str(prefix),
                      source=str(source.resolve()),
                      archive_sha256=digest(prefix/'lib'/('lib'+name+'.a')),
                      headers_sha256=manifest(prefix, (prefix/'include').rglob('*')),
                      configs_sha256=manifest(prefix, (prefix/'lib/cmake').rglob('*')),
                      cache_sha256=digest(build/'CMakeCache.txt'),
                      compile_commands_sha256=digest(build/'compile_commands.json'),
                      install_manifest_sha256=digest(build/'install_manifest.txt'),
                      cxx_flags=c['CMAKE_CXX_FLAGS:STRING'],
                      mode_flags=c['CMAKE_CXX_FLAGS_'+mode.upper()+':STRING'],
                      cxx_standard=17,
                      ctest_exit=int((tc/'logs'/(name+'-'+mode+'-tests.exit')).read_text()))
        (prefix/'provenance.json').write_text(json.dumps(record, indent=2)+'\n')
        packages[mode][name] = record
compiler_root = Path('/opt/spack/linux-rocky9-zen2/gcc-11.4.1/gcc-14.1.0-yd5t5aj4bggpy5swlajrqvsiohzyyfbh')
roots = {'gsl': tc/'install/gsl-2.7.1', 'sysroot': tc/'sysroot',
         'compiler': compiler_root, 'tools': tc/'bin'}
closure = {}
for name, root in roots.items():
    paths = sorted(p for p in root.rglob('*') if p.is_file())
    inventory = ''.join(f'{digest(p)}  {p.relative_to(root)}\n' for p in paths)
    (tc/(name+'-files.sha256')).write_text(inventory)
    closure[name] = {'root': str(root), 'sha256': hashlib.sha256(inventory.encode()).hexdigest()}
lock = dict(schema='compactstar-linux-dependencies-v1', status='CANDIDATE_UNQUALIFIED',
            platform={'system': platform.system(), 'machine': platform.machine()},
            gsl_source_sha256='dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b',
            toolchain={'compiler': str(tc/'bin/g++'), 'compiler_sha256': digest(tc/'bin/g++'),
                       'compiler_version': '14.1.0', 'cxx_flags': '-ffp-contract=off',
                       'mode_flags': {'Debug': '-g', 'Release': '-O3 -DNDEBUG'}},
            link_inputs={'gsl_include': str(roots['gsl']/'include'),
                         'gsl_library': str(roots['gsl']/'lib/libgsl.a'),
                         'gsl_cblas': str(roots['gsl']/'lib/libgslcblas.a'),
                         'zlib_include': str(tc/'sysroot/usr/include'),
                         'zlib_library': str(tc/'sysroot/usr/lib64/libz.so')},
            closure=closure, packages=packages)
a.output.write_text(json.dumps(lock, indent=2)+'\n')
print('Recorded candidate lock:', a.output, digest(a.output))
