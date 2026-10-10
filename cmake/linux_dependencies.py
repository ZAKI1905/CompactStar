"""Fail-closed authentication of a separately built native Linux candidate."""
import json
import platform
import subprocess
from pathlib import Path

SOURCES = {'Zaki': 'e263a6e180c5c417198e7778bd21fc9c0a32dc33',
           'CONFIND': 'b0cbd510fd3fd0c772fa50499cd749287cb39e7b'}


def check_linux(args, digest, manifest):
    def require(condition, message):
        if not condition:
            raise SystemExit('Linux dependency authentication: ' + message)

    require(platform.system() == 'Linux', 'native Linux required')
    require(args.linux_lock and args.linux_lock.is_absolute(), 'absolute Linux lock required')
    lock = json.loads(args.linux_lock.read_text())
    require(lock['schema'] == 'compactstar-linux-dependencies-v1', 'unknown lock schema')
    require(lock['platform'] == {'system': 'Linux', 'machine': platform.machine()}, 'platform mismatch')
    require(args.system == 'Linux' and args.processor == platform.machine(), 'cross compilation unsupported')
    require(lock['status'] == 'CANDIDATE_UNQUALIFIED', 'lock cannot confer scientific qualification')
    tool = lock['toolchain']
    require(args.compiler == tool['compiler'], 'compiler path mismatch')
    require(digest(Path(args.compiler)) == tool['compiler_sha256'], 'compiler bytes mismatch')
    require(args.compiler_version == tool['compiler_version'], 'compiler version mismatch')
    require(set(tool['executables']) == {'cmake', 'as', 'ld', 'ar', 'ranlib', 'git'}, 'tool inventory mismatch')
    for name, executable in tool['executables'].items():
        require(digest(Path(executable['path'])) == executable['sha256'], name + ' executable mismatch')
    require(args.cxx_flags == tool['cxx_flags'], 'compiler flags mismatch')
    require(args.mode_flags == tool['mode_flags'][args.mode], 'mode flags mismatch')
    require(lock['gsl_source_sha256'] == 'dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b', 'GSL source mismatch')
    require(set(lock['closure']) == {'gsl', 'sysroot', 'compiler', 'tools'}, 'unknown/missing closure')
    for name, item in lock['closure'].items():
        root = Path(item['root'])
        require(root.is_absolute() and root.is_dir(), name + ' root missing')
        require(manifest(root, root.rglob('*')) == item['sha256'], name + ' closure mismatch')
    requested = {'gsl_include': args.gsl_include, 'gsl_library': args.gsl_library,
                 'gsl_cblas': args.gsl_cblas, 'zlib_include': args.zlib_include,
                 'zlib_library': args.zlib_library}
    require(requested == lock['link_inputs'], 'GSL/zlib selection mismatch')
    require(set(lock['packages']) == {'Debug', 'Release'}, 'mode inventory mismatch')
    expected_mode = lock['packages'][args.mode]
    require(set(expected_mode) == set(SOURCES), 'package inventory mismatch')
    records = []
    for name, supplied in [('Zaki', args.zaki), ('CONFIND', args.confind)]:
        expected = expected_mode[name]
        require(supplied.is_absolute(), name + ' prefix must be absolute')
        require(str(supplied.resolve(strict=True)) == expected['prefix'], name + ' prefix mismatch')
        require(expected['source_sha'] == SOURCES[name] and expected['mode'] == args.mode
                and expected['dirty'] is False, name + ' source/mode mismatch')
        prefix = supplied.resolve()
        actual = {'archive_sha256': digest(prefix/'lib'/('lib'+name+'.a')),
                  'headers_sha256': manifest(prefix, (prefix/'include').rglob('*')),
                  'configs_sha256': manifest(prefix, (prefix/'lib/cmake').rglob('*'))}
        require(all(expected[k] == v for k, v in actual.items()), name + ' artifact mismatch')
        require(json.loads((prefix/'provenance.json').read_text()) == expected,
                name + ' provenance mismatch')
        records.append(expected)
    return {'platform': lock['platform'], 'lock_sha256': digest(args.linux_lock),
            'status': lock['status'], 'packages': records, 'link_inputs': requested,
            'toolchain': tool, 'closure': lock['closure']}
