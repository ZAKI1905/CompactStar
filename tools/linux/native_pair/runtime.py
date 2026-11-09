"""Authenticate native ELF linkage through PT_INTERP, never the host ldd wrapper."""
import os
import pathlib
import re
import subprocess


def interpreter_from_headers(headers):
    matches = re.findall(r"\[Requesting program interpreter: ([^\]\n]+)\]", headers)
    if len(matches) != 1 or not matches[0].startswith('/'):
        raise ValueError('missing/ambiguous absolute ELF interpreter')
    return matches[0]


def resolved_paths(listing):
    paths = []
    for line in listing.splitlines():
        line = line.strip()
        if not line:
            continue
        if re.fullmatch(r'linux-vdso\.so\.1 \(0x[0-9a-fA-F]+\)', line):
            continue  # Kernel mapping, not an on-disk shared library.
        match = re.fullmatch(r'(?:[^\s/]+ => )?(/[^\s]+) \(0x[0-9a-fA-F]+\)', line)
        if not match:
            raise ValueError('unresolved/unrecognized ELF dependency: ' + line)
        paths.append(match[1])
    if not paths:
        raise ValueError('empty ELF dependency listing')
    return sorted(set(paths))


def runtime_binding(binary):
    # Injected libraries/auditors are outside the declared native experiment.
    for key in ('LD_PRELOAD', 'LD_AUDIT'):
        if os.environ.get(key):
            raise ValueError('undeclared runtime injection: ' + key)
    preload = pathlib.Path('/etc/ld.so.preload')
    if preload.exists() and any(line.split('#', 1)[0].strip()
                                for line in preload.read_text().splitlines()):
        raise ValueError('undeclared system runtime preload')
    headers = subprocess.check_output(['readelf', '-lW', str(binary)], text=True)
    interpreter = interpreter_from_headers(headers)
    listing = subprocess.check_output([interpreter, '--list', str(binary)], text=True)
    paths = resolved_paths(listing)
    paths.append(interpreter)
    dependencies = []
    for path in sorted(set(paths)):
        resolved = pathlib.Path(path).resolve(strict=True)
        if not resolved.is_file():
            raise ValueError('runtime dependency is not a file: ' + path)
        dependencies.append({'path': path, 'resolved': str(resolved)})
    return {'interpreter': interpreter, 'dependencies': dependencies}


def verify_runtime(bindings, bound):
    for binary, expected in bindings.items():
        if binary not in bound:
            raise ValueError('unbound runtime executable: ' + binary)
        actual = runtime_binding(binary)
        if actual != expected:
            raise ValueError('runtime dependency resolution mismatch: ' + binary)
        for dependency in actual['dependencies']:
            if dependency['resolved'] not in bound:
                raise ValueError('unbound runtime dependency: ' + dependency['resolved'])
