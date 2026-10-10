#!/usr/bin/env python3
"""Stage installed RPM development/runtime bytes without modifying the host.

Only byte-verified regular files and relative symlinks are copied. This is an
installed-RPM integrity check, not an independent package-signature attestation.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=False)
records=[]
for pkg in ('glibc','glibc-devel','glibc-headers','kernel-headers','zlib','zlib-devel'):
    version=subprocess.check_output(['rpm','-q',pkg],text=True).strip()
    result=subprocess.run(['rpm','-V',pkg],capture_output=True,text=True)
    (a.output/(pkg+'-rpm-verify.log')).write_text(result.stdout+result.stderr)
    if result.stderr:raise RuntimeError(result.stderr)
    for line in result.stdout.splitlines():
        if line[:5]!='.....' or line[7:9]!='..':raise RuntimeError('RPM content verification: '+line)
    for n in subprocess.check_output(['rpm','-ql',pkg],text=True).splitlines():
        src=Path(n)
        if not n.startswith(('/usr/include/','/usr/lib64/')) or src.is_dir():continue
        dst=a.output/'sysroot'/n.lstrip('/');dst.parent.mkdir(parents=True,exist_ok=True)
        if src.is_symlink():
            target=src.readlink()
            if target.is_absolute():raise RuntimeError('Absolute symlink: '+n)
            if not dst.is_symlink():dst.symlink_to(target)
            records.append(dict(package=version,path=n,symlink=str(target)))
        elif src.is_file():
            shutil.copy2(src,dst)
            records.append(dict(package=version,path=n,sha256=hashlib.sha256(dst.read_bytes()).hexdigest()))
(a.output/'sysroot/lib64').symlink_to('usr/lib64')
(a.output/'sysroot-manifest.json').write_text(json.dumps(records,indent=2)+'\n')
