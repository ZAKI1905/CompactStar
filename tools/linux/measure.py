#!/usr/bin/env python3
"""Run one command and report Linux child resource usage without GNU time."""
import json
import resource
import subprocess
import sys
import time

started = time.monotonic()
result = subprocess.run(sys.argv[1:])
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
print('RESOURCE_USAGE_JSON ' + json.dumps(dict(
    command=sys.argv[1:], exit_code=result.returncode,
    elapsed_seconds=time.monotonic()-started,
    user_cpu_seconds=usage.ru_utime, system_cpu_seconds=usage.ru_stime,
    max_rss_kib=usage.ru_maxrss,
    note='Linux RUSAGE_CHILDREN; maximum child RSS, not simultaneous aggregate RSS')),
    file=sys.stderr)
sys.exit(result.returncode)
