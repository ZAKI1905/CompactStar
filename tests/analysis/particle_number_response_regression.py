"""Regenerate the governed Phase-5B response and compare exact artifact bytes."""

import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile


parser = argparse.ArgumentParser()
parser.add_argument("--producer", required=True, type=Path)
parser.add_argument("--executable", required=True, type=Path)
parser.add_argument("--source-root", required=True, type=Path)
parser.add_argument("--baseline", required=True, type=Path)
parser.add_argument("--output-root", required=True, type=Path)
args = parser.parse_args()

output_root = args.output_root.resolve()
output_root.mkdir(parents=True, exist_ok=True)
run_root = Path(tempfile.mkdtemp(prefix="regression-", dir=output_root))
generated_root = run_root / "generated"
command = [
    sys.executable,
    str(args.producer.resolve()),
    "--executable",
    str(args.executable.resolve()),
    "--source-root",
    str(args.source_root.resolve()),
    "--output-dir",
    str(generated_root),
]
result = subprocess.run(command, capture_output=True, text=True)
(run_root / "producer.log").write_text(result.stdout + result.stderr)
print(result.stdout, end="")
if result.returncode:
    raise SystemExit(f"Phase-5B producer failed; see {run_root / 'producer.log'}")

generated = (generated_root / "structural-response.json").resolve()
baseline = args.baseline.resolve()
if generated == baseline:
    raise SystemExit("Refusing a self-comparison: generated artifact is the governed baseline")

generated_bytes = generated.read_bytes()
baseline_bytes = baseline.read_bytes()
generated_sha = hashlib.sha256(generated_bytes).hexdigest()
baseline_sha = hashlib.sha256(baseline_bytes).hexdigest()
print(f"generated_sha256={generated_sha}")
print(f"baseline_sha256={baseline_sha}")
print(f"isolated_evidence={run_root}")
if generated_bytes != baseline_bytes:
    raise SystemExit("Fresh Phase-5B structural response differs from governed baseline")
