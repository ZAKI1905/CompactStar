"""Run one Phase-5B gate with isolated, preserved output on every invocation."""
import argparse
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("executable")
parser.add_argument("stage")
parser.add_argument("output_root")
args = parser.parse_args()
root = Path(args.output_root)
root.mkdir(parents=True, exist_ok=True)
run = Path(tempfile.mkdtemp(prefix=args.stage + "-", dir=root))
print(f"Phase-5B isolated evidence: {run}", flush=True)
result = subprocess.run([args.executable, args.stage, str(run / "data")],
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
(run / "run.log").write_text(result.stdout)
print(result.stdout, end="")
raise SystemExit(result.returncode)
