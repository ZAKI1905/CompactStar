#!/usr/bin/env python3
"""Test an inert restart-selection proposal; never call a trajectory/RHS.

Run on a compute allocation. Applies the proposed patch only to a fresh scratch
repository containing two harness files, then compiles its standalone selector
test. The active runner and runtime snapshots are not modified.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from decimal import Decimal, localcontext

ROOT = Path(__file__).resolve().parents[3]
DOC = ROOT / 'docs/validation/results-first-bnv-pilot'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out):
    out.mkdir(parents=True, exist_ok=False)
    result = dict(status='FAIL', trajectory_calls=0, proposal_applied_to_active_runner=False)
    try:
        affinity = sorted(os.sched_getaffinity(0))
        assert len(affinity) == 1, affinity
        result['affinity'] = affinity
        patch = DOC / 'restart-coverage/proposed-selector.patch'
        old = DOC / 'stage-a-recovery/evidence'
        history = old / 'main.accepted_states.tsv'
        proof = json.loads((old / 'restart-gate.json').read_text())
        assert proof['status'] == 'FAIL_REQUIRED_RESTART_COVERAGE'
        assert proof['quarter_cap_executed'] is False
        for name, digest in proof['hashes'].items():
            assert sha(Path(name)) == digest, name
        rows = list(csv.DictReader(history.open(), delimiter='\t'))
        half, end = 15778800., 31557600.
        original = next(r for r in rows if float(r['t_after_s']) >= half)
        proposed = [r for r in rows if 0 < float(r['t_after_s']) <= half][-1]
        assert float(original['t_after_s']) == end
        assert proposed['sequence'] == '11'
        packet_path = Path(next(p for p in proof['hashes'] if p.endswith('/main.accepted_states.tsv') and '/endpoint-schedule-repeat/' not in p)).parent / 'checkpoints-240/point-120/checkpoint.json'
        packet = json.loads(packet_path.read_text())
        assert packet['time_s'] == half and packet['main_left_time_s'] < half < packet['main_right_time_s']
        assert 'reconstructed' in packet['main_at_checkpoint']
        shutil.copyfile(packet_path, out / 'midpoint-checkpoint.json')
        with localcontext() as ctx:
            ctx.prec = 80
            t = Decimal.from_float(float(proposed['t_after_s']))
            duration = Decimal.from_float(end) - t
            result['selection'] = dict(original_ordinal=12, original_tail_s='0',
                proposed_ordinal=11, time_s=str(t), time_binary64_hex=float(t).hex(),
                tail_s=str(duration), tail_fraction=str(duration/Decimal.from_float(end)),
                state_exact_binary64=[str(Decimal.from_float(float(proposed[k]))) for k in ('x_state','eta_e_MeV','eta_mu_MeV')])
        result['quarter_cap'] = dict(executed=False, coverage_inferred_from_authenticated_inactive_cap_proof=True,
            proof_sha256=sha(old/'restart-gate.json'), minimum_cap_s=proof['minimum_planned_cap_s'])
        result['midpoint'] = dict(packet_index=120, accepted_endpoint=False,
            method=packet['method'], bracket=[packet['main_left_time_s'],packet['main_right_time_s']])
        scratch = out / 'candidate'
        scratch.mkdir()
        subprocess.run(['git','init','-q',str(scratch)],check=True)
        for name in ('tests/native/stage.cpp','tests/CMakeLists.txt'):
            target = scratch/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/name,target)
        subprocess.run(['git','-C',str(scratch),'apply','--check',str(patch)],check=True)
        subprocess.run(['git','-C',str(scratch),'apply',str(patch)],check=True)
        compiler = os.environ['CXX']
        command = [compiler,'-std=c++17','-O2','-ffp-contract=off','-Wall','-Wextra','-Werror',
                   str(scratch/'tests/native/restart_schedule_test.cpp'),'-o',str(out/'selector-test')]
        subprocess.run(command,check=True)
        tested = subprocess.run([str(out/'selector-test'),str(history)],check=True,capture_output=True,text=True)
        (out/'tests.log').write_text(tested.stdout+tested.stderr)
        result['test_command'] = command
        result['compiler_version'] = subprocess.check_output([compiler,'--version'],text=True)
        result['test_result'] = tested.stdout
        result['scope'] = 'Pure C++ selector/refusal/serialization tests; patched stage executable not built or run. No accuracy qualification.'
        result['input_hashes'] = {str(p):sha(p) for p in (Path(__file__),patch,history,old/'restart-gate.json',packet_path,ROOT/'tests/native/stage.cpp',ROOT/'tests/CMakeLists.txt')}
        result['candidate_hashes'] = {str(p.relative_to(scratch)):sha(p) for p in sorted((scratch/'tests').rglob('*')) if p.is_file()}
        result['status'] = 'PASS_PROPOSED_SELECTION_ONLY'
    except Exception as e:
        result['failure'] = repr(e)
        raise
    finally:
        (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fresh_output',type=Path)
    main(parser.parse_args().fresh_output)
