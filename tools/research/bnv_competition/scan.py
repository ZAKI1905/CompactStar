#!/usr/bin/env python3
"""Run independent S1 roots in a Slurm allocation; store every matched corner."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from model import Model, ROOT


def write_csv(path, data):
    keys = list(dict.fromkeys(k for row in data for k in row))
    with Path(path).open('w') as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(data)


def worker(task):
    drive, betas = task
    m = Model()
    critical = -2*m.I[0]*drive/m.t[0]/m.B*m.year
    rc = m.diagnostics(drive, 0.)
    out = []
    observables = ['Tinf_K', 'Tsurface_inf_K', 'eta_e_MeV', 'eta_mu_MeV',
                   'Lgamma_erg_s', 'LH_erg_s', 'DeltaLnu_erg_s', 'net_weak_erg_s']
    for beta in betas:
        gamma = float(beta)*critical
        both = m.diagnostics(drive, gamma, beta_e=float(beta))
        bnv = m.diagnostics(0., gamma)
        row = dict(both)
        row['Gamma_critical_e_yr'] = critical
        row['surface_ratio'] = both['Tsurface_inf_K']/rc['Tsurface_inf_K']
        row['luminosity_ratio'] = both['Lgamma_erg_s']/rc['Lgamma_erg_s']
        row['q_equals_RC_gross_keV'] = rc['LH_erg_s']/(m.erg*both['D_count_s'])*1000
        row['q_equals_RC_net_keV'] = rc['net_weak_erg_s']/(m.erg*both['D_count_s'])*1000
        row['q_offsets_luminosity_change_keV'] = abs(both['Lgamma_erg_s']-rc['Lgamma_erg_s'])/(m.erg*both['D_count_s'])*1000
        for key in observables:
            row['RC_'+key], row['BNV_'+key] = rc[key], bnv[key]
            row['difference_'+key] = both[key]-rc[key]
            row['interaction_'+key] = both[key]-rc[key]-bnv[key]  # 00=0.
            den = abs(rc[key])+abs(bnv[key])
            row['normalized_interaction_'+key] = row['interaction_'+key]/den if den else np.nan
        # Two-channel large-imbalance prediction, not the one-beta approximation.
        _, sp, bn, r = m.drivers(drive, gamma, beta_e=float(beta))
        e = np.sign(r)*(abs(r)*m.kberg/m.L/m.h[-1])**(1/7)*m.kb
        e0 = (-sp*m.kberg/m.L/m.h[-1])**(1/7)*m.kb
        row['asym_surface_ratio'] = (np.dot(e, r)/np.dot(e0, -sp))**.25
        row['asym_Lgamma_erg_s'] = .625*m.erg*np.dot(e, r)
        out.append((row, bnv))
    return rc, out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('out', type=Path)
    p.add_argument('--workers', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 1)))
    p.add_argument('--smoke', action='store_true')
    a = p.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'Scientific computations require the requested Slurm allocation'
    a.out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    m = Model()
    null_mu = m.I[1]*m.t[0]/(m.I[0]*m.t[1])
    betas = np.unique(np.r_[np.logspace(-3, 2, 101), np.linspace(.8, 1.4, 81),
                            1., null_mu, .1, .3, 2., 3., 10.])
    drives = np.logspace(-16, -9, 29)
    fiducial = 4*np.pi**2*1e-20/.005**3
    drives = np.unique(np.r_[drives, fiducial])
    if a.smoke:
        drives, betas = [fiducial], [.1, 1., null_mu, 3.]
    table, controls, bnv = [], [], []
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
        for rc, solutions in pool.map(worker, [(float(d), betas) for d in drives]):
            controls.append(rc)
            for row, single in solutions:
                table.append(row)
                bnv.append(single)
    controls.insert(0, m.diagnostics(0., 0.))
    write_csv(a.out/'stationary-full.csv', table)
    write_csv(a.out/'bnv-only.csv', bnv)
    write_csv(a.out/'controls.csv', controls)
    representative = [r for r in table if r['drive_rad2_s3'] == fiducial]
    write_csv(a.out/'fiducial.csv', representative)
    compact_keys = ['drive_rad2_s3','Gamma_yr','beta_e','beta_mu','Tinf_K','Tsurface_inf_K',
                    'eta_e_MeV','eta_mu_MeV','Lgamma_erg_s','LH_erg_s','DeltaLnu_erg_s',
                    'net_weak_erg_s','Lnu_eq_erg_s','A_B_erg_s','A_spin_erg_s','surface_ratio',
                    'normalized_interaction_Lgamma_erg_s','tau_slow_yr','build_time_yr',
                    'depletion_reach_estimate','temperature_in_certificate','reachable_in_certificate',
                    'q_equals_RC_net_keV','stable','chemical_residual','thermal_residual']
    write_csv(a.out/'map.csv', [{k:r[k] for k in compact_keys} for r in table])
    maximum = lambda key: max(abs(r[key]) for r in table+bnv)
    summary = dict(status='SCIENTIFIC_CANDIDATE', job=os.environ['SLURM_JOB_ID'],
                   host=platform.node(), workers=a.workers, wall_seconds=time.monotonic()-start,
                   python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                   combined_points=len(table), bnv_points=len(bnv), controls=len(controls),
                   all_stable=all(r['stable'] for r in table+bnv+controls),
                   max_chemical_residual=maximum('chemical_residual'),
                   max_thermal_residual=maximum('thermal_residual'),
                   max_subsystem_residual=maximum('subsystem_residual'),
                   muon_null_beta_e=null_mu,
                   within_temperature_gate=sum(r['temperature_in_certificate'] for r in table),
                   within_temperature_and_reach_depletion=sum(r['reachable_in_certificate'] for r in table),
                   minimum_fiducial=min(representative,key=lambda r:r['surface_ratio']),
                   input_manifest_sha256=hashlib.sha256((m.inputs/'manifest.json').read_bytes()).hexdigest(),
                   source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in Path(__file__).parent.glob('*.py')})
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['minimum_fiducial','source_sha256']},indent=2),flush=True)
    assert summary['all_stable']
    assert max(summary['max_chemical_residual'], summary['max_thermal_residual'], summary['max_subsystem_residual']) < 1e-9


if __name__ == '__main__':
    main()
