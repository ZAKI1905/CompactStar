"""ADR-0025 S1: independent frozen Structure-1 stationary model (candidate).

Infinity time/temperature/energies. Counts are dimensionless. No production
integration, certificate promotion, microscopic source, or finite-T BNV model.
"""
from pathlib import Path
import csv
import hashlib
import json
import math

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / 'docs/research/bnv-rotochemical-competition/inputs'


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


class Model:
    def __init__(self, inputs=INPUTS):
        self.inputs = Path(inputs)
        for name, item in json.loads((self.inputs/'manifest.json').read_text()).items():
            assert hashlib.sha256((self.inputs/name).read_bytes()).hexdigest() == item['sha256'], name
        c = {k: float(v) for k, v in rows(self.inputs/'coefficients.tsv')[0].items()}
        t = {k: float(v) for k, v in rows(self.inputs/'tangent.tsv')[0].items()}
        a = {r['quantity']: float(r['value']) for r in rows(self.inputs/'anchors.tsv')}
        self.structure = json.loads((self.inputs/'structure.json').read_text())
        self.Z = np.array([[c['Z00'], c['Z01']], [c['Z10'], c['Z11']]])
        self.Zinv = np.linalg.inv(self.Z)
        self.L = np.array([c['Ltilde_Me'], c['Ltilde_Mmu']])
        self.t = np.array([t['t_e'], t['t_mu']])
        self.I = np.array(self.structure['I_count_s2'])
        self.N = np.array(self.structure['N_leptons'])
        self.B, self.mu = c['B_solved'], c['mu_B_inf']
        self.kb, self.erg = self.structure['kb_MeV_K'], self.structure['MeV_erg']
        self.kberg = self.kb*self.erg
        self.year = self.structure['year_s']
        self.Lg0, self.C0 = a['Lgamma'], a['Cstar']
        self.surface_area = 4*math.pi*(self.structure['R_km']*1e5)**2
        self.redshift = math.exp(self.structure['nu_surface'])
        self.sigmaSB = self.structure['sigmaSB_cgs']
        self.power_exponent = 2.42  # Actual production Iron envelope, fixed geometry.
        cap = rows(self.inputs/'capacity.tsv')
        self.logT = np.log([float(r['T_MeV'])/self.kb for r in cap])
        self.cap = np.array([float(r['C_erg_K']) for r in cap])
        self.h = np.array([14680, 7560, 840, 24])/11513/np.pi**np.array([2, 4, 6, 8])
        self.f = np.array([22020, 5670, 420, 9])/11513/np.pi**np.array([2, 4, 6, 8])
        self.odd = np.array([1, 3, 5, 7])
        self.even = self.odd+1
        assert abs(t['t_n']+sum(self.t)-1) < 3e-16
        assert np.array_equal(self.Z, self.Z.T) and min(np.linalg.eigvalsh(self.Z)) > 0

    def capacity(self, T):
        if self.logT[0] <= math.log(T) <= self.logT[-1]:
            return float(np.interp(math.log(T), self.logT, self.cap))
        # Times only: explicitly declared cold-degenerate extrapolation, not used
        # to determine a stationary root and never admitted by production.
        return self.C0*T/1e8

    def reactions(self, T, eta):
        eta = np.asarray(eta)
        xi = eta/(self.kb*T)
        H = np.sum(self.h*xi[:, None]**self.odd, axis=1)
        Fm1 = np.sum(self.f*xi[:, None]**self.even, axis=1)
        R = self.L/self.kberg*T**7*H
        return R, float(self.erg*np.dot(eta, R)), float(np.sum(self.L*T**8*Fm1)), float(sum(self.L)*T**8)

    def photon(self, T):
        return self.Lg0*(T/1e8)**self.power_exponent

    def surface(self, T):
        return (self.photon(T)*self.redshift**2/(self.surface_area*self.sigmaSB))**.25

    def drivers(self, omega_abs_omegadot, gamma_yr, beta_e=None):
        D = gamma_yr*self.B/self.year
        spin = 2*self.I*omega_abs_omegadot  # s_Omega = -2 I Omega Omegadot < 0.
        bnv = D*self.t
        if beta_e is not None:
            # The scan prescribes beta, including exact mathematical nulls.
            # Factoring 1-beta avoids subtracting ~1e36 count/s at a null.
            # Gamma remains the equivalent rounded display coordinate.
            null_mu = self.I[1]*self.t[0]/(self.I[0]*self.t[1])
            beta = np.array([beta_e, beta_e/null_mu])
            bnv = -spin*beta
            return D, spin, bnv, -spin*(1-beta)
        return D, spin, bnv, -(spin+bnv)  # Required R at stationarity.

    def eta_at(self, T, target, xtol=2e-13):
        values = []
        for j, rate in enumerate(target):
            if rate == 0:
                values.append(0.0)
                continue
            # Monotone odd polynomial: exactly one signed root in each channel.
            scale = self.L[j]/self.kberg*T**7
            def equation(logxi):
                xi = 10**logxi
                return math.log10(scale*np.dot(self.h, xi**self.odd))-math.log10(abs(rate))
            logxi = brentq(equation, -40, 12, xtol=xtol)
            values.append(math.copysign(self.kb*T*10**logxi, rate))
        return np.array(values)

    def solve(self, drive, gamma, xtol=2e-12, beta_e=None):
        _, spin, bnv, target = self.drivers(drive, gamma, beta_e)
        if not np.any(target):
            return 0.0, np.zeros(2)  # Only boundary equilibrium; no positive-T root.
        def balance(logT):
            T = 10**logT
            eta = self.eta_at(T, target, xtol/10)
            _, heat, excess, equilibrium = self.reactions(T, eta)
            photon = self.photon(T)
            return (heat-excess-equilibrium-photon)/max(heat, excess+equilibrium+photon)
        logT = brentq(balance, 0, 10, xtol=xtol)
        T = 10**logT
        return T, self.eta_at(T, target, xtol/10)

    def jacobian(self, T, eta):
        """Full T,eta vector-field Jacobian, evaluated at a stationary solution.

        d(C^-1)/dT drops out because Pnet=0. No chemical-only stability claim.
        """
        eta = np.asarray(eta)
        xi = eta/(self.kb*T)
        H = np.sum(self.h*xi[:, None]**self.odd, axis=1)
        Hp = np.sum(self.odd*self.h*xi[:, None]**(self.odd-1), axis=1)
        F = 1+np.sum(self.f*xi[:, None]**self.even, axis=1)
        Fp = np.sum(self.even*self.f*xi[:, None]**(self.even-1), axis=1)
        R = self.L/self.kberg*T**7*H
        Re = self.L/self.kberg/self.kb*T**6*Hp
        Rt = self.L/self.kberg*T**6*(7*H-xi*Hp)
        nT = self.L*T**7*(8*F-xi*Fp)
        ne = self.L/self.kb*T**7*Fp
        J = np.zeros((3, 3))
        J[0, 0] = (self.erg*np.dot(eta, Rt)-sum(nT)-self.power_exponent*self.photon(T)/T)/self.capacity(T)
        J[0, 1:] = (self.erg*(R+eta*Re)-ne)/self.capacity(T)
        J[1:, 0] = -self.Z@Rt
        J[1:, 1:] = -self.Z@np.diag(Re)
        return J

    def diagnostics(self, drive, gamma, P_s=.005, beta_e=None):
        D, spin, bnv, target = self.drivers(drive, gamma, beta_e)
        T, eta = self.solve(drive, gamma, beta_e=beta_e)
        omega = 2*np.pi/P_s
        beta = np.divide(bnv, -spin, out=np.full(2, np.nan), where=spin != 0)
        rotation = abs(self.I)*omega**2/self.N if drive else np.zeros(2)
        # Newtonian dynamical frequency is a conservative explicitly named
        # admission estimate, not a relativistic mass-shedding certificate.
        omega_K = 299792.458*math.sqrt(self.structure['M_km']/self.structure['R_km']**3)
        d = dict(drive_rad2_s3=drive, Gamma_yr=gamma, P_s=P_s, Pdot=drive*P_s**3/(4*np.pi**2),
                 beta_e=beta[0], beta_mu=beta[1], D_count_s=D,
                 Tinf_K=T, Tsurface_inf_K=self.surface(T) if T else 0.,
                 eta_e_MeV=eta[0], eta_mu_MeV=eta[1], q_dep_MeV=0., P_direct_erg_s=0.,
                 rotation_e=rotation[0], rotation_mu=rotation[1], omega_over_dynamic_K=omega/omega_K if drive else 0.)
        if not T:
            d.update(Lgamma_erg_s=0., LH_erg_s=0., DeltaLnu_erg_s=0., net_weak_erg_s=0.,
                     Lnu_eq_erg_s=0., A_B_erg_s=0., A_spin_erg_s=0., Echem_erg=0.,
                     stable=True, stability='boundary nonlinear attracting; zero linear rates',
                     tau_slow_yr=math.inf, build_time_yr=math.inf, depletion_build=0.,
                     thermal_residual=0., chemical_residual=0., subsystem_residual=0.,
                     temperature_in_certificate=False, reachable_in_certificate=False)
            return d
        R, LH, excess, eq = self.reactions(T, eta)
        photon = self.photon(T)
        AB, AO = -self.erg*np.dot(eta, bnv), -self.erg*np.dot(eta, spin)
        J = self.jacobian(T, eta)
        # Balance coordinates to avoid unit-induced numerical conditioning.
        S = np.array([T, .01, .01])
        balanced = J*S[None, :]/S[:, None]
        eig = np.linalg.eigvals(balanced)
        eigen_precision = 'binary64'
        if min(abs(eig)) < 1e-10*max(abs(eig)):
            # Exact channel cancellation can separate the thermal and chemical
            # eigenvalues by >20 decades. Resolve the slow eigenvalue instead
            # of interpreting a rounded zero as instability.
            import mpmath as mp
            with mp.workdps(55):
                eig = np.array([complex(v) for v in mp.eig(mp.matrix(balanced.tolist()),left=False,right=False)])
            eigen_precision = '55-digit eigenanalysis of binary64 Jacobian'
        stable = bool(max(eig.real) < 0)
        tau = -1/max(eig.real)/self.year if stable else math.inf
        # A lower bound on attaining eta from eta=0 for constant drives: the
        # chemical reservoir can never grow faster than source-only work.
        # Per-component initial-slope estimates can diverge when Zs cancels;
        # retain them explicitly, and use the energy-metric norm estimate too.
        initial = self.Z@target
        with np.errstate(divide='ignore', invalid='ignore'):
            build = np.abs(eta/initial)/self.year
        energy_build = math.sqrt(float(eta@self.Zinv@eta)/float(target@self.Z@target))/self.year
        # Constant dipole PPdot -> |d ln(Omega Omegadot)/dt|^-1 = Omega²/(4 drive).
        tdrive = omega**2/(4*drive)/self.year if drive else math.inf
        ell = -self.Zinv@eta
        displacement = abs(ell)/self.N
        chem_scale = max(np.max(abs(spin)), np.max(abs(bnv)), np.max(abs(R)), 1.)
        d.update(xi_e=eta[0]/(self.kb*T), xi_mu=eta[1]/(self.kb*T),
                 R_e_count_s=R[0], R_mu_count_s=R[1],
                 Lgamma_erg_s=photon, LH_erg_s=LH, DeltaLnu_erg_s=excess,
                 net_weak_erg_s=LH-excess, Lnu_eq_erg_s=eq,
                 A_B_erg_s=AB, A_spin_erg_s=AO, Echem_erg=.5*self.erg*float(eta@self.Zinv@eta),
                 P0_erg_s=AB, chemical_relative_P0_erg_s=AO-excess,
                 Eescape_star_MeV=self.mu+float(self.t@eta),
                 Lescape_star_erg_s=self.erg*D*(self.mu+float(self.t@eta)),
                 a_signed_keV=1000*float(self.t@eta),
                 q_equals_net_keV=(LH-excess)/(self.erg*D)*1000 if D else math.nan,
                 q_equals_work_keV=abs(float(self.t@eta))*1000,
                 stable=stable, stability='stable' if stable else 'UNSTABLE',
                 eigen_precision=eigen_precision,
                 eig_max_real_s=max(eig.real), tau_slow_yr=tau,
                 build_e_yr=build[0], build_mu_yr=build[1], build_time_yr=energy_build,
                 depletion_build=gamma*energy_build, depletion_relax=gamma*tau,
                 tdrive_yr=tdrive, relax_over_drive=tau/tdrive,
                 reach_time_estimate_yr=max(energy_build, 3*tau),
                 depletion_reach_estimate=gamma*max(energy_build, 3*tau),
                 composition_eta_e=displacement[0], composition_eta_mu=displacement[1],
                 composition_B_e=self.B*self.t[0]/self.N[0]*gamma*energy_build,
                 composition_B_mu=self.B*self.t[1]/self.N[1]*gamma*energy_build,
                 chemical_residual=float(max(abs(R-target))/chem_scale),
                 thermal_residual=(LH-excess-eq-photon)/max(LH, excess+eq+photon),
                 subsystem_residual=(AB+AO-excess-eq-photon)/max(LH, excess+eq+photon),
                 chemical_energy_residual=(AB+AO-LH)/max(LH, abs(AB)+abs(AO)),
                 temperature_in_certificate=1e6 <= T <= 1e9,
                 reachable_in_certificate=bool(1e6 <= T <= 1e9 and gamma*max(energy_build, 3*tau) <= 1e-6),
                 smooth_finiteT_scale_keV=(np.pi**2/6)*(self.kb*T)**2/100*1000,
                 shell_finiteT_scale_keV=self.kb*T*1000)
        return d
