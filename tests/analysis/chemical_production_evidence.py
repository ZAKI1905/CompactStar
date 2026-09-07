#!/usr/bin/env python3
"""Test-only certificate assembly and production fixture runner.

Expected EOS/refinement/tail mathematics remain independent of production.
No baseline file is read. Goal/envelope authority is the immutable predeclaration.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import platform
import struct
import copy
import sys
sys.dont_write_bytecode = True
import numpy as np
import mpmath as mp
from chemical_trackr_budget import Model, matrix, weight, tail, integrate


def tail_certificate(d, model, measured):
    r, m, nu, nb = map(float, d[-1, :4])
    with mp.workdps(80):
        R, M, N = map(mp.mpf, (r, m, nb))
        p, e = map(mp.mpf, model.mass[1:3]); hc = mp.mpf(model.hc)
        pf = hc * (3*mp.pi**2*N)**(mp.mpf(1)/3)
        excess = sum(pf*pf/(mp.sqrt(ma*ma+pf*pf)+ma) for ma in (p,e))
        h = mp.log1p(excess/(p+e))
        hu = float(np.nextafter(float(h), np.inf))
        # Positive pressure and m(r)>=M imply H >= M(1/Rcut-1/R0).
        # This first upper radius is independent of the unknown total mass.
        bootstrap = R/(1-mp.mpf(hu)*R/M)
        rb = float(np.nextafter(float(bootstrap), np.inf))
        _, eps, _, _ = model.pe_at_h(hu)
        eps = float(np.nextafter(eps*(1+256*np.finfo(float).eps), np.inf))
        mu = M+4*mp.pi/3*mp.mpf(eps)*(mp.mpf(rb)**3-R**3)
        mu = float(np.nextafter(float(mu), np.inf))
        # d[nu - log(1-2m/r)/2]/dr = 4*pi*r²*(eps+P)/(r-2m)>=0.
        # Thus R0 <= 2*M_total/[1-(1-2*M_cut/Rcut)*exp(2*Hcut)].
        # M_total is replaced by its UPPER enclosure, not M_cut. M_cut in
        # the denominator is the known boundary mass, a different quantity.
        denominator = 1-(1-2*M/R)*mp.exp(2*mp.mpf(hu))
        assert denominator > 0
        ru = float(np.nextafter(float(2*mp.mpf(mu)/denominator), np.inf))
        assert r < ru <= rb and r > 2*mu
        chip = mp.sqrt(p*p+pf*pf)*pf/(mp.pi**2*hc**3)
        chie = mp.sqrt(e*e+pf*pf)*pf/(mp.pi**2*hc**3)
        cp = mp.mpf(float(np.nextafter(float(chip*chie/(chip+chie)),np.inf)))
        lapse_error = mp.mpf(float(np.nextafter(float((mp.mpf(mu)-M)/(R*(1-2*mp.mpf(mu)/R))),np.inf)))
        volume = 4*mp.pi/3*(mp.mpf(ru)**3-R**3)
        bound = volume*mp.exp(-mp.mpf(nu)+lapse_error)*cp/mp.sqrt(1-2*mp.mpf(mu)/R)*10**54
        bound = float(np.nextafter(float(bound)*(1+256*np.finfo(float).eps), np.inf))
    assert measured['direct_R_vacuum_km'] <= ru
    assert measured['direct_mass_increment_km'] <= mu-m
    assert measured['direct_Gee_count_per_MeV']+measured['direct_numerical_error'] < bound
    assert measured['independent_60_digit_exterior_shell_integral']+measured['exterior_self_gravity_allowance'] < bound
    result = {'classification':'analytic comparison certificate under positive monotone pe source hypotheses',
        'R_cut_km':r,'M_cut_km':m,'nu_cut':nu,'susceptibility_upper':float(cp),'h_upper':hu,'bootstrap_radius_upper_km':rb,
        'M_total_upper_km':mu,'R_upper_km':ru,'epsilon_upper_km_minus2':eps,
        'bound_Gee':bound,'lapse_normalization_error':float(lapse_error),
        'radius_inequality':'R0 <= 2 M_total_upper / [1-(1-2 M_cut/R_cut) exp(2 h_upper)]',
        'mass_inequality':'M_total <= M_cut + (4 pi/3) epsilon_upper (R_bootstrap^3-R_cut^3)',
        'source_hypotheses':'pe; positive pressure/energy; monotone energy; authenticated cut; no horizon',
        'direct_same_cut':measured['direct_Gee_count_per_MeV'],
        'independent_high_precision_shell':measured['independent_60_digit_exterior_shell_integral'],
        'direct_numerical_error':measured['direct_numerical_error'],
        'direct_and_independent_enclosed':True}
    error=np.zeros((3,3));error[1,1]=bound
    return result,error


def certificates(directory, characterized, goals):
    d=np.loadtxt(directory/'profile.tsv',skiprows=1);model=Model(directory)
    r,m,nu,nb=d[:,:4].T
    onsets=np.interp(model.onsets,nb[::-1],r[::-1])
    cuts=np.unique(np.r_[0,r,onsets])
    intervals=[(a,b,int(a in onsets),int(b in onsets)) for a,b in zip(cuts[:-1],cuts[1:])]
    raw=np.atleast_2d(np.loadtxt(directory/'windows.tsv',skiprows=1))
    neutron_upper,nn_guard,root_error=model.neutron_source_upper()
    mu_lower,mu_upper,mu_margin=model.muon_source_interval()
    # Use a source-derived ULP neighbourhood on the otherwise available pe
    # side as well: interpolation may round a mapped density to exact onset.
    bands=[(float(np.nextafter(model.onsets[0]-64*abs(np.spacing(model.onsets[0])),-np.inf)),neutron_upper,'neutron'),(mu_lower,mu_upper,'muon')]
    entries=[]
    for lo,hi,name in bands:
        measured=raw[raw[:,0]==model.onsets[0 if name=='neutron' else 1]]
        assert np.all(measured[:,2:4]>=lo) and np.all(measured[:,2:4]<=hi)
        left,right=np.interp([hi,lo],nb[::-1],r[::-1])
        # Outward representational radius margin, 64 radius ULPs independent
        # of observed provider answers. This covers interpolation/node rounding.
        left=float(np.nextafter(left-64*abs(np.spacing(left)),-np.inf))
        right=float(np.nextafter(right+64*abs(np.spacing(right)),np.inf))
        first=int(np.searchsorted(cuts,left,side='right')-1)
        last=int(np.searchsorted(cuts,right,side='left')-1)
        assert 0<=first<=last<len(cuts)-1
        containing_left,containing_right=cuts[first],cuts[last+1]
        assert containing_left<=left<right<=containing_right
        geometry_nodes=np.unique(np.r_[containing_left,containing_right,r[(r>=containing_left)&(r<=containing_right)]])
        mv=np.interp(geometry_nodes,r,m);nv=np.interp(geometry_nodes,r,nu)
        radius_upper=float(np.nextafter(containing_right,np.inf));mass_upper=float(np.nextafter(max(mv),np.inf));nu_lower=float(np.nextafter(min(nv),-np.inf));mass_lower=float(np.nextafter(min(mv),0))
        assert left>2*mass_upper and mass_lower>0
        lower_n=float(np.nextafter(np.interp(right,r,nb),-np.inf));upper_n=float(np.nextafter(np.interp(left,r,nb),np.inf))
        hlo,_=model.mp_state(lower_n);hhi,chi=model.mp_state(upper_n)
        with mp.workdps(70):dh=float(np.nextafter(float(hhi-hlo),np.inf))
        chi=np.array([0. if x==0 else np.nextafter(float(x),np.inf) for x in chi])
        cmax=np.diag(chi[[0,2,3]]);cmax[1,2]=cmax[2,1]=min(chi[2],chi[3])
        physical_width=dh*radius_upper*(radius_upper-2*mass_lower)/mass_lower
        width=max(right-left,physical_width)
        maximum_weight=4*np.pi*radius_upper**2*np.exp(-nu_lower)/np.sqrt(1-2*mass_upper/left)
        error=np.nextafter(cmax*width*maximum_weight*1e54*(1+256*np.finfo(float).eps),np.inf)
        error[cmax==0]=0
        entries.append({'name':name,'radial_left_branch':'npe' if name=='neutron' else 'npemu','radial_right_branch':'pe' if name=='neutron' else 'npe','original_profile_first_cell':int(np.searchsorted(r,left,side='right')-1),'original_profile_last_cell':int(np.searchsorted(r,right,side='left')-1),'left':left,'right':right,'containing_left':float(containing_left),'containing_right':float(containing_right),'first_cell':first,'last_cell':last,'radius_upper':radius_upper,'mass_upper':mass_upper,'nu_lower':nu_lower,'source_density_interval':[lo,hi],'actual_expanded_density_interval':[lower_n,upper_n],'provider_brackets':measured.tolist(),'whole_interval_C_upper':cmax.tolist(),'physical_width_bound':physical_width,'error':error.tolist(),'availability':'source inequalities, physical branches retained; finite refusal interval excluded','source_guard':{'neutron_ulps':2**30,'nn_guard':nn_guard,'root_density_error_upper':root_error} if name=='neutron' else {'potential_margin':mu_margin}})
    tail,et=tail_certificate(d,model,characterized['tail'])
    errors=characterized['error_G_components']
    bg=sum((np.array(errors[k]) for k in ['background_characterization','provider_nodal_comparison','equilibrium_anchor']),np.zeros((3,3)))
    # Account for normalization of the added positive tail mass separately.
    bg+=abs(np.array(characterized['G']))*tail['lapse_normalization_error']
    center=np.array(errors['centre'])
    return {'intervals':intervals,'refusals':entries,'tail':tail,'center_error':center.tolist(),'background_error':bg.tolist(),'tail_error':et.tolist(),'goals':goals},d


def transport(cert, path, fixed):
    lines=[str(len(cert['intervals']))]
    lines.extend(' '.join(map(str,row)) for row in cert['intervals'])
    lines.append(str(len(cert['refusals'])))
    for f in cert['refusals']:
        lines.append(' '.join(str(f[k]) for k in ['left','right','containing_left','containing_right','first_cell','last_cell','radius_upper','mass_upper','nu_lower']))
        lines.extend(' '.join(map(str,row)) for row in f['error'])
    tail=cert['tail'];lines.append(' '.join(str(tail[k]) for k in ['R_cut_km','M_cut_km','nu_cut','h_upper','epsilon_upper_km_minus2','bootstrap_radius_upper_km','M_total_upper_km','R_upper_km','susceptibility_upper','lapse_normalization_error','bound_Gee']))
    for m in [cert['center_error'],cert['tail_error'],cert['background_error'],fixed['goals']['G'],fixed['goals']['Q'],fixed['goals']['Z']]:
        lines.extend(' '.join(map(str,row)) for row in m)
    for v in [fixed['V_I_validation'],fixed['goals']['W_numerical'],fixed['goals']['W_validation']]:lines.append(' '.join(map(str,v)))
    path.write_text('\n'.join(lines)+'\n')


def run():
    parser=argparse.ArgumentParser();parser.add_argument('executable',type=Path);parser.add_argument('characterization',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--global-only',action='store_true');parser.add_argument('--fixture',default='t8192-r80000');args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).resolve().parents[2]
    fixed=json.loads((source/'docs/validation/phase5c2_preproduction_evidence.json').read_text())
    char=json.loads(args.characterization.read_text());directory=args.characterization.parent/args.fixture
    if args.fixture!='t8192-r80000':
        dd=np.loadtxt(directory/'profile.tsv',skiprows=1);mm=Model(directory);char=copy.deepcopy(char);char['tail']=tail(dd,mm,matrix(integrate(dd,mm,16)[0]))[0]
    cert,d=certificates(directory,char,fixed['goals']);transport(cert,args.output/'certificate.txt',fixed)
    (args.output/'certificates.json').write_text(json.dumps(cert,indent=2,sort_keys=True,allow_nan=False)+'\n')
    with (args.output/'production.log').open('w') as log:
        result=subprocess.run([str(args.executable.resolve()),str(directory),str(args.output/'certificate.txt'),str(args.output/'fresh-production'), args.fixture.split('-r')[1],'global' if args.global_only else 'full'],stdout=log,stderr=subprocess.STDOUT,cwd=source)
    print('production raw rc',result.returncode,flush=True)
    if result.returncode:raise SystemExit(result.returncode)
    values={};runtime={}
    matrix_names={'G','E_G','E_quadrature','Q','E_Q','E_Schur_arithmetic','Z','E_Z','E_global_solve','E_Z_arithmetic'}
    for line in (args.output/'production.log').read_text().splitlines():
        if line.startswith('PROVENANCE '):
            _,key,value=line.split(' ',2);runtime[key]=value;continue
        if not line.startswith('RESULT '):continue
        _,key,n,*data=line.split();n=int(n);data=list(map(float,data));values[key]=np.array(data).reshape(n,n).tolist() if key in matrix_names or key.startswith(('G_ladder_', 'E_local_', 'E_center', 'E_tail', 'E_background', 'E_refusal')) else data
    assert 'Z' in values
    z=np.array(values['Z']);ez=np.array(values['E_Z']);old=np.array(char['old_Z'])
    # Old F2005 remains independently integrated intrinsic test mathematics.
    old_error=np.array(char['correction_gate']['combined_uncertainties'])-np.array([char['E_Z'][0][0],char['E_Z'][0][1],char['E_Z'][1][1]])
    delta=abs(z-old);combined=np.array([ez[0,0],ez[0,1],ez[1,1]])+old_error
    separations=delta[np.triu_indices(2)];margins=separations/combined
    assert np.all(margins>1)
    independent=np.array(char['independent_background']['G']);eg=np.array(values['E_G'])
    assert np.all(abs(independent-np.array(values['G']))<=eg)
    if not args.global_only:
        w=np.array(values['W']);iw=np.array(values['I']);ew=np.array(values['E_W_numerical']);vw=np.array(values['V_W_validation'])
        assert np.all(ew<=fixed['goals']['W_numerical']) and np.all(vw<=fixed['goals']['W_validation'])
        action=np.array(values['spin_action']);c=299792.458
        mutants={'flip_Omega_dot':-action,'geometric_without_c_minus2':action*c*c,'double_c_minus2':action/(c*c),'omit_e':-4*z@np.array([0,iw[1]]),'omit_mu':-4*z@np.array([iw[0],0]),'swap_channels':-4*z@iw[::-1]}
        assert all(np.linalg.norm(v-action)>4*np.linalg.norm(ew) for v in mutants.values())
        assert np.all(abs(w-z@iw)<=np.array(values['E_W_arithmetic']))
        assert np.all(vw>=abs(z)@np.array(fixed['V_I_validation']))
    partition=values.pop('partition')
    provenance_files=['CompactStar/Analysis/ChemicalResponse.hpp','CompactStar/Analysis/src/ChemicalResponse.cpp','CompactStar/Analysis/ParticleNumberResponse.hpp','CompactStar/Analysis/src/ParticleNumberResponse.cpp','CompactStar/EOS/TrackRFreeGasThermodynamics.hpp','CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp','CompactStar/EOS/LocalThermodynamics.hpp','CompactStar/EOS/src/LocalThermodynamics.cpp','tests/analysis/chemical_production_evidence.py','tests/analysis/chemical_production_fixture.cpp','tests/analysis/chemical_production_contract.cpp','tests/analysis/chemical_production_validation.py','tests/analysis/chemical_exact_oracles.py','tests/analysis/chemical_curved_gc9.py','tests/analysis/chemical_trackr_budget.py','tests/analysis/chemical_trackr_fixture.cpp','docs/validation/phase5c2_preproduction_evidence.json']
    provenance={'particle_constants_and_provider_source':(directory/'model.txt').read_text(),'particle_constants_source_sha256':hashlib.sha256((directory/'model.txt').read_bytes()).hexdigest(),'source_sha256':{p:hashlib.sha256((source/p).read_bytes()).hexdigest() for p in provenance_files},'profile_sha256':hashlib.sha256((directory/'profile.tsv').read_bytes()).hexdigest(),'table_sha256':hashlib.sha256((directory/'freegas.tsv').read_bytes()).hexdigest(),'partition_sha256':hashlib.sha256(b''.join(struct.pack('>d',x) for x in partition)).hexdigest(),'partition':partition,'basis':['Neutron','Electron','Muon'],'channels':['Npe','NpMu'],'model':char['model_identity'],'metric':'nu=Phi; exp(-nu) exactly once; proper volume exactly once','quadrature_orders':[8,16,32,32],'last_ladder_partition':'bisected','accumulation':'Neumaier ordered segment/node sums','default_order':16,'toolchain':{'compiler':runtime['compiler'],'cplusplus':int(runtime['cplusplus']),'configuration':runtime['configuration'],'platform':platform.system(),'architecture':platform.machine(),'GSL':runtime['gsl']},'profile_versions':values['profile_versions'],'structural_zeros':np.array(values['structural_zero_indices'],dtype=int).reshape(-1,2).tolist(),'local_H_semantics':'declared binary neutral Hessian; solve/congruence numerical error plus separately characterized source/provider errors'}
    correction={'named_order':['Z_npe','Z_np','Z_npmu'],'old_Z':old.tolist(),'absolute_separations':separations.tolist(),'combined_uncertainties':combined.tolist(),'margins':margins.tolist(),'minimum_margin':float(min(margins)),'matrix_frobenius_separation':float(np.linalg.norm(z-old)),'matrix_margin':float(np.linalg.norm(z-old)/np.linalg.norm([[combined[0],combined[1]],[combined[1],combined[2]]])),'PASS':True,'A18_substitute':False}
    artifact={'candidate_only':True,'A18_validated':False,'GC13_realistic':False,'INV11_resolved':False,'BNV':False,'predeclaration_sha':'a87f0212c2bd7bfba92db91dfac82447a6561334','fixture':'Structure-1 midpoint, whole star, rho_c=1.10e15 g/cm^3','values':values,'provenance':provenance,'correction_gate':correction,'independent_background_difference':abs(independent-np.array(values['G'])).tolist(),'certificate':cert,'units':{'G':'count/MeV','Q':'count/MeV','Z':'MeV/count','I':'count s^2','W':'MeV s^2'},'numerical_error':{'classification':'propagated numerical uncertainty for the declared computation','G':values['E_G'],'Q':values['E_Q'],'Z':values['E_Z'],'I':values.get('E_I_numerical'),'W':values.get('E_W_numerical')},'validation_envelope':{'classification':'predeclared empirical structural stability envelope','certified_continuum_bound':False,'confidence_interval':False,'probability_statement':False,'I':values.get('V_I_validation'),'W':values.get('V_W_validation')},'goals':fixed['goals']}
    (args.output/'candidate.json').write_text(json.dumps(artifact,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(hashlib.sha256((args.output/'candidate.json').read_bytes()).hexdigest(),args.output/'candidate.json')


if __name__=='__main__':run()
