"""Frozen predeclaration acceptance; candidate evidence only, never a baseline."""
import hashlib,json,re,subprocess
from pathlib import Path
import numpy as np
YEAR=365.25*86400
ENTRY='f7116c1408c06f976527f86d4397ad6d4540dedf'
CANONICAL='d019ae390be4f5e3daba05039903485cb497e397'

def check(ok,message):
    if not ok:raise RuntimeError(message)
def table(path):return np.atleast_1d(np.genfromtxt(path,names=True))
def plain(row):return {n:float(row[n]) for n in row.dtype.names}
def crossings(x,root):
    result=[]
    for channel in ['xi_e','xi_mu']:
        exact=np.flatnonzero(x[channel]==root)
        for i in exact:result.append({'channel':channel,'time_bracket_yr':[float(x['t_s'][i]/YEAR)]*2,'estimated_time_yr_log_interpolation':float(x['t_s'][i]/YEAR),'xi_bracket':[root,root],'root':root,'exact_checkpoint':True})
        for i in np.flatnonzero((x[channel][:-1]-root)*(x[channel][1:]-root)<0):
            fraction=float((root-x[channel][i])/(x[channel][i+1]-x[channel][i]))
            estimate=float(np.exp(np.log(x['t_s'][i])*(1-fraction)+np.log(x['t_s'][i+1])*fraction)/YEAR) if x['t_s'][i]>0 else float(x['t_s'][i+1]*fraction/YEAR)
            result.append({'estimated_time_yr_log_interpolation':estimate,'channel':channel,'time_bracket_yr':[float(x['t_s'][i]/YEAR),float(x['t_s'][i+1]/YEAR)],'xi_bracket':[float(x[channel][i]),float(x[channel][i+1])],'root':root,'exact_checkpoint':False})
    return result

def validate(out):
    base=table(out/'trajectory.tsv');ref=table(out/'refined.tsv')
    check(len(base)==len(ref) and np.array_equal(base['t_s'],ref['t_s']),'ODE checkpoint grids differ')
    check(base['t_s'][0]==0 and base['t_s'][-1]==1e10*YEAR,'run interval changed')
    expected_times=np.concatenate(([0],YEAR*10**(np.arange(401)/40)))
    check(len(base)==402 and np.allclose(base['t_s'][:402],expected_times,rtol=1e-12,atol=0),'output cadence changed')
    ledger_residuals={}
    for label,data in [('baseline',base),('refined',ref)]:
        check(data['Tinf_K'][0]==1e8 and data['eta_e_MeV'][0]==0 and data['eta_mu_MeV'][0]==0,'initial physical state changed: '+label)
        for name in data.dtype.names:check(np.all(np.isfinite(data[name])),'nonfinite '+name)
        check(np.all(data['Tinf_K']>0) and np.all(data['eta_e_MeV'][1:]>0) and np.all(data['eta_mu_MeV'][1:]>0),'temperature/imbalance sign: '+label)
        check(np.all(data['R_e']*data['eta_e_MeV']>=0) and np.all(data['R_mu']*data['eta_mu_MeV']>=0),'rate sign: '+label)
        check(np.all(data['LH']>=0) and np.all(data['DeltaLnu']>=0),'power sign: '+label)
        expected={'Pnet':data['LH']-data['DeltaLnu']-data['Lnu_eq']-data['Lgamma']-data['Lother_neutrino'],'DeltaPbeta':data['LH']-data['DeltaLnu'],'Lnu_full':data['Lnu_eq']+data['DeltaLnu'],'x_dot':data['Pnet']/(data['Tinf_K']*data['Cstar'])}
        scale=np.maximum(data['LH']+data['DeltaLnu']+data['Lnu_eq']+data['Lgamma']+data['Lother_neutrino'],1e-300)
        ledger_residuals[label]={}
        for name,expected_values in expected.items():
            denominator=scale/(data['Tinf_K']*data['Cstar']) if name=='x_dot' else scale
            residual=float(np.max(np.abs(data[name]-expected_values)/denominator))
            check(residual<=5e-14,'checkpoint ledger roundoff identity failed: '+name)
            ledger_residuals[label][name]=residual
    check(np.all(base['Tinf_K']>0),'nonpositive temperature')
    check(np.all(base['eta_e_MeV'][1:]>0)&np.all(base['eta_mu_MeV'][1:]>0),'spin-down eta sign')
    check(np.all(base['LH']>=0)&np.all(base['DeltaLnu']>=0),'nonnegative luminosities')
    check(np.all(base['R_e']*base['eta_e_MeV']>=0)&np.all(base['R_mu']*base['eta_mu_MeV']>=0),'reaction sign inconsistency')
    check(base['x_dot'][0]<0 and np.max(base['LH'])>0 and np.max(base['DeltaLnu'])>0,'missing cooling/activation')
    check(np.min(base['DeltaPbeta'])<0<np.max(base['DeltaPbeta']),'incremental sign history missing')
    comparison={}
    for name in base.dtype.names:
        if name=='t_s':continue
        floor=1e-10 if name.startswith('eta_') and name.endswith('MeV') else 1e-300
        comparison[name]=float(np.max(np.abs(base[name]-ref[name])/np.maximum(np.abs(ref[name]),floor)))
    for name in ['Tinf_K','eta_e_MeV','eta_mu_MeV']:check(comparison[name]<=2e-4,'predeclared ODE convergence failed: '+name)
    initial={}
    expected_variants={'initial-T1e7.tsv':(1e7,0),'initial-T1e9.tsv':(1e9,0),'initial-xi1.tsv':(1e8,1),'initial-xi20.tsv':(1e8,20)}
    check({p.name for p in out.glob('initial-*.tsv')}==set(expected_variants),'incorrect IC variants')
    for path in sorted(out.glob('initial-*.tsv')):
        x=table(path)
        for name in x.dtype.names:check(np.all(np.isfinite(x[name])),'nonfinite IC '+name)
        T0,xi0=expected_variants[path.name]
        check(abs(x['Tinf_K'][0]/T0-1)<1e-14 and x['t_s'][-1]==base['t_s'][-1],'IC initial temperature or endpoint changed')
        check(abs(x['xi_e'][0]-xi0)<1e-12 and abs(x['xi_mu'][0]-xi0)<1e-12,'IC chemical input changed')
        errors={n:abs(float(x[n][-1]/base[n][-1]-1)) for n in ['Tinf_K','eta_e_MeV','eta_mu_MeV']}
        check(max(errors.values())<=.01,'initial-condition endpoint convergence failed: '+path.name);initial[path.name]=errors
    check(len(initial)==4,'missing initial-condition validation runs')
    console=(out/'console.log').read_text();coeff={};vectors={}
    for line in console.splitlines():
        if line.startswith('LTILDE '):
            _,process,value,refined,_,spread=line.split();coeff[process]={'value':float(value),'refined':float(refined),'relative':float(spread)}
        if line.startswith('RESULT ') and line.split()[1] in ['W','I']:
            _,name,count,*data=line.split();vectors[name]=list(map(float,data))
    check(set(coeff)=={'2','3'} and set(vectors)=={'W','I'},'missing semantic coefficient evidence')
    eligible=(base['t_s']>=1e9*YEAR*(1-1e-12))&(base['t_s']<=1e10*YEAR)&(np.abs(base['xi_e'])>=100)&(np.abs(base['xi_mu'])>=100)
    qs={'window_yr':[1e9,1e10],'minimum_abs_xi_both':100,'eligible_count':int(sum(eligible))}
    if sum(eligible)==0:
        qs['status']='NOT REACHED';qs['authority']='Phase-5D-1C2 section34 permits NOT REACHED; frozen fit window retained'
    else:
        check(sum(eligible)>=10,'asymptotic regime entered but fewer than10 eligible samples')
        selected=base[eligible];drive=np.abs(selected['Omega']*selected['Omega_dot']);kb_erg=1.380649e-16;kb_mev=kb_erg/1.602176634e-6;ch=24/(11513*np.pi**8)
        qs['status']='PASS';qs['channels']={}
        for i,name in enumerate(['eta_e_MeV','eta_mu_MeV']):
            expected=(2*kb_erg*kb_mev**7*np.abs(vectors['I'][i])*drive/(ch*coeff[str(i+2)]['value']))**(1/7)
            discrepancy=float(np.max(np.abs(selected[name]/expected-1)));slope=float(np.polyfit(np.log(drive),np.log(selected[name]),1)[0]);check(discrepancy<=.05 and abs(slope-1/7)<=.015,'predeclared quasi-steady criterion failed: '+name)
            qs['channels'][name]={'max_relative_asymptote_error':discrepancy,'log_slope':slope,'slope_error':abs(slope-1/7)}
        qs['max_relative_thermal_residual']=float(np.max(np.abs(selected['x_dot']*selected['Tinf_K']*selected['Cstar'])/selected['LH']))
        qs['max_relative_photon_vs_5over8_heating']=float(np.max(np.abs(selected['Lgamma']/(.625*selected['LH'])-1)))
    jac=table(out/'trajectory.tsv.jacobian');jacobian=[]
    for row in jac:
        matrix=np.array([row[n] for n in row.dtype.names[1:]]).reshape(3,3);ev=np.linalg.eigvals(matrix);mod=np.abs(ev)
        jacobian.append({'time_yr':float(row['t_s']/YEAR),'matrix':matrix.tolist(),'eigenvalues_real_s':ev.real.tolist(),'eigenvalues_imag_s':ev.imag.tolist(),'eigenvalue_magnitude_ratio':float(max(mod)/min(mod)),'modal_magnitude_timescales_s':(1/mod).tolist(),'signed_real_part_timescales_s':(1/ev.real).tolist()})
    check(len(jacobian)==3 and np.allclose([j['time_yr'] for j in jacobian],[1e6,1e8,1e10],rtol=1e-10),'missing preselected stiffness epochs')
    metrics={}
    for name in ['trajectory','refined']:
        lines=(out/(name+'.tsv.steps')).read_text().splitlines();v=list(map(float,lines[1].split()));late=np.atleast_2d(np.loadtxt(lines[3:]));metrics[name]=dict(zip(lines[0].split(),v));metrics[name]['late_outputs']=late[-10:].tolist();metrics[name]['rejection_fraction']=float(v[1]/(v[0]+v[1]));metrics[name]['RHS_definition']='all derivative evaluations, including accepted-endpoint validation'
        increments=np.diff(np.r_[0,late[:,1]])
        check(np.all(increments>=0) and np.all(increments==increments.astype(int)),'invalid cumulative accepted-step table')
        metrics[name]['maximum_accepted_steps_per_checkpoint']=int(np.max(increments))
        metrics[name]['last10_accepted_steps_per_checkpoint']=increments[-10:].astype(int).tolist()
    result={'classification':'CONTROLLED MATHEMATICAL / ARCHITECTURE BENCHMARK; CANDIDATE ONLY; NOT GOVERNED BASELINE','canonical':CANONICAL,'frozen_context_plan_sha':ENTRY,'candidate_only':True,'governed_baseline':False,'branch':'physics/phase5d-controlled-rotochemical-evolution','fixture':{'radial_resolution':80000,'EOS_resolution':8192,'rho_c_g_cm3':1.10e15},'normalizations':{'SMe':1e-51,'SMmu':2e-51,'units':'erg cm^-3 s^-1 K^-8','enabled':['Me','Mmu'],'disabled':['De','Dmu']},'initial':{'Tinf_K':1e8,'eta_MeV':[0,0]},'spin':{'B_G':1e8,'P0_s':.001,'PPdot':(1e8/3.2e19)**2},'solver':{'name':'GSL RKF45','rtol':1e-7,'atol':[1e-12,1e-18,1e-18],'refined_rtol':1e-9,'refined_atol':[1e-14,1e-20,1e-20],'control_mapping':{'eps_abs':1,'eps_rel':'rtol','a_y':1,'a_dydt':0,'scale_abs':'component_atol','error_level':'atol_i + rtol*abs(y_i)'}},'Ltilde':coeff,'semantic_W_I':vectors,'checkpoints':[plain(row) for row in base],'final':plain(base[-1]),'ranges':{n:[float(min(base[n])),float(max(base[n]))] for n in base.dtype.names},'ODE_convergence':comparison,'initial_condition_convergence':initial,'quasi_steady':qs,'incremental_root_crossings':crossings(base,4.909710028924132),'full_root_crossings':crossings(base,5.633717467648343),'stiffness':{'step_statistics':metrics,'jacobians':jacobian},'protected_baseline_installed':False}
    root=Path(__file__).resolve().parents[2]
    result['source_hashes']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for directory in ['CompactStar/Physics/Rotochemical','CompactStar/Physics/Evolution','CompactStar/Physics/Driver/Thermal','CompactStar/Analysis','tests/rotochemical'] for p in (root/directory).rglob('*') if p.is_file() and '__pycache__' not in str(p)}
    result['protected_gate']=json.loads((out/'pretrajectory-gate.json').read_text())
    result['controlled_ledger']='Pnet = LH - DeltaLnu - Lnu_eq_controlled - Lgamma - Lnu_other; one From conversion'
    result['run_provenance']=json.loads((out/'run-provenance.json').read_text())
    result['ledger_residuals']=ledger_residuals
    result['physical_Ltilde_radial_comparison']=[{'process_index':int(line.split()[1]),'radial20000_value':float(line.split()[2]),'relative_difference':float(line.split()[4]),'acceptance_goal':5e-3} for line in console.splitlines() if line.startswith('RADIAL_LTILDE ')]
    check(len(result['physical_Ltilde_radial_comparison'])==2,'missing physical radial Ltilde comparison')
    result['thermal_source_hashes']={p:h for p,h in result['run_provenance']['inputs'].items() if '/thermal/eos.' in p}
    result['normalization_classification']='MATHEMATICAL / ARCHITECTURE BENCHMARK; no FR2005 absolute normalization claim'
    result['stiffness']['conclusion']='RKF45 completed this controlled benchmark within the predeclared step budget; no DU or superfluid adequacy claim'
    result['semantic_Z']=[[float(v) for v in line.split()[3:5]] for line in console.splitlines() if line.startswith('FROZEN_ROW ')]
    Z=np.asarray(result['semantic_Z']);timescales=[]
    for age in [1e6,1e8,1e10]:
        row=base[np.argmin(abs(base['t_s']/YEAR-age))];reaction=Z@np.array([row['R_e'],row['R_mu']])
        def ratio(a,b):return float(abs(a/b)) if b!=0 else None
        timescales.append({'time_yr':float(row['t_s']/YEAR),'thermal_T_over_abs_Tdot_s':ratio(1,row['x_dot']),'spin_Omega_over_abs_OmegaDot_s':ratio(row['Omega'],row['Omega_dot']),'reaction_eta_over_abs_ZR_s':[ratio(row['eta_e_MeV'],reaction[0]),ratio(row['eta_mu_MeV'],reaction[1])],'definition':'instantaneous magnitude scales; reaction excludes compensating spin drive; null denotes zero denominator'})
    result['stiffness']['physical_timescales']=timescales
    result['power_convergence_absolute']={n:float(np.max(np.abs(base[n]-ref[n]))) for n in ['LH','DeltaLnu','DeltaPbeta']}
    aggregate=[]
    for i in np.flatnonzero(base['DeltaPbeta'][:-1]*base['DeltaPbeta'][1:]<0):
        fraction=float(-base['DeltaPbeta'][i]/(base['DeltaPbeta'][i+1]-base['DeltaPbeta'][i]))
        aggregate.append({'time_bracket_yr':[float(base['t_s'][i]/YEAR),float(base['t_s'][i+1]/YEAR)],'power_bracket_erg_s':[float(base['DeltaPbeta'][i]),float(base['DeltaPbeta'][i+1])],'estimated_time_yr_log_interpolation':float(np.exp((1-fraction)*np.log(base['t_s'][i])+fraction*np.log(base['t_s'][i+1]))/YEAR),'xi_e_bracket':[float(base['xi_e'][i]),float(base['xi_e'][i+1])],'xi_mu_bracket':[float(base['xi_mu'][i]),float(base['xi_mu'][i+1])]})
    result['aggregate_incremental_power_crossings']=aggregate
    relative=np.abs(base['Pnet']-ref['Pnet'])/np.maximum(np.abs(ref['Pnet']),1e-300);index=int(np.argmax(relative));absolute=float(abs(base['Pnet'][index]-ref['Pnet'][index]));gross=float(base['LH'][index]+base['DeltaLnu'][index]+base['Lnu_eq'][index]+base['Lgamma'][index]+base['Lother_neutrino'][index])
    result['net_power_cancellation_diagnostic']={'time_yr':float(base['t_s'][index]/YEAR),'maximum_relative_residual_difference':float(relative[index]),'absolute_difference_at_that_checkpoint_erg_s':absolute,'difference_over_gross_ledger_power':absolute/gross,'maximum_absolute_Pnet_difference_all_checkpoints_erg_s':float(np.max(np.abs(base['Pnet']-ref['Pnet'])))}
    print('CONTROLLED TRAJECTORY VALIDATION COMPLETED',json.dumps({'final':result['final'],'ODE':comparison,'quasi_steady':qs}),flush=True)
    return result

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);args=parser.parse_args()
    result=validate(args.directory)
    (args.directory/'candidate.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
