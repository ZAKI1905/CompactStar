#!/usr/bin/env python3
"""Measured matched-reference comparison; never promotes it to qualification."""
import csv,json,pathlib,sys,math
from decimal import Decimal as D,localcontext
from run import ROOT,BASE,digest,save
sys.path.insert(0,str(ROOT/'tools/linux/native_qualification'))
from verify_chemical_reference import rows,lift
SOURCE=BASE/'stage-a-parallel-20261110-r1/source-stage/run/trajectory/source'

def analyze(out,destination):
 destination.mkdir(parents=True,exist_ok=True)
 inputs=[SOURCE/'reference-4-observations.tsv',out/'control/reference-4-observations.tsv']
 source,control=map(rows,inputs);assert len(source)==len(control)==1921
 assert [r['t'] for r in source]==[r['t'] for r in control]
 channels=['P_dir_actual','LH','Lnu_full','Lgamma','Pnet','Lout','DeltaLnu']
 with localcontext() as ctx:
  ctx.prec=110;energies={};cumulative={};convergence={}
  for role,data in [('source',source),('control',control)]:
   cumulative[role]={k:[D(0)] for k in channels}
   for a,b in zip(data,data[1:]):
    dt=lift(b['t'])-lift(a['t'])
    for k in channels:cumulative[role][k].append(cumulative[role][k][-1]+dt*(lift(a[k])+lift(b[k]))/2)
   energies[role]={k:v[-1] for k,v in cumulative[role].items()}
   energies[role]['Echem_final_erg']=lift(data[-1]['Echem'])*lift(1.602176634e33/1e39)
   values=[rows((SOURCE if role=='source' else out/'control')/f'reference-{i}-observations.tsv') for i in range(5)]
   convergence[role]={}
   for k in ('x','e','mu','Tinf'):
    maxima=[max(abs(lift(a[k])-lift(b[k])) for a,b in zip(values[i],values[j])) for i,j in ((0,1),(1,2),(3,4),(2,4))]
    convergence[role][k]={'maximum_pairwise_differences':maxima,'twice_maximum_measured_difference':2*max(maxima),'is_complete_uncertainty':False}
  summary=dict(status='PRELIMINARY — NOT NUMERICALLY QUALIFIED',source_Tinf_K=lift(source[-1]['Tinf']),control_Tinf_K=lift(control[-1]['Tinf']),source_minus_control_K=lift(source[-1]['Tinf'])-lift(control[-1]['Tinf']),source_rise_K=lift(source[-1]['Tinf'])-D('1e8'),control_change_K=lift(control[-1]['Tinf'])-D('1e8'),energies_erg=energies,reference_convergence=convergence,final_chemical_imbalance_MeV={role:[data[-1]['e'],data[-1]['mu']] for role,data in [('source',source),('control',control)]},exact_zero_control_chemistry=all(float(r[k])==0 for r in control for k in ('e','mu','LH','Echem','P_dir_actual','P_dir_eq','DeltaLnu')),input_sha256={str(p):digest(p) for p in inputs},physical_scope='Controlled P2 mathematical neutron sink on a frozen certified-input candidate background; not a realistic decay-rate prediction.',uncertainty_status='Measured convergence and conditional high-precision estimates are separate from unresolved formal uncertainty bounds.')
  supplemental=out/'supplement.json'
  if supplemental.exists():
   q=json.loads(supplemental.read_text());assert not q['partial_matrix']
   summary['conditional_numerical_temperature_uncertainty_K']={k:v for k,v in q['paired'].items() if 'temperature_uncertainty' in k}
   summary['input_sha256'][str(supplemental)]=digest(supplemental)
   summary['uncertainty_status']='Conditional analytic enclosure with fixed accepted inputs; independent review of the derivation and finite-precision assumptions remains required.'
  save(destination/'comparison.json',summary)
  with (destination/'temperature.csv').open('w') as f:
   writer=csv.writer(f);writer.writerow(['time_s','time_Julian_year','source_Tinf_K','control_Tinf_K','source_minus_control_K'])
   for a,b in zip(source,control):writer.writerow([a['t'],lift(a['t'])/31557600,a['Tinf'],b['Tinf'],lift(a['Tinf'])-lift(b['Tinf'])])
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 import numpy as np
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160})
 t=np.array([float(r['t'])/31557600 for r in source]);Ts=np.array([float(r['Tinf']) for r in source]);Tc=np.array([float(r['Tinf']) for r in control])
 label='PRELIMINARY — NOT NUMERICALLY QUALIFIED'
 fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
 axes[0].plot(t,Ts/1e6,label='P2 source',color='#b45309');axes[0].plot(t,Tc/1e6,label='Zero-BNV control',color='#0369a1')
 axes[0].set(xlabel='Time [Julian years]',ylabel=r'Redshifted temperature $T_\infty$ [$10^6$ K]');axes[0].ticklabel_format(useOffset=False);axes[0].legend()
 axes[1].plot(t,Ts-Tc,color='#7c3aed');axes[1].set(xlabel='Time [Julian years]',ylabel=r'$T_{\infty,\mathrm{source}}-T_{\infty,\mathrm{control}}$ [K]')
 axes[1].annotate(f'{Ts[-1]-Tc[-1]:,.6f} K',(.98,Ts[-1]-Tc[-1]),xytext=(-8,-22),textcoords='offset points',ha='right',color='#6d28d9')
 for ax in axes:ax.grid(alpha=.18)
 fig.suptitle(label,color='#991b1b',fontsize=12,fontweight='bold')
 if 'conditional_numerical_temperature_uncertainty_K' in summary:
  u=float(summary['conditional_numerical_temperature_uncertainty_K']['plotted_reference_maximum_temperature_uncertainty_K'])
  unit=10.0**(math.floor(math.log10(u))-1);display_bound=math.ceil(u/unit)*unit
  fig.supxlabel(f'Conditional numerical uncertainty in the difference ≤ {display_bound:.2g} K (below plot resolution); assumptions in report.',fontsize=8)
 for extension in ('png','pdf'):fig.savefig(destination/f'temperature-comparison.{extension}')
 plt.close(fig)
 fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
 for k,name,color in [('P_dir_actual','Direct BNV deposition','#b45309'),('Lgamma','Photon cooling','#0369a1'),('Lnu_full','Neutrino cooling','#15803d')]:
  axes[0,0].plot(t,[float(r[k]) for r in source],label=name,color=color)
  if k!='P_dir_actual':axes[0,0].plot(t,[float(r[k]) for r in control],ls='--',color=color)
 axes[0,0].set(yscale='log',ylabel='Power [erg/s]',title='Solid: source; dashed: control');axes[0,0].legend(fontsize=8)
 axes[0,1].plot(t,[float(r['LH']) for r in source],color='#be123c');axes[0,1].set(ylabel='Chemical-reaction heating [erg/s]',title='Control chemical heating is exactly zero')
 for k,name in [('e','Electron'),('mu','Muon')]:axes[1,0].plot(t,[float(r[k])*1e10 for r in source],label=name)
 axes[1,0].axhline(0,color='#0369a1',ls='--',label='Control');axes[1,0].set(ylabel=r'Chemical imbalance [$10^{-10}$ MeV]');axes[1,0].legend()
 for k,name in [('P_dir_actual','Direct deposition'),('Lgamma','Photon losses'),('Lnu_full','Neutrino losses'),('Pnet','Net thermal input')]:axes[1,1].plot(t,[float(v)/1e42 for v in cumulative['source'][k]],label=name)
 axes[1,1].set(ylabel=r'Source cumulative energy [$10^{42}$ erg]');axes[1,1].legend(fontsize=8)
 for ax in axes.flat:ax.set_xlabel('Time [Julian years]');ax.grid(alpha=.18)
 fig.suptitle(label,color='#991b1b',fontsize=12,fontweight='bold')
 for extension in ('png','pdf'):fig.savefig(destination/f'physics-ledger.{extension}')
 plt.close(fig)
 print(json.dumps({k:str(summary[k]) for k in ('source_Tinf_K','control_Tinf_K','source_minus_control_K')},indent=2))

if __name__=='__main__':analyze(pathlib.Path(sys.argv[1]),pathlib.Path(sys.argv[2]))
