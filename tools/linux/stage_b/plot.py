#!/usr/bin/env python3
"""Scientific figures from saved curves; no numerical recomputation."""
import csv,json,pathlib,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=pathlib.Path(sys.argv[1]);plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
qualification=json.loads((out/'qualification.json').read_text())
status='Conditional verification candidate' if qualification['status']=='CONDITIONAL_NUMERICAL_VERIFICATION_CANDIDATE' else 'Preliminary: unresolved numerical gates'
def savefig(fig,name):
 fig.savefig(out/(name+'.png'));fig.savefig(out/(name+'.pdf'));plt.close(fig)
def time_axis(a):
 a.set_xlim(0,2048);a.set_xticks([0,512,1024,1536,2048]);a.grid(alpha=.18)
def read(name):
 with (out/name).open() as f:return {k:np.array([float(r[k]) for r in rs]) for rs in [list(csv.DictReader(f))] for k in rs[0]}
p=read('temperature.csv');s=read('curve-source.csv');c=read('curve-control.csv');t=p['year']
fig,ax=plt.subplots(2,1,figsize=(8,6),sharex=True,constrained_layout=True)
ax[0].plot(t,p['T_source']/1e6,label='P2 source',color='#b3422b');ax[0].plot(t,p['T_control']/1e6,label='Zero-source control',color='#236398')
ax[0].set_ylabel(r'$T_\infty$ [million K]')
ax[1].plot(t,p['delta_T']/1e6,color='#68439a');ax[1].set_ylabel(r'$T_{source}-T_{control}$ [million K]');ax[1].set_xlabel('Time [Julian years]')
hp=json.loads((out/'hp-source-1000000.json').read_text())
for a in ax:
 for i,e in enumerate(hp['crossings']):a.axvline(float(e['left_s'])/31557600,color='.5',ls=':',lw=.8,label='Source Cstar knots' if i==0 and a is ax[0] else None)
 time_axis(a)
ax[0].legend()
fig.suptitle('2,048-year mathematical P2 benchmark\n'+status+'; fixed accepted stellar inputs',fontsize=11)
savefig(fig,'temperature')
fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
for key,label,color,style in [('Pdir','Direct BNV','#b3422b','-'),('Lgamma','Photon','#236398','-'),('Lnu','Neutrino','#69623b','-'),('LH','Chemical heating','#68439a','-'),('DeltaLnu','Excess neutrino','#a76622',':')]:
 line,=ax[0,0].semilogy(t,s[key],label=label+' (source)',color=color,ls=style)
 if key in ('Lgamma','Lnu'):ax[0,0].semilogy(t,c[key],ls='--',color=line.get_color(),label=label+' (control)')
ax[0,0].set_ylabel('Power [erg/s]');ax[0,0].legend(fontsize=7,ncol=2,loc='center left',bbox_to_anchor=(0,.57))
ax[0,1].plot(t,s['e']*1e6,label=r'$\eta_e$',color='#236398');ax[0,1].plot(t,s['mu']*1e6,label=r'$\eta_\mu$',color='#68439a',ls='--');ax[0,1].set_ylabel('Chemical imbalance [eV]');ax[0,1].legend()
ax[1,0].plot(t,s['LH']/s['Pdir'],label='Heating / direct BNV',color='#68439a');ax[1,0].plot(t,(s['LH']-s['DeltaLnu'])/s['Pdir'],label='Net beta / direct BNV',color='#a76622',ls='--');ax[1,0].set_ylabel('Power ratio');ax[1,0].legend(fontsize=8)
ax[1,1].semilogy(t,s['Echem'],color='#68439a');ax[1,1].set_ylabel('Stored chemical energy [erg]')
ax[1,1].text(.03,.97,'Control chemistry is exactly zero',transform=ax[1,1].transAxes,va='top',fontsize=8)
for a in ax.flat:a.set_xlabel('Time [Julian years]');time_axis(a)
fig.suptitle('Stage B coupled response — '+status.lower(),fontsize=12)
savefig(fig,'physics')

summary=json.loads((out/'final-summary.json').read_text())
fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
for role,color,marker in [('source','#b3422b','o'),('control','#236398','s')]:
 r=summary['roles'][role]
 values=[r['accuracy_maxima']['known_terms_state'],r['accepted_restart_state_ratio'],r['accuracy_maxima']['power'],r['accepted_power_ratio']]
 ax[0,0].plot(range(4),list(map(float,values)),label=role,color=color,marker=marker,ls='none')
ax[0,0].set_xticks(range(4),['Scientific\nstates','Accepted /\nrestart states','Scientific\npowers','Accepted /\nrestart powers'],fontsize=8)
ax[0,0].set_yscale('log');ax[0,0].axhline(1,color='.3',ls='--');ax[0,0].set_ylabel('Maximum error / allowed budget');ax[0,0].legend(fontsize=8)
configs=summary['matched_configurations'];idx=np.arange(12)
for role,color,marker in [('source','#b3422b','o'),('control','#236398','s')]:
 ax[0,1].plot(idx,[float(r[role+'_R20_over_N'])/2e-4 for r in configs],color=color,marker=marker,label=role)
ax[0,1].axhline(1,color='.3',ls='--',label='Acceptance limit');ax[0,1].set_ylim(0,1.08);ax[0,1].set_ylabel('R20 bound / allowed budget');ax[0,1].legend(fontsize=8)
ax[1,0].semilogy(idx,[float(r['temperature']['total_error_bound_about_native_K']) for r in configs],color='#68439a',marker='o');ax[1,0].set_ylabel('Matched endpoint error bound [K]')
labels=[f'{tier}{cap}' for tier in ('B','R','U') for cap in range(4)]
for a in (ax[0,1],ax[1,0]):a.set_xticks(idx,labels,fontsize=8);a.set_xlabel('Tier B/R/U; cap 0/1/2/3')
for i,key in enumerate(('LH','net_beta_loss')):
 lo,hi=np.array(summary['chemical_temperature_comparison'][key]['temperature_effect_K'],dtype=float)*1e6
 mid=(lo+hi)/2;ax[1,1].errorbar(mid,i,xerr=[[mid-lo],[hi-mid]],fmt='o',color=('#68439a','#a76622')[i],capsize=5)
ax[1,1].axvline(0,color='.4',lw=.8);ax[1,1].set_yticks([0,1],['Heating only','Heating − excess ν']);ax[1,1].set_ylim(-.6,1.6);ax[1,1].set_xlabel('Isolated temperature contribution [μK]')
ax[1,1].set_title('Analytic comparison bounds; not confidence intervals',fontsize=8)
for a in ax.flat:a.grid(alpha=.18)
fig.suptitle('Stage B numerical diagnostics — '+status.lower(),fontsize=11)
savefig(fig,'numerics')
