#!/usr/bin/env python3
"""Scientific figures from saved curves; no numerical recomputation."""
import csv,json,pathlib,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=pathlib.Path(sys.argv[1]);plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
def read(name):
 with (out/name).open() as f:return {k:np.array([float(r[k]) for r in rs]) for rs in [list(csv.DictReader(f))] for k in rs[0]}
p=read('temperature.csv');s=read('curve-source.csv');c=read('curve-control.csv');t=p['year']
fig,ax=plt.subplots(2,1,figsize=(8,6),sharex=True,constrained_layout=True)
ax[0].plot(t,p['T_source']/1e6,label='P2 source',color='#b3422b');ax[0].plot(t,p['T_control']/1e6,label='Zero-source control',color='#236398')
ax[0].set_ylabel(r'$T_\infty$ [million K]');ax[0].legend()
ax[1].plot(t,p['delta_T']/1e6,color='#68439a');ax[1].set_ylabel(r'$T_{source}-T_{control}$ [million K]');ax[1].set_xlabel('Time [Julian years]')
hp=json.loads((out/'hp-source-1000000.json').read_text())
for a in ax:
 for e in hp['crossings']:a.axvline(float(e['left_s'])/31557600,color='.5',ls=':',lw=.8)
 a.grid(alpha=.18)
fig.suptitle('2,048-year mathematical P2 benchmark — preliminary Stage B\nPiecewise high-precision pair; native matrix qualification remains separate',fontsize=11)
fig.savefig(out/'temperature.png');plt.close(fig)
fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
for key,label in [('Pdir','Direct BNV'),('Lgamma','Photon'),('Lnu','Neutrino'),('LH','Chemical heating'),('DeltaLnu','Excess neutrino')]:ax[0,0].semilogy(t,s[key],label=label)
ax[0,0].set_ylabel('Source power [erg/s]');ax[0,0].legend(fontsize=8)
ax[0,1].plot(t,s['e']*1e6,label=r'$\eta_e$');ax[0,1].plot(t,s['mu']*1e6,label=r'$\eta_\mu$');ax[0,1].set_ylabel('Chemical imbalance [eV]');ax[0,1].legend()
ax[1,0].plot(t,s['LH']/s['Pdir'],label='Heating / direct BNV');ax[1,0].plot(t,(s['LH']-s['DeltaLnu'])/s['Pdir'],label='Net beta / direct BNV');ax[1,0].set_ylabel('Power ratio');ax[1,0].legend(fontsize=8)
ax[1,1].semilogy(t,s['Echem']);ax[1,1].set_ylabel('Stored chemical energy [erg]')
for a in ax.flat:a.set_xlabel('Time [Julian years]');a.grid(alpha=.18)
fig.suptitle('Stage B coupled response — preliminary, fixed accepted stellar inputs',fontsize=12)
fig.savefig(out/'physics.png');plt.close(fig)
