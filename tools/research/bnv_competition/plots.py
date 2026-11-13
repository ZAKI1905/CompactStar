#!/usr/bin/env python3
"""Standalone publication-candidate PDF/PNG figures from retained S1 CSVs."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from model import Model

BLUE='#28618c'; ORANGE='#c76a2b'; GOLD='#b49436'; PURPLE='#846399'; GREY='#555b60'
plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.labelsize':11,
                     'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.dpi':200,'axes.axisbelow':True})


def main():
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    m=Model(); destination=a.out/'figures';destination.mkdir(exist_ok=True)
    table=list(csv.DictReader((a.out/'stationary-full.csv').open()))
    ds=np.unique([float(r['drive_rad2_s3']) for r in table])
    # beta computed as a ratio can differ by one ulp across drive rows; row
    # ordering is fixed by scan.py's prescribed shared beta grid.
    blocks=[[r for r in table if float(r['drive_rad2_s3'])==d] for d in ds]
    bs=np.array([float(r['beta_e']) for r in blocks[0]])
    array=lambda key:np.array([[float(r[key]) for r in b] for b in blocks]).T
    X,Y=np.meshgrid(np.log10(ds),np.log10(bs))
    null=m.I[1]*m.t[0]/(m.I[0]*m.t[1])
    note='Structure-1: 0.62 M$_\odot$ free-gas benchmark · prescribed modified Urca · $q_{dep}=0$ · scientific candidate'
    def save(fig,name):
        fig.get_layout_engine().set(rect=(0,.09,1,.91))
        fig.text(.5,.015,note+'\nControlled mathematical stationary model; applicability and reachability limits are shown separately.',ha='center',fontsize=7.5,color=GREY)
        fig.savefig(destination/(name+'.pdf'),bbox_inches='tight')
        fig.savefig(destination/(name+'.png'),bbox_inches='tight')
        plt.close(fig)
    def axes_map(ax):
        ax.set_xlabel(r'$\log_{10}(\Omega|\dot\Omega|\,/\,\mathrm{s}^{-3})$')
        ax.set_ylabel(r'$\log_{10}\beta_e$; $\Gamma=\beta_e\,\Gamma_{crit,e}$')
        ax.axhline(0,color='k',lw=.8,ls='--')
        ax.axhline(np.log10(null),color='k',lw=.8,ls=':')

    fig,ax=plt.subplots(figsize=(7.7,5.5),layout='constrained')
    xx,yy=np.meshgrid(np.linspace(-16,-9,180),np.linspace(-18,-8,200))
    critical=(-2*m.I[0]/m.t[0]/m.B*m.year)*10**xx
    value=yy-np.log10(critical)
    im=ax.pcolormesh(xx,yy,value,cmap='PuOr_r',norm=TwoSlopeNorm(0,-4,4),shading='auto',rasterized=True)
    fig.colorbar(im,ax=ax,label=r'$\log_{10}\beta_e$',extend='both')
    xx1=np.linspace(-16,-9,100);gg=np.log10((-2*m.I[0]/m.t[0]/m.B*m.year)*10**xx1)
    ax.plot(xx1,gg,'k--',lw=1.1,label=r'electron null: $\beta_e=1$')
    ax.plot(xx1,gg+np.log10(null),'k:',lw=1.2,label=r'muon null: $\beta_e=1.16508$')
    threshold=json.loads((a.out/'verification.json').read_text())['thresholds']['return_to_RC_beta']
    ax.plot(xx1,gg+np.log10(threshold),'k-',lw=.7,label='fiducial return to RC luminosity (approx. locus)')
    ax.text(-14.8,-16.2,'spin-dominated\n$\eta_e,\eta_\mu>0$',color='k')
    ax.text(-13.7,-9.9,'BNV-dominated\n$\eta_e,\eta_\mu<0$',color='k')
    ax.set(xlabel=r'$\log_{10}(\Omega|\dot\Omega|\,/\,\mathrm{s}^{-3})$',ylabel=r'$\log_{10}(\Gamma\,/\,\mathrm{yr}^{-1})$',title='A · Opposing chemical drives and their cancellation loci')
    ax.legend(loc='lower right',fontsize=8,frameon=True)
    save(fig,'A_competition')

    fig,ax=plt.subplots(figsize=(7.7,5.5),layout='constrained')
    ratio=array('surface_ratio')
    im=ax.pcolormesh(X,Y,np.log10(ratio),cmap='PuOr_r',norm=TwoSlopeNorm(0,-.55,.6),shading='auto',rasterized=True)
    fig.colorbar(im,ax=ax,label=r'$\log_{10}(T_{s,RC+BNV}^{\infty}/T_{s,RC}^{\infty})$')
    c=ax.contour(X,Y,ratio,levels=[.5,.9,1.,1.5,2.,3.],colors='k',linewidths=.7)
    ax.clabel(c,levels=[.9,1.,1.5,2.,3.],fmt='%g',fontsize=8)
    axes_map(ax);ax.set_title('B · Stationary surface-temperature modification')
    save(fig,'B_temperature')

    fid=list(csv.DictReader((a.out/'fiducial.csv').open()))
    beta=np.array([float(r['beta_e']) for r in fid])
    f=lambda k:np.array([float(r[k]) for r in fid])
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for key,label,color,style in [('LH_erg_s','gross reaction heating',ORANGE,'-'),('DeltaLnu_erg_s','excess neutrinos',PURPLE,'--'),('net_weak_erg_s','net weak thermal power',BLUE,'-'),('Lgamma_erg_s','photon luminosity',GREY,':')]:
        axs[0].loglog(beta,f(key),label=label,color=color,ls=style,lw=1.7)
    axs[0].axhline(f('RC_Lgamma_erg_s')[0],color='k',ls='-.',lw=.9,label='RC-only photon luminosity')
    axs[0].set(xlabel=r'$\beta_e$',ylabel='Power at infinity [erg s$^{-1}$]',title='C · Heating and cooling at P=5 ms, Pdot=$10^{-20}$')
    axs[0].legend(fontsize=8)
    axs[1].semilogx(beta,f('eta_e_MeV')*1000,color=BLUE,label=r'$\eta_e$')
    axs[1].semilogx(beta,f('eta_mu_MeV')*1000,color=ORANGE,ls='--',label=r'$\eta_\mu$')
    axs[1].axhline(0,color=GREY,lw=.7);axs[1].axvspan(1,null,color=GOLD,alpha=.15)
    axs[1].set(xlabel=r'$\beta_e$',ylabel='Stationary imbalance [keV]',title='Imbalance reversal and the mixed-sign interval')
    axs[1].legend()
    for ax in axs: ax.axvline(1,color='k',ls='--',lw=.6);ax.axvline(null,color='k',ls=':',lw=.6)
    save(fig,'C_channels')

    fig,axs=plt.subplots(1,2,figsize=(11,4.9),layout='constrained')
    im=axs[0].pcolormesh(X,Y,np.log10(array('Tinf_K')),cmap='cividis',shading='auto',rasterized=True)
    fig.colorbar(im,ax=axs[0],label=r'$\log_{10}(T_\infty/\mathrm{K})$')
    c=axs[0].contour(X,Y,array('Tinf_K'),levels=[1e6],colors='white',linewidths=1.4);axs[0].clabel(c,fmt={1e6:'certificate floor'},fontsize=8)
    axes_map(axs[0]);axs[0].set_title('D · Temperature applicability')
    periods=np.geomspace(1,50,250);omegas=2*np.pi/(periods*.001)
    for j,label,color in [(0,'electron',BLUE),(1,'muon',ORANGE)]:
        axs[1].loglog(periods,abs(m.I[j])*omegas**2/m.N[j],label=label,color=color)
    axs[1].axhline(1,color='k',ls='--',lw=.8,label='order-unity composition change')
    axs[1].axhline(.1,color=GREY,ls=':',lw=.8,label='10% diagnostic (not a certificate)')
    axs[1].axvline(5,color=GREY,ls='-.',lw=.8)
    axs[1].set(xlabel='Prescribed spin period [ms]',ylabel=r'$|I_{\Omega,\ell}|\Omega^2/N_\ell$',title='Rotation: sub-Kepler alone is insufficient')
    axs[1].legend(fontsize=8)
    save(fig,'D_applicability')

    fig,axs=plt.subplots(1,2,figsize=(11,4.9),layout='constrained')
    for ax,key,label in [(axs[0],'tau_slow_yr',r'$\log_{10}(\tau_{slow}/\mathrm{yr})$'),(axs[1],'depletion_reach_estimate',r'$\log_{10}(\Gamma\,\max[t_{build},3\tau_{slow}])$')]:
        limits=(5,12) if key=='tau_slow_yr' else (-7,-2)
        im=ax.pcolormesh(X,Y,np.log10(array(key)),cmap='cividis',vmin=limits[0],vmax=limits[1],shading='auto',rasterized=True)
        fig.colorbar(im,ax=ax,label=label,extend='both');axes_map(ax)
    axs[0].set_title('E · Local equilibration time')
    axs[0].text(.03,.96,'Exact channel nulls: extremely slow modes\n(up to $10^{35}$ yr); narrow cores unresolved here.',transform=axs[0].transAxes,va='top',fontsize=8,bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
    axs[1].set_title('Reachability: constant-drive depletion estimate')
    c=axs[1].contour(X,Y,array('depletion_reach_estimate'),levels=[1e-6],colors='white',linewidths=1.2);axs[1].clabel(c,fmt={1e-6:'existing gate'},fontsize=8)
    save(fig,'E_reachability')

    fig,axs=plt.subplots(1,2,figsize=(11,4.9),layout='constrained')
    axs[0].semilogx(beta,f('surface_ratio'),color=BLUE,label='full polynomials')
    axs[0].semilogx(beta,f('asym_surface_ratio'),color=ORANGE,ls='--',label='two-channel large-imbalance prediction')
    axs[0].semilogx(beta,abs(1-beta)**(2/7),color=GREY,ls=':',label='one-channel cancellation scaling')
    axs[0].set(xlabel=r'$\beta_e$',ylabel=r'$T_{s,RC+BNV}^{\infty}/T_{s,RC}^{\infty}$',title='F · Full equations versus asymptotic scaling')
    axs[0].legend(fontsize=8)
    for key,label,color,ls in [('q_equals_RC_net_keV','direct equals RC-only net power',BLUE,'-'),('q_equals_net_keV','direct equals combined net power',ORANGE,'--'),('q_offsets_luminosity_change_keV','direct equals luminosity difference',PURPLE,'-.'),('shell_finiteT_scale_keV',r'$k_B T_\infty$ (shell weighting scale)',GREY,':')]:
        axs[1].loglog(beta,f(key),label=label,color=color,ls=ls)
    axs[1].set(xlabel=r'$\beta_e$',ylabel='Deposit per destroyed baryon [keV]',title='Direct deposition scales at the fiducial drive')
    axs[1].legend(fontsize=8)
    save(fig,'F_scaling_deposition')
    (a.out/'plot-receipt.json').write_text(json.dumps(dict(job=os.environ['SLURM_JOB_ID'],figures=sorted(p.name for p in destination.iterdir()),data='stationary-full.csv; fiducial.csv',status='rendered; visual inspection required'),indent=2)+'\n')


if __name__=='__main__':main()
