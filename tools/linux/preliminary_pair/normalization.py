#!/usr/bin/env python3
import pathlib,shlex,subprocess,sys,json
from decimal import Decimal as D,localcontext
from run import ROOT,BASE,RETAINED,SNAP,load_pair,digest,save
from high_precision import Coupled,rows,lift,COEFF
from thermal_enclosure import constants
from verify_chemical_reference import PI
out=pathlib.Path(sys.argv[1]);build=out/'normalization-build';build.mkdir(exist_ok=False)
s=load_pair(SNAP)[0]
flags={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in (RETAINED/'CMakeFiles/native_p2_stage.dir/flags.make').read_text().splitlines() if ' = ' in l}
link=shlex.split((RETAINED/'CMakeFiles/native_p2_stage.dir/link.txt').read_text());obj=build/'normalization.o';binary=build/'normalization'
commands=[[link[0],*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],'-c',str(ROOT/'tools/linux/preliminary_pair/normalization.cpp'),'-o',str(obj)],
[str(obj) if x=='CMakeFiles/native_p2_stage.dir/native/stage.cpp.o' else str(binary) if x=='native_p2_stage' else x for x in link]]
for cmd in commands:subprocess.run(cmd,cwd=RETAINED,check=True)
operands=out/'normalization-operands.tsv'
args=[str(binary),s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(out/'normalization-assembly'),str(operands)]
subprocess.run(args,cwd=ROOT,check=True)
with localcontext() as ctx:
 ctx.prec=110;g={r['quantity']:lift(r['value']) for r in rows(operands)}
 radius=g['R_km']*lift(1e5);mass=g['M_km']*lift(1e5)
 gravity=lift(2.99792458e10)**2*mass/(radius**2*g['expnu_surface'])/lift(1e14)
 photon=4*PI*radius**2*g['exp2nu_surface']*lift(5.670374419e-5)*lift(1e24)*gravity*(lift(2.42)*(lift(1.81)/g['expnu_boundary']).ln()).exp()
 result=dict(geometry=g,exact_frozen_operand_photon_normalization=photon,roles={})
 for role in ('source','control'):
  model=Coupled(out,role);c=constants(model)
  direct=-model.bdot*lift(rows(COEFF)[0]['mu_B_inf'])*model.erg if role=='source' else D(0)
  ddir=abs(model.pdir-direct)+D('1e-60');dph=abs(model.photon-photon)+D('1e-60')
  rhs=(ddir+dph*(lift(2.42)*c['real_box']['x_max']).exp())/(c['real_box']['T_min']*c['real_box']['C_min'])
  result['roles'][role]=dict(exact_direct_normalization=direct,direct_normalization_error=ddir,photon_normalization_error=dph,uniform_scaled_rhs_perturbation=rhs,one_year_scaled_state_bound=rhs*D(31557600)*(c['L']*D(31557600)).exp(),L=c['L'],scale=c['scale'])
 result.update(status='CONDITIONAL_NORMALIZATION_CLOSURE',scope='Exact frozen binary64 geometry-cache operands, original envelope constants and source coefficients; no static stellar/background uncertainty requalification.',assumptions=['110-digit elementary functions and PI approximation are enclosed by conservative 1e-60 erg/s padding.','Cached geometry and heat-capacity cell values remain the accepted fixed numerical inputs.'])
 result['files']={str(p):digest(p) for p in [ROOT/'tools/linux/preliminary_pair/normalization.cpp',pathlib.Path(__file__),binary,operands]};result['commands']=commands
 save(out/'normalization.json',result)
