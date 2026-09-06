// ADR-0012 U9: accepted D-prime checks on the coherent fixed-epsilon_c CMF
// star. Independent (mhat,hhat) and Stieltjes oracles; no production RHS edits.
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "tests/rotation/eos_table_knots.hpp"
#include "tests/rotation/hartle_monopole_reference.hpp"
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <stdexcept>
using namespace CompactStar::Core;
using namespace hartle_mono_ref;
void need(bool x,const char* why){if(!x)throw std::runtime_error(why);}
int main(int argc,char**argv){try{
 need(argc==3,"usage: unit_monopole EOS_TABLE FRESH_WORK_DIR");gsl_set_error_handler_off();
 std::filesystem::path out=argv[2];need(!std::filesystem::exists(out),"fresh directory");std::filesystem::create_directories(out);
 auto knots=eos_knots::Read(argv[1]);need(knots.ok,"EOS knots");
 std::cout<<std::setprecision(17);std::vector<double> masses,radii,envelopes;
 // Same ratified 1.6 Msun central density; no fit to corrected deltaM.
 for(int N:{5000,10000,20000,40000,80000}){
  TOVSolver tov;tov.SetWrkDir(out.string());tov.ImportEOS(argv[1],true);tov.SetRadialRes(N);std::vector<TOVPoint> pts;
  need(tov.SingleStarSolveToTOVPoints(731253342677476.125,pts)>0&&tov.LastSolveStatus()==TOVSolveStatus::SURFACE_REACHED,"complete CMF star");
  NStar ns(pts);need(ns.ComputeHartleMonopoleResponse(),"production monopole completion");const NStar&cns=ns;
  const auto& p=cns.Profile();const auto& fo=cns.RotationResponse();const auto& prod=*cns.MonopoleResponse();Background2 b;
  for(size_t j=0;j<pts.size();j++){
   b.r.push_back((*p.GetRadius())[j]);b.m.push_back((*p.GetMass())[j]);b.eps.push_back((*p.GetEnergyDensity())[j]);b.p.push_back((*p.GetPressure())[j]);b.nu.push_back((*p.GetMetricNu())[j]);b.nup.push_back((*p.GetMetricNuPrime())[j]);b.dedp.push_back((*p.GetEosDEdP())[j]);b.s.push_back(fo.omega_bar_over_Omega[j]);b.sp.push_back(fo.domega_bar_over_Omega_dr[j]);
  }
  MHOptions o;o.I_exterior=fo.I;StieltjesOptions sp;sp.refine=4;auto prof=SolveStieltjes(b,o,sp);need(prof.ok,"profile oracle");
  StieltjesOptions sk;sk.refine=2;sk.knot_p=&knots.p;sk.knot_eps=&knots.eps;auto knot=SolveStieltjes(b,o,sk);need(knot.ok,"knot oracle");
  auto nodal=Solve(b,o);need(nodal.ok,"nodal negative source");
  double envelope=0;
  auto W=[&](size_t i){return 4*M_PI*b.r[i]*b.r[i]*prod.xi0_over_Omega2[i];};
  for(size_t i=0;i+1<pts.size();i++)if(pts[i].e<1e14){size_t a=i?i-1:0,z=std::min(i+2,pts.size()-1);double wp=(W(z)-W(a))/(b.r[z]-b.r[a]);envelope+=std::abs(b.eps[i+1]-b.eps[i])*std::abs(wp)*(b.r[i+1]-b.r[i])/2;}
  const double d=prod.deltaM_over_Omega2;double ep=std::abs(d/prof.deltaM_hat-1),ek=std::abs(d/knot.deltaM_hat-1),er=std::abs(d-knot.deltaM_hat)/envelope;
  std::cout<<"DPRIME "<<N<<' '<<d<<' '<<prof.deltaM_hat<<' '<<knot.deltaM_hat<<' '<<ep<<' '<<ek<<' '<<envelope<<' '<<er<<' '<<b.r.back()<<' '<<(1-nodal.deltaM_hat/d)<<' '<<fo.I<<' '<<prod.m0_over_Omega2[pts.size()-1]<<' '<<prod.p0star_over_Omega2[pts.size()-1]<<' '<<prod.delta_p0_over_Omega2[pts.size()-1]<<' '<<prod.xi0_over_Omega2[pts.size()-1]<<' '<<prod.surface_shell_mass_over_Omega2<<std::endl;
  need(ep<=1e-6,"Dprime2 same-representation disagreement");need(ek<=1e-4,"Dprime3 knot disagreement");need(er<=1,"Dprime4 crust envelope");need(1-nodal.deltaM_hat/d>=.03,"Dprime5 nodal deficit");
  masses.push_back(d);radii.push_back(b.r.back());envelopes.push_back(envelope);
 }
 auto spread=[](const auto&v){auto mm=std::minmax_element(v.begin(),v.end());return (*mm.second-*mm.first)/v.back();};
 need(spread(masses)<=1e-4,"Dprime1 spread");need(spread(radii)<=1e-8,"Dprime5 surface spread");
 for(size_t i=1;i<envelopes.size();i++){double ratio=envelopes[i-1]/envelopes[i];std::cout<<"ENVELOPE_RATIO "<<ratio<<std::endl;need(ratio>1,"Dprime4 contraction");}
 std::cout<<"DPRIME_PASS spread="<<spread(masses)<<" surface="<<spread(radii)<<"; sequence criterion remains in full CMF independent test"<<std::endl;return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<std::endl;return 1;}}
