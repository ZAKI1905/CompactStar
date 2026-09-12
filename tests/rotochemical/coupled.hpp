#include "fixture.hpp"
#include <CompactStar/Physics/Rotochemical/SecularEvolutionDriver.hpp>
#include <CompactStar/Physics/Rotochemical/ScaledRKF45.hpp>
#include <CompactStar/Physics/Evolution/StatePacking.hpp>
#include <CompactStar/Physics/Evolution/Observers/IObserver.hpp>
namespace P=CompactStar::Physics;
namespace EV=P::Evolution;
namespace TH=P::Driver::Thermal;
using Tag=P::State::StateTag;
constexpr double Year=365.25*86400;
std::shared_ptr<const RC::GlobalUrcaChannelCoefficient> Channels(const ControlledFixture& f,unsigned order=16,bool half=false,RC::UrcaProcessSelection selection=RC::UrcaProcessSelection::ControlledModifiedOnly()) {
 RC::UrcaIntegrationRequest r;r.chemical_domain=f.g;r.metric=f.metric;r.radial_partition_km=f.g->Partition();r.domain_identity=f.g->Lifetime()->revision->domain;r.metric_identity="qualified Structure-1 radial80000 canonical nu/lambda";r.order=order;r.selection=selection;
 if(half){auto v=r.radial_partition_km;r.radial_partition_km.clear();for(size_t i=1;i<v.size();++i){r.radial_partition_km.push_back(v[i-1]);r.radial_partition_km.push_back((v[i-1]+v[i])/2);}r.radial_partition_km.push_back(v.back());}
 auto edge=[&](double nb){double lo=0,hi=r.radial_partition_km.back();for(int i=0;i<80;++i){double mid=(lo+hi)/2;if(f.interpolate(mid)[2]>nb)lo=mid;else hi=mid;}return (lo+hi)/2;};
 r.normalizations.emplace_back(RC::UrcaProcess::Me,[](double){return 1e-51;},std::vector<RC::UrcaSupportInterval>{{0,edge(f.provider->NeutronOnsetBaryonDensityFm3())}},"predeclared mathematical benchmark SMe=1e-51 erg cm^-3 s^-1 K^-8");
 r.normalizations.emplace_back(RC::UrcaProcess::Mmu,[](double){return 2e-51;},std::vector<RC::UrcaSupportInterval>{{0,edge(f.provider->MuonOnsetBaryonDensityFm3())}},"predeclared mathematical benchmark SMmu=2e-51 erg cm^-3 s^-1 K^-8");
 return std::make_shared<const RC::GlobalUrcaChannelCoefficient>(RC::GlobalUrcaChannelCoefficient::Compute(r));
}

struct TestSpin final:RC::PrescribedSpinHistory {
 const double omega,slope;const std::string identity="analytic owned linear spin";
 TestSpin(double o,double d):omega(o),slope(d){}
 RC::SpinHistorySample Sample(double t)const override{return {omega+slope*t,slope};}
 void RequireCurrent()const override{}
 const std::string& Identity()const override{return identity;}
};
struct RunState {
 const std::shared_ptr<const RC::FrozenRotochemicalRunContext> context;
 EV::DriverContext ctx;P::State::ThermalState thermal;P::State::ChemState chem;
 EV::StateVector state;EV::StateLayout layout;EV::RHSAccumulator rhs;
 explicit RunState(std::shared_ptr<const RC::FrozenRotochemicalRunContext> c):context(std::move(c)),ctx(context->DriverContext()){
  thermal.Resize(1);thermal.SetTinf(1e8);chem.Resize(2);state.Register(Tag::Thermal,thermal);state.Register(Tag::Chem,chem);
  layout.Configure(state,{Tag::Thermal,Tag::Chem});rhs.Configure(Tag::Thermal,1);rhs.Configure(Tag::Chem,2);
 }
 std::vector<double> Pack(){std::vector<double> y(3);EV::PackStateVector(state,layout,y.data());return y;}
};
std::shared_ptr<const RC::FrozenRotochemicalRunContext> Context(const ControlledFixture& f,std::shared_ptr<const RC::GlobalUrcaChannelCoefficient> c,std::shared_ptr<const RC::FrozenThermalSource> thermal,std::shared_ptr<const RC::PrescribedSpinHistory> spin,RC::RunQualification q,std::shared_ptr<const RC::RunDependencyToken> token){return std::make_shared<const RC::FrozenRotochemicalRunContext>(f.z,f.fixed,f.vi,f.gn,f.gv,c,thermal,spin,q,token);}
template<class F>void MustRefuse(F f,const char* label){bool refused=false;try{f();}catch(const std::exception& e){refused=true;std::cout<<"REFUSAL "<<label<<" : "<<e.what()<<'\n';}require(refused,label);}
void Near(double a,double b,double tolerance,const char* label){require(std::abs(a-b)<=tolerance*std::max(std::abs(b),1e-300),label);}
std::vector<double> LinearTimes(double end,unsigned n=16){std::vector<double> t;for(unsigned i=1;i<=n;++i)t.push_back(end*i/n);return t;}
