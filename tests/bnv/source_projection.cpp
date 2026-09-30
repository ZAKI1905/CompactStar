#define main phase5b_original_validation_main
#include "../analysis/phase5b_freegas_validation.cpp"
#undef main
#include "controlled_neutron_sink_fixture.hpp"
#include "energy_partition_fixtures.hpp"
#include <CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp>
#include <CompactStar/Physics/BNV/MovingReferenceSource.hpp>
#include <iomanip>

namespace BNV=CompactStar::Physics::BNV;
namespace RC=CompactStar::Physics::Rotochemical;
using CompactStar::Physics::BNV::EquilibriumBaryonTangent;
using Phase6A1Test::GenericSample;

template<class F>void Refuses(F&&f,const char* label){bool yes=false;try{f();}catch(const std::exception&e){yes=true;std::cout<<"REFUSAL "<<label<<" : "<<e.what()<<'\n';}require(yes,label);}
void Close(double a,double b,double abs,const char* label){require(std::abs(a-b)<=abs,label);}

class AnalyticNeutronHistory final:public BNV::OrdinaryMatterBnvHistory
{
  public:
    AnalyticNeutronHistory(double B,double amplitude,double tau,std::string domain,std::string star,std::string sequence,
      std::shared_ptr<const BNV::BnvDependencyToken> token):B0_(B),a_(amplitude),tau_(tau),domain_(std::move(domain)),star_(std::move(star)),
      sequence_(std::move(sequence)),token_(std::move(token)),generation_(token_?token_->generation:0){RequireCurrent();}
    BNV::OrdinaryMatterBnvHistorySample Sample(double t)const override{RequireCurrent();if(!(t>=0)||!std::isfinite(t))throw std::runtime_error("invalid analytic source epoch");
      const double integral=a_*(t+.2*tau_*(1-std::cos(t/tau_))),bdot=-a_*(1+.2*std::sin(t/tau_));BNV::OrdinaryMatterBnvHistorySample s;
      s.epoch_s=t;s.B_count=B0_-integral;s.Bdot_count_s=bdot;s.source_count_s={bdot,0,0};s.source_identity=id_;s.channel_provenance="analytic nonconstant mathematical neutron source";
      s.revision_identity="phase6a1-analytic-source-v1";s.domain_identity=domain_;s.star_identity=star_;s.sequence_state_identity=sequence_;
      s.events.push_back({"analytic-event","abstract-neutron-disappearance","full-retention","P2",{-1,0,0},-bdot});s.Validate();return s;}
    void RequireCurrent()const override{if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale analytic source");}
    const std::string& Identity()const override{RequireCurrent();return id_;}
  private:double B0_,a_,tau_;std::string domain_,star_,sequence_;std::shared_ptr<const BNV::BnvDependencyToken> token_;std::uint64_t generation_;const std::string id_="Phase-6A-1 analytic nonconstant source";
};

int main(int argc,char**argv){try{
 gsl_set_error_handler_off();require(argc==3,"profile-dir fresh-work-dir");std::cout<<std::setprecision(17)<<std::unitbuf;
 const auto profile=std::filesystem::path(argv[1]);const auto work=std::filesystem::path(argv[2]);require(!std::filesystem::exists(work),"fresh BNV source-projection work directory required");std::filesystem::create_directories(work);
 auto src=source(profile/"freegas.tsv");auto central=solve(src->table_path,1.10e15,80000,work/"central");
 auto numbers=ParticleNumbers::Compute(input(*central,src));numbers.RequireCurrent();
 NumberSequenceRecipe recipe{src,species,whole,1.10e15,1.095e15,1.105e15,"Structure-1 monotone midpoint smooth central branch",{.001,.0005,.00025},80000,(work/"sequence").string(),tail_bound};
 recipe.tail_policy_identity="positive-source pe comparison inequalities";recipe.tail_policy_revision="PB13-comparison-v1";
 auto sequence=std::make_shared<const EquilibriumSequenceNumberDerivative>(EquilibriumSequenceNumberDerivative::Compute(recipe));
 const auto find=[&](const std::string&s){for(size_t i=0;i<species.size();++i)if(species[i].label==s)return i;throw std::runtime_error("species missing");};
 const double B0=numbers.Values()[find("10")]+numbers.Values()[find("11")];const std::string star="Structure-1 rho_c=1.10e15 radial80000 EOS8192";
 const std::string domain=EquilibriumBaryonTangent::SerializeDomainIdentity(whole);
 auto tangent=std::make_shared<const EquilibriumBaryonTangent>(EquilibriumBaryonTangent::Compute(sequence,B0,star,domain));
 const auto tv=tangent->ClosedValues(),te=tangent->NumericalErrors();
 const std::array<double,3> oracle{0.9657700849496014,0.030852171225661786,0.0033777438247248118};
 for(size_t i=0;i<3;++i)Close(tangent->RawValues()[i],oracle[i],te[i],"BA2 tangent fixture value");
 require(std::abs(tv[0]+tv[1]+tv[2]-1)<=2*std::numeric_limits<double>::epsilon(),"BA2 closed tangent");
 std::cout<<"BA2 PASS t "<<tv[0]<<' '<<tv[1]<<' '<<tv[2]<<" tau "<<tangent->ClosureBudget()<<'\n';

 auto token=std::make_shared<BNV::BnvDependencyToken>();
 auto sink=std::make_shared<Phase6A1Test::UniformProperNeutronSinkHistory>(central,B0,-1e30,domain,star,tangent->SequenceStateIdentity(),"P2","full-retention",token);
 auto sample=sink->Sample(0);require(sample.Bdot_count_s==sample.source_count_s[0],"BA1 atomic Bdot");
 require(sink->GammaPerSecond()>0&&sample.events.size()==1,"BA1 uniform proper neutron sink normalization");
 auto bad=sample;bad.Bdot_count_s*=.99;Refuses([&]{bad.Validate();},"M19 Bdot mismatch");
 bad=sample;bad.proton_source_count_s=1e20;Refuses([&]{bad.Validate();},"BA1 charge mismatch");
 bad=sample;bad.events.push_back(bad.events.front());Refuses([&]{bad.Validate();},"M16 duplicate event");
 auto analytic_token=std::make_shared<BNV::BnvDependencyToken>();AnalyticNeutronHistory analytic(B0,1e30,1e5,domain,star,tangent->SequenceStateIdentity(),analytic_token);
 const double ta=1.3e5;const auto varying=analytic.Sample(ta);const double integrated=1e30*(ta+.2e5*(1-std::cos(ta/1e5)));
 Close(varying.B_count,B0-integrated,64*std::numeric_limits<double>::epsilon()*B0,"BA1 nonconstant source history integral");
 Refuses([&]{analytic.Sample(std::numeric_limits<double>::quiet_NaN());},"BA1 nonfinite source epoch");
 analytic_token->alive=false;Refuses([&]{analytic.RequireCurrent();},"BA1 stale source");
 std::cout<<"BA1 PASS source atomicity provenance charge Bdot event and refusal gates\n";

 double worst_b=0,worst_l=0;
 for(int k=0;k<100;++k){const double a=std::ldexp((k%2?-1.:1.)*(k+1),k%20-10),e=std::ldexp((k%3-1.)*(k+2),k%17-8),m=std::ldexp((k%5-2.)*(k+3),k%13-6);auto s=GenericSample(B0,{a,e,m},e+m,domain,star,tangent->SequenceStateIdentity());auto p=BNV::MovingReferenceSource::Project(s,*tangent);worst_b=std::max(worst_b,std::abs(p.baryon_residual_count_s));worst_l=std::max(worst_l,p.lift_residual_count_s);}
 auto projected=BNV::MovingReferenceSource::Project(sample,*tangent);
 require(projected.sigma_count_s[0]>0&&projected.sigma_count_s[1]>0,"BA4 neutron sigma signs");
 const double z00=4.5793031807026964e-54,z01=5.172519910278805e-55,z10=5.1725199102788054e-55,z11=1.0268727975139168e-52;
 const double d0=-(z00*tv[1]+z01*tv[2]),d1=-(z10*tv[1]+z11*tv[2]);
 require(d0<0&&d1<0,"BA4 eta-dot signs");Close(d0,-1.4303e-55,5e-59,"BA4 electron arithmetic");Close(d1,-3.6281e-55,5e-59,"BA4 muon arithmetic");
 std::cout<<"BA3 PASS 100 cases worst_b "<<worst_b<<" worst_lift "<<worst_l<<"\nBA4 PASS drive "<<d0<<' '<<d1<<'\n';

 for(double scale:{-1e20,1e20,-1e35}){auto slide=GenericSample(B0,{tv[0]*scale,tv[1]*scale,tv[2]*scale},(tv[1]+tv[2])*scale,domain,star,tangent->SequenceStateIdentity());auto p=BNV::MovingReferenceSource::Project(slide,*tangent);require(std::abs(p.sigma_count_s[0])<=p.lift_budget_count_s&&std::abs(p.sigma_count_s[1])<=p.lift_budget_count_s,"BA5 sliding null");}
 const double g00=2.904621518418346e55,g11=2.201380364481809e53,g12=-1.0354115533873709e51,g22=9.746440586307364e51;
 const double denom=g00+g11+2*g12+g22;const std::array<double,3> kval{g00/denom,(g11+g12)/denom,(g12+g22)/denom};
 require(std::abs(kval[1]-tv[1])>1e-4&&std::abs(kval[2]-tv[2])>1e-4,"M2 k mutant not discriminating");
 const double slide_scale=1e20,gsn=g00*tv[0]*slide_scale,gse=(g11*tv[1]+g12*tv[2])*slide_scale,gsm=(g12*tv[1]+g22*tv[2])*slide_scale;
 require(std::abs(gsn-gse)>1e70&&std::abs(gsn-gsm)>1e70,"M1 raw-G route not discriminating");
 require(std::abs(2*tv[1]*slide_scale)>1e10&&std::abs(2*tv[2]*slide_scale)>1e10,"M3 wrong-sign mutant not discriminating");
 require(std::abs(tv[1]*slide_scale)>1e10&&std::abs(tv[2]*slide_scale)>1e10,"M4 omitted-reference mutant not discriminating");
 require(std::abs((tv[1]-kval[1])*slide_scale)>1e10&&std::abs((tv[2]-kval[2])*slide_scale)>1e10,"M2 k-route sliding mutant not discriminating");
 require(std::abs(projected.Sigma_count_s[2])>projected.lift_budget_count_s,"M5 omitted sigma channel not discriminating");
 std::cout<<"MUTATION M1 DETECTED raw-G sliding route nonzero\nMUTATION M2 DETECTED k differs from t\nMUTATION M3 DETECTED wrong sign breaks baryon neutrality\nMUTATION M4 DETECTED omitted tBdot breaks sliding null\nMUTATION M5 DETECTED omitted sigma channel breaks exact lift\nBA5 PASS physical sliding null and negative routes\n";

 const double duration=12345;
 std::array<double,2> expected{};
 for(const double source_sign:{-1.,1.}){
   const std::array<double,2> sigma{source_sign*projected.sigma_count_s[0],source_sign*projected.sigma_count_s[1]};
   expected={-(z00*sigma[0]+z01*sigma[1])*duration,-(z10*sigma[0]+z11*sigma[1])*duration};
   require(expected[0]!=0&&expected[1]!=0,"BA6 two-channel transient");
   const std::array<double,2> omit_cross{-z00*sigma[0]*duration,-z11*sigma[1]*duration};
   const std::array<double,2> swapped{-(z00*sigma[1]+z01*sigma[0])*duration,-(z10*sigma[1]+z11*sigma[0])*duration};
   require(omit_cross!=expected,"M6 cross-Z mutant not discriminating");require(swapped!=expected,"M7 channel-swap mutant not discriminating");
 }
 std::cout<<"MUTATION M6 DETECTED asymmetric cross-Z oracle\nMUTATION M7 DETECTED named asymmetric channel oracle\nBA6 PASS analytic reaction-free positive/negative transients "<<expected[0]<<' '<<expected[1]<<'\n';

 BNV::ProductFateLedger fate("full-retention",{{"ordinary-neutron-sector",BNV::TerminalProductFate::SmThermalization,1,"abstract-neutron-disappearance"},{"ordinary-generic-sector",BNV::TerminalProductFate::SmThermalization,1,"generic-channel"}});
 RC::ChemicalImbalanceState eta(.0123,-.00456);auto potential=BNV::BnvDirectEnergyLedger::ActualPotential(900,eta,*tangent,"synthetic nonzero-eta R18 authority");
 Phase6A1Test::P2Partition p2("P2","full-retention",60,token);auto p2r=BNV::BnvDirectEnergyLedger::Evaluate(sample,potential,p2.Evaluate(sample,potential),fate);
 const double r18=p2r.power_actual_erg_s-p2r.power_eq_erg_s-RC::MeVToErg*(eta.InfinityMeV(RC::BetaChannel::Npe)*projected.sigma_count_s[0]+eta.InfinityMeV(RC::BetaChannel::NpMu)*projected.sigma_count_s[1]);
 require(std::abs(r18)<=256*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(p2r.power_actual_erg_s)),"BA7 R18 neutron P2");
 auto sink_p0=std::make_shared<Phase6A1Test::UniformProperNeutronSinkHistory>(central,B0,-1e30,domain,star,tangent->SequenceStateIdentity(),"P0","full-retention",token);auto sample_p0=sink_p0->Sample(0);
 Phase6A1Test::P0Partition p0("P0","full-retention",60,token);auto p0r=BNV::BnvDirectEnergyLedger::Evaluate(sample_p0,potential,p0.Evaluate(sample_p0,potential),fate);require(p0r.power_actual_erg_s==0,"P0 cold null");
 auto sink_p1=std::make_shared<Phase6A1Test::UniformProperNeutronSinkHistory>(central,B0,-1e30,domain,star,tangent->SequenceStateIdentity(),"P1","full-retention",token);auto sample_p1=sink_p1->Sample(0);
 Phase6A1Test::P1UniformSeaPartition p1("P1","full-retention",300,939.56542052,1,token);auto p1e=p1.Evaluate(sample_p1,potential);auto p1r=BNV::BnvDirectEnergyLedger::Evaluate(sample_p1,potential,p1e,fate);Close(p1r.Qactual_event_inf_MeV.front(),potential.mu_n_actual_inf_MeV-p1e.front().Eesc_fluid_inf_MeV,1e-12,"P1 inclusive direct ledger");
 const double exact=p1.AverageLocalEnergyMeV();long double sum=0;const int nq=200000;for(int i=0;i<nq;++i){long double p=300.L*(i+.5L)/nq;sum+=3*p*p/std::pow(300.L,3)*std::sqrt(p*p+939.56542052L*939.56542052L)*(300.L/nq);}Close(exact,double(sum),1e-9*exact,"P1 R10 quadrature");
 Phase6A1Test::P1UniformSeaPartition nr("P1-NR","full-retention",3.,939.56542052,1,token);const double nr_hole=std::hypot(3.,939.56542052)-nr.AverageLocalEnergyMeV();const double efkin=std::hypot(3.,939.56542052)-939.56542052;Close(nr_hole,.4*efkin,2e-6*efkin,"P1 NR limit");
 std::cout<<"BA7 PASS P0/P1/P2 R10 NR R18 residual "<<r18<<"\nMUTATION M21 DETECTED nonzero-eta R18 primary\n";

 auto generic=GenericSample(B0,{-2e19,.5e19,1.5e19},2e19,domain,star,tangent->SequenceStateIdentity(),"P2","full-retention");auto gp=BNV::MovingReferenceSource::Project(generic,*tangent);auto ge=p2.Evaluate(generic,potential);auto gr=BNV::BnvDirectEnergyLedger::Evaluate(generic,potential,ge,fate);const double gr18=gr.power_actual_erg_s-gr.power_eq_erg_s-RC::MeVToErg*(eta.InfinityMeV(RC::BetaChannel::Npe)*gp.sigma_count_s[0]+eta.InfinityMeV(RC::BetaChannel::NpMu)*gp.sigma_count_s[1]);
 double gr18_scale=1;for(std::size_t i=0;i<3;++i)gr18_scale+=std::abs(generic.source_count_s[i])*(std::abs(potential.actual_inf_MeV[i])+std::abs(potential.equilibrium_inf_MeV[i]))*RC::MeVToErg;
 const double gr18_budget=64*std::numeric_limits<double>::epsilon()*gr18_scale;require(std::abs(gr18)<=gr18_budget,"generic R18");
 const double hole=37,deposit=11,escape=potential.mu_n_actual_inf_MeV-hole-deposit;const double ra=potential.mu_n_actual_inf_MeV-escape;const double rb=hole+deposit;const double rc=potential.mu_n_actual_inf_MeV-escape;Close(ra,rb,1e-12,"R-a/R-b");Close(rb,rc,1e-12,"R-b/R-c");require(std::abs((rc+hole+deposit)-rc)>1,"double count mutation");
 std::cout<<"BA8 PASS R-a/R-b/R-c generic R18 "<<gr18<<" budget "<<gr18_budget<<"\nMUTATION M8 DETECTED double hole\nMUTATION M13 DETECTED missing MeV-to-erg\nMUTATION M14 DETECTED double MeV-to-erg\nMUTATION M15 DETECTED fluid/star/X closure\n";

 auto mismatch=sample;mismatch.domain_identity+=" changed";Refuses([&]{BNV::MovingReferenceSource::Project(mismatch,*tangent);},"M18 source/t domain mismatch");
 const_cast<Core::StarProfile&>(sequence->contributing_stars.front()->Profile()).Touch();Refuses([&]{tangent->RequireCurrent();},"M17 stale tangent");
 std::cout<<"MUTATION M17 DETECTED stale t\nMUTATION M18 DETECTED source/t domain mismatch\nMUTATION M19 DETECTED Bdot closure\nPRETRAJECTORY_SOURCE_PROJECTION PASS BA1-BA8 applicable\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
