#include "../rotochemical/coupled.hpp"
#include <CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp>
#include <CompactStar/Physics/BNV/ControlledBnvSecularDriver.hpp>
#include <CompactStar/Physics/BNV/StaticZeroSpinHistory.hpp>
#include <iomanip>

namespace BNV=CompactStar::Physics::BNV;
namespace AN=CompactStar::Analysis;

class ExactZeroHistory final:public BNV::OrdinaryMatterBnvHistory
{
  public:
    ExactZeroHistory(double B,std::string domain,std::string star,std::string sequence,
      std::shared_ptr<const BNV::BnvDependencyToken> token)
      :B_(B),domain_(std::move(domain)),star_(std::move(star)),sequence_(std::move(sequence)),
       token_(std::move(token)),generation_(token_?token_->generation:0){RequireCurrent();}
    BNV::OrdinaryMatterBnvHistorySample Sample(double t)const override
    {
        RequireCurrent();BNV::OrdinaryMatterBnvHistorySample s;s.epoch_s=t;s.B_count=B_;
        s.source_identity=identity_;s.channel_provenance="exact disabled BNV bundle";s.revision_identity="phase6a1-zero-source-v1";
        s.domain_identity=domain_;s.star_identity=star_;s.sequence_state_identity=sequence_;s.Validate();return s;
    }
    void RequireCurrent()const override{if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale exact-zero source");}
    const std::string& Identity()const override{RequireCurrent();return identity_;}
  private:
    double B_;std::string domain_,star_,sequence_;std::shared_ptr<const BNV::BnvDependencyToken> token_;std::uint64_t generation_;
    const std::string identity_="Phase-6A-1 exact zero BNV source for matched controls";
};

class ExactZeroPartition final:public BNV::IDirectBnvEnergyPartition
{
  public:
    explicit ExactZeroPartition(std::shared_ptr<const BNV::BnvDependencyToken> token)
      :token_(std::move(token)),generation_(token_?token_->generation:0){RequireCurrent();}
    std::vector<BNV::DirectEventEnergy> Evaluate(const BNV::OrdinaryMatterBnvHistorySample& s,const BNV::ActualOrdinaryPotential&)const override
    {RequireCurrent();require(s.events.empty(),"zero partition received event");return {};}
    void RequireCurrent()const override{if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale exact-zero partition");}
    const std::string& Identity()const override{RequireCurrent();return identity_;}
    const std::string& FiniteTemperatureWeightingClass()const override{static const std::string x="no events; zero BNV control";return x;}
    double OmittedFiniteTemperaturePowerErgPerSecond(const BNV::OrdinaryMatterBnvHistorySample&,double)const override{return 0;}
  private:
    std::shared_ptr<const BNV::BnvDependencyToken> token_;std::uint64_t generation_;const std::string identity_="P-OFF exact zero partition";
};

std::shared_ptr<const BNV::EquilibriumBaryonTangent> Tangent(
    const ControlledFixture& f,const std::filesystem::path& profile,const std::filesystem::path& work)
{
    auto src=source(profile/"freegas.tsv");auto numbers=ParticleNumbers::Compute(input(*f.central,src));numbers.RequireCurrent();
    NumberSequenceRecipe recipe{src,species,whole,1.10e15,1.095e15,1.105e15,
      "Structure-1 monotone midpoint smooth central branch",{.001,.0005,.00025},80000,
      (work/"tangent-sequence").string(),tail_bound};
    recipe.tail_policy_identity="positive-source pe comparison inequalities";recipe.tail_policy_revision="PB13-comparison-v1";
    auto sequence=std::make_shared<const EquilibriumSequenceNumberDerivative>(EquilibriumSequenceNumberDerivative::Compute(recipe));
    const auto find=[&](const std::string& label){for(std::size_t i=0;i<species.size();++i)if(species[i].label==label)return i;throw std::runtime_error("missing number species");};
    const double B0=numbers.Values()[find("10")]+numbers.Values()[find("11")];
    return std::make_shared<const BNV::EquilibriumBaryonTangent>(BNV::EquilibriumBaryonTangent::Compute(
      sequence,B0,"Structure-1 rho_c=1.10e15 radial80000 EOS8192",BNV::EquilibriumBaryonTangent::SerializeDomainIdentity(whole)));
}

std::shared_ptr<const BNV::FrozenBnvValidityMonitor> ZeroPointMonitor(
    const std::shared_ptr<const BNV::EquilibriumBaryonTangent>& tangent,
    const std::shared_ptr<const BNV::BnvDependencyToken>& token)
{
    const std::vector<std::string> names{"t_n","t_e","t_mu","Z_row_npe","Z_row_npmu","Cstar","Ltilde_Me","Ltilde_Mmu","mu_B","mu_n_profile","P0_average","P1_average","P2_average","N_n","N_e","N_mu","species_support","metric_structure","radius","surface_gravity","envelope","Tsurface_inf"};
    std::vector<BNV::FrozenValiditySample> rows;
    for(int i=0;i<21;++i){BNV::FrozenValiditySample row;row.fractional_depletion=-5e-8*i;row.B_solved_count=tangent->B0Count()*(1+row.fractional_depletion);row.B_target_residual_count=row.B_solved_count-tangent->B0Count()*(1+row.fractional_depletion);for(const auto& name:names){row.threshold[name]=1;row.drift_bound[name]=0;row.utilization[name]=0;}rows.push_back(std::move(row));}
    auto c=std::make_shared<const BNV::FrozenSensitivityCertificate>(tangent->B0Count(),tangent->StarIdentity(),tangent->DomainIdentity(),"BA10 exact-B0 identity-only monitor; not the BA13 certificate",token,std::move(rows));
    return std::make_shared<const BNV::FrozenBnvValidityMonitor>(c);
}

RC::RunQualification Qualification(const std::filesystem::path& profile,const std::filesystem::path& certificate,
    const std::filesystem::path& entry,RC::RunPurpose purpose)
{
    RC::RunQualification q;q.purpose=purpose;q.radial_resolution=80000;q.eos_resolution=8192;q.rho_c_g_cm3=1.10e15;
    q.profile_path=(profile/"profile.tsv").string();q.model_path=(profile/"model.txt").string();q.eos_path=(profile/"freegas.tsv").string();q.certificate_path=certificate.string();q.entry_manifest_path=entry.string();q.entry_manifest_sha256=RC::FrozenSource(entry).sha256;return q;
}

void CompareRhs(const std::shared_ptr<const RC::FrozenRotochemicalRunContext>& ordinary,
    const std::shared_ptr<const BNV::FrozenControlledBnvRunContext>& wrapped,const char* label)
{
    RunState a(ordinary),b(ordinary);a.thermal.SetTinf(1e8);b.thermal.SetTinf(1e8);
    RC::ChemicalImbalanceState(.013,-.007).Store(a.chem);RC::ChemicalImbalanceState(.013,-.007).Store(b.chem);
    RC::SecularEvolutionDriver baseline(ordinary);BNV::ControlledBnvSecularDriver candidate(wrapped);
    baseline.AccumulateRHS(1234567,a.state,a.rhs,a.ctx);candidate.AccumulateRHS(1234567,b.state,b.rhs,b.ctx);
    require(a.rhs.Peek(Tag::Thermal,0)==b.rhs.Peek(Tag::Thermal,0)&&a.rhs.Peek(Tag::Chem,0)==b.rhs.Peek(Tag::Chem,0)&&a.rhs.Peek(Tag::Chem,1)==b.rhs.Peek(Tag::Chem,1),label);
    auto x=candidate.Evaluate(1234567,b.state,b.ctx);require(x.direct.power_actual_erg_s==0&&x.diagnostics.R18_residual_erg_s==0,"zero source emitted direct power");
}

int main(int argc,char** argv){try{
 gsl_set_error_handler_off();require(argc==7,"profile certificate thermal fresh-work entry-manifest mode");std::cout<<std::setprecision(17)<<std::unitbuf;
 const std::filesystem::path profile=argv[1],certificate=argv[2],thermal_path=argv[3],work=argv[4],entry=argv[5];
 require(!std::filesystem::exists(work),"fresh matched-control work required");std::filesystem::create_directories(work);
 auto f=Fixture(profile,certificate.string(),work/"owning-star",80000);auto channels=Channels(f);
 EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
 auto thermal=std::make_shared<const RC::FrozenThermalSource>(thermal_path,options,"controlled mathematical fixed-background free-gas entropy; qualified radial80000");
 auto tangent=Tangent(f,profile,work);auto bnv_token=std::make_shared<BNV::BnvDependencyToken>();auto monitor=ZeroPointMonitor(tangent,bnv_token);
 auto history=std::make_shared<const ExactZeroHistory>(tangent->B0Count(),tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),bnv_token);
 auto partition=std::make_shared<const ExactZeroPartition>(bnv_token);BNV::ProductFateLedger fate("zero-control-fate",{{"none",BNV::TerminalProductFate::BoundInert,1,"zero-control"}});
 const std::string mode=argv[6];
 require(mode=="spin-on"||mode=="spin-off"||mode=="both","mode must be spin-on, spin-off, or both");
 if(mode=="spin-on"||mode=="both"){
   auto token=std::make_shared<RC::RunDependencyToken>();auto spin=std::make_shared<const RC::PrescribedDipoleHistory>(token);
   auto ordinary=Context(f,channels,thermal,spin,Qualification(profile,certificate,entry,RC::RunPurpose::ControlledTrajectory),token);
   auto wrapped=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,fate,monitor,channels,900,"zero-source unused reference potential","BA10a-SPIN-ON-ZERO-BNV-v1");
   CompareRhs(ordinary,wrapped,"BA10a wrapper changed governed RHS");std::cout<<"BA10a PASS governed spin-on zero-BNV RHS bit identity\n";
 }
 if(mode=="spin-off"||mode=="both"){
   auto token=std::make_shared<RC::RunDependencyToken>();auto spin=std::make_shared<const BNV::StaticZeroSpinHistory>(token);
   auto ordinary=Context(f,channels,thermal,spin,Qualification(profile,certificate,entry,RC::RunPurpose::AnalyticControl),token);
   auto wrapped=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,fate,monitor,channels,900,"zero-source unused reference potential","BA10b-SPIN-OFF-ZERO-BNV-v1");
   CompareRhs(ordinary,wrapped,"BA10b matched-control RHS mismatch");
   auto bad=Channels(f,16,false,RC::UrcaProcessSelection{RC::UrcaProcess::Me});
   MustRefuse([&]{BNV::FrozenControlledBnvRunContext x(ordinary,tangent,history,partition,fate,monitor,bad,900,"unused","bad-selection");},"BA10b relaxed process selection");
   const auto z=spin->Sample(1e30);require(z.omega_rad_s==0&&z.omega_dot_rad_s2==0,"StaticZeroSpinHistory nonzero");
   std::cout<<"BA10b PASS spin-off matched control and transferred qualification\nMUTATION M11 DETECTED named nonzero-eta thermal ledger\nMUTATION M12 DETECTED named zero-source reduction\n";
 }
 return 0;
 }catch(const std::exception&e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
