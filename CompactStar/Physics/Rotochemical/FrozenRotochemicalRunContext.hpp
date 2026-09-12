#pragma once
#include <CompactStar/Physics/Rotochemical/FrozenThermalSource.hpp>
#include <CompactStar/Physics/Rotochemical/PrescribedSpinHistory.hpp>
#include <CompactStar/Physics/Rotochemical/RotochemicalReactionResponse.hpp>
#include <CompactStar/Physics/Evolution/EvolutionConfig.hpp>
#include <CompactStar/Physics/Evolution/StarContext.hpp>
#include <CompactStar/Physics/Evolution/GeometryCache.hpp>
#include <CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp>
#include <CompactStar/Physics/Driver/Thermal/PhotonCooling_Details.hpp>
#include <CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp>
#include <CompactStar/Physics/Driver/Thermal/NeutrinoCooling_Details.hpp>
#include <CompactStar/Physics/State/ThermalState.hpp>

namespace CompactStar::Physics::Rotochemical
{
// Layer 4 frozen evaluator: every snapshot retains its single semantic authority.
class FrozenReactionEvaluator final
{
    friend class FrozenRotochemicalRunContext;
  private:
    explicit FrozenReactionEvaluator(std::shared_ptr<const GlobalUrcaChannelCoefficient> owner):owner_(std::move(owner))
    {if(!owner_)throw std::runtime_error("missing Ltilde authority");owner_->RequireCurrent();for(auto p:UrcaProcesses)entries_[ProcessIndex(p)]=owner_->Entry(p);}
    RotochemicalReactionResult Evaluate(double T,const ChemicalImbalanceState& eta)const
    {
        RotochemicalReactionResult r;
        for(auto p:UrcaProcesses){const size_t i=ProcessIndex(p);const auto l=LeptonChannel(p);const double L=entries_[i].luminosity_erg_s_Kq;
          if(L==0)continue;
          const unsigned q=TemperatureExponent(p);const double xi=eta.Xi(l,T);
          const double Fm1=q==8?UrcaImbalanceFunctions::ModifiedIncrement(xi):UrcaImbalanceFunctions::DirectIncrement(xi);
          const double H=q==8?UrcaImbalanceFunctions::HM(xi):UrcaImbalanceFunctions::HD(xi);
          r.equilibrium_erg_s[i]=L*std::pow(T,q);r.increment_erg_s[i]=r.equilibrium_erg_s[i]*Fm1;r.full_erg_s[i]=r.equilibrium_erg_s[i]+r.increment_erg_s[i];
          const double rate=L/BoltzmannErgPerK*std::pow(T,q-1)*H;r.rate_count_s[ChannelIndex(l)]+=rate;r.chemical_power_MeV_s+=eta.InfinityMeV(l)*rate;
          if(!std::isfinite(r.full_erg_s[i])||!std::isfinite(rate)||!std::isfinite(r.chemical_power_MeV_s))throw std::runtime_error("nonfinite frozen reaction");
        }return r;
    }
    double Ltilde(UrcaProcess p)const{return entries_[ProcessIndex(p)].luminosity_erg_s_Kq;}
    const UrcaCoefficientEntry& Entry(UrcaProcess p)const{return entries_[ProcessIndex(p)];}
  private:
    const std::shared_ptr<const GlobalUrcaChannelCoefficient> owner_;
    std::array<UrcaCoefficientEntry,4> entries_;
};
enum class RunPurpose {ControlledTrajectory,AnalyticControl};
struct RunQualification
{
    RunPurpose purpose=RunPurpose::ControlledTrajectory;
    size_t radial_resolution=0,eos_resolution=0;
    double rho_c_g_cm3=0;
    // Authenticated paths bound both to semantic source ownership and exact bytes.
    std::string profile_path,model_path,certificate_path,eos_path,entry_manifest_path;
    std::string entry_manifest_sha256;
};
struct SecularEvaluation
{
    RotochemicalReactionResult reaction;
    RotochemicalThermalPower beta;
    SpinHistorySample spin;
    std::array<double,2> eta_dot_MeV_s{};
    double x_dot_s=0,equilibrium_x_dot_s=0;
    double Cstar_erg_K=0,Lgamma_erg_s=0,Lother_neutrino_erg_s=0,Tsurface_infinity_K=0,Pnet_erg_s=0;
};
// Sealed v1 run. Semantic owners are never replaced by public bare coefficient arrays.
class FrozenRotochemicalRunContext final
{
    friend struct FrozenContextTestAccess;
  public:
    FrozenRotochemicalRunContext(std::shared_ptr<const Analysis::ChemicalImbalanceResponse> z,
      std::shared_ptr<const Analysis::FixedBaryonNumberResponse> fixed,
      std::vector<double> validation_I,std::vector<double> numerical_goal,std::vector<double> validation_goal,
      std::shared_ptr<const GlobalUrcaChannelCoefficient> channels,
      std::shared_ptr<const FrozenThermalSource> thermal,std::shared_ptr<const PrescribedSpinHistory> spin,
      RunQualification qualification,std::shared_ptr<const RunDependencyToken> source_token)
      :z_(std::move(z)),fixed_(CopyFixed(z_,fixed)),channels_(std::move(channels)),thermal_(std::move(thermal)),spin_(std::move(spin)),qualification_(std::move(qualification)),source_token_(std::move(source_token)),reaction_(channels_)
    {
        Need(z_&&thermal_&&spin_&&source_token_,"missing run dependency");if(qualification_.purpose==RunPurpose::ControlledTrajectory)Need(dynamic_cast<const PrescribedDipoleHistory*>(spin_.get())!=nullptr,"unqualified production spin implementation");z_->RequireCurrent();fixed_->RequireCurrent();channels_->RequireCurrent();
        global_=z_->Global();Need(channels_->ChemicalDomain()==global_,"foreign Z/Ltilde owner");
        lifetime_=global_->Lifetime();provider_=lifetime_->provider;revision_=lifetime_->revision;stars_=lifetime_->stars;source_paths_=lifetime_->source_files;
        Need(provider_&&revision_&&!stars_.empty(),"incomplete chemical lifetime");revision_snapshot_=revision_->Serialize();provider_snapshot_=ProviderIdentity(*provider_);source_generation_=source_token_->generation;
        for(const auto& s:stars_){Need(bool(s),"null star owner");profile_versions_.push_back(s->Profile().Version());}
        const std::vector<BetaChannel> order{BetaChannel::Npe,BetaChannel::NpMu};Need(z_->Channels()==order,"wrong semantic Z channel order");const auto& matrix=z_->Values();ValidateSemanticZ(z_->Channels(),matrix);
        for(auto row:order)for(auto col:order)z_values_[ChannelIndex(row)][ChannelIndex(col)]=z_->PaperZ(row,col);
        Need(z_values_[0][1]!=0&&z_values_[1][0]!=0,"missing cross Z");
        // Private copy removes caller's mutable structural-result alias; Compute binds exact Z.
        w_=std::make_unique<const Analysis::RotochemicalSpinDrive>(Analysis::RotochemicalSpinDrive::Compute(z_,fixed_,validation_I,numerical_goal,validation_goal));
        Need(w_->Values().size()==2&&w_->IPhysical().size()==2,"wrong W channel count");for(auto c:order){w_values_[ChannelIndex(c)]=w_->Values()[ChannelIndex(c)];i_values_[ChannelIndex(c)]=w_->IPhysical()[ChannelIndex(c)];}
        selection_=channels_->Selection();Need(!selection_.Enabled(UrcaProcess::De)&&!selection_.Enabled(UrcaProcess::Dmu),"controlled DU disabled");
        Need(qualification_.radial_resolution==80000&&qualification_.eos_resolution==8192&&qualification_.rho_c_g_cm3==1.10e15,"unqualified background resolution/density");
        Need(fixed_->metadata.sequence_radial_resolution==80000,"structural radial provenance mismatch");
        const std::array<std::pair<std::string,std::string>,5> sources{{
          {qualification_.profile_path,"e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890"},
          {qualification_.model_path,"3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4"},
          {qualification_.certificate_path,"7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453"},
          {qualification_.eos_path,"7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a"},
          {qualification_.entry_manifest_path,qualification_.entry_manifest_sha256}}};
        for(size_t i=0;i<sources.size();++i){const auto& s=sources[i];Need(!s.first.empty()&&s.second.size()==64,"missing qualification identity");qualification_sources_.emplace_back(s.first);Need(qualification_sources_.back().sha256==s.second,"qualification source hash mismatch");if(i<4)Need(std::any_of(source_paths_.begin(),source_paths_.end(),[&](const auto& p){return std::filesystem::absolute(p)==std::filesystem::absolute(s.first);}),"qualification not in semantic source owner");}
        domain_=channels_->DomainIdentity();metric_=channels_->MetricIdentity();partition_=channels_->PartitionKm();Need(domain_==revision_->domain&&partition_==global_->Partition()&&!metric_.empty(),"foreign domain/partition/metric");
        star_=std::make_unique<Evolution::StarContext>(stars_.front()->Profile());geo_=std::make_unique<Evolution::GeometryCache>(*star_);cfg_.n_eta=2;
        Driver::Thermal::PhotonCooling::Options po;po.surface_model=Driver::Thermal::PhotonCooling::Options::SurfaceModel::EnvelopeTbTs;photon_=std::make_unique<const Driver::Thermal::PhotonCooling>(po);
        Driver::Thermal::NeutrinoCooling::Options no;no.include_direct_urca=false;no.include_modified_urca=false;no.include_pair_breaking=false;other_=std::make_unique<const Driver::Thermal::NeutrinoCooling>(no);
        Need(!spin_->Identity().empty(),"empty spin identity");
        if(qualification_.purpose==RunPurpose::ControlledTrajectory){
          Need(dynamic_cast<const PrescribedDipoleHistory*>(spin_.get())!=nullptr,"unqualified production spin implementation");
          Need(selection_.Enabled(UrcaProcess::Me)&&selection_.Enabled(UrcaProcess::Mmu),"production MU selection mismatch");
          Need(metric_=="qualified Structure-1 radial80000 canonical nu/lambda","production metric identity mismatch");
          Need(reaction_.Entry(UrcaProcess::Me).normalization_identity=="predeclared mathematical benchmark SMe=1e-51 erg cm^-3 s^-1 K^-8"&&reaction_.Entry(UrcaProcess::Mmu).normalization_identity=="predeclared mathematical benchmark SMmu=2e-51 erg cm^-3 s^-1 K^-8","production normalization identity mismatch");
        }
        spin_identity_=spin_->Identity();RequireFullCurrent();
    }
    FrozenRotochemicalRunContext(const FrozenRotochemicalRunContext&)=delete;
    FrozenRotochemicalRunContext& operator=(const FrozenRotochemicalRunContext&)=delete;
    void RequireFullCurrent()const {RequireCheapCurrent();global_->RequireCurrent();z_->RequireCurrent();fixed_->RequireCurrent();w_->RequireCurrent();channels_->RequireCurrent();thermal_->RequireDiskCurrent();for(const auto& s:qualification_sources_)s.RequireDiskCurrent();}
    void RequireCheapCurrent()const
    {
        Need(source_token_&&source_token_->alive&&source_token_->generation==source_generation_,"stale source identity token");spin_->RequireCurrent();Need(spin_->Identity()==spin_identity_,"changed spin identity");
        Need(lifetime_->provider==provider_&&lifetime_->revision==revision_&&revision_->alive&&revision_->Serialize()==revision_snapshot_&&ProviderIdentity(*provider_)==provider_snapshot_,"stale chemical revision/provider");
        Need(lifetime_->stars==stars_&&lifetime_->source_files==source_paths_,"changed lifetime owner coverage");
        for(size_t i=0;i<stars_.size();++i)Need(stars_[i]->Profile().Version()==profile_versions_[i],"stale profile");
        for(const auto& p:fixed_->metadata.sources){Need(p.eos&&SameEOS(*p.eos,p.eos_snapshot),"stale structural EOS owner");Need(&p.star->Profile()==p.profile&&p.profile->Version()==p.profile_version,"stale structural profile");if(p.first_order)Need(&p.star->RotationResponse()==p.first_order&&p.first_order->MatchesSource(p.profile,p.profile_version),"stale Hartle first order");if(p.monopole)Need(p.star->MonopoleResponse()==p.monopole&&p.monopole->MatchesSource(p.profile,p.profile_version),"stale Hartle monopole");}
        Need(star_&&geo_&&geo_->Matches(*star_)&&star_->Provenance().source==&stars_.front()->Profile(),"foreign thermal geometry");
    }
    Evolution::DriverContext DriverContext()const {RequireCheapCurrent();Evolution::DriverContext c;c.star=star_.get();c.geo=geo_.get();c.cfg=&cfg_;c.thermo=&thermal_->Table();return c;}
    void RequireOwners(const Evolution::DriverContext& c,const PrescribedSpinHistory* spin)const {RequireCheapCurrent();Need(spin==spin_.get(),"foreign shared spin owner");Need(c.star==star_.get()&&c.geo==geo_.get()&&c.thermo==&thermal_->Table()&&c.cfg==&cfg_,"foreign run context owner");}
    const std::shared_ptr<const PrescribedSpinHistory>& SpinOwner()const{return spin_;}
    double Z(BetaChannel row,BetaChannel col)const{RequireCheapCurrent();return z_values_[ChannelIndex(row)][ChannelIndex(col)];}
    double W(BetaChannel c)const{RequireCheapCurrent();return w_values_[ChannelIndex(c)];}
    double I(BetaChannel c)const{RequireCheapCurrent();return i_values_[ChannelIndex(c)];}
    double Ltilde(UrcaProcess p)const{RequireCheapCurrent();return reaction_.Ltilde(p);}
    const FrozenThermalSource& ThermalSource()const{return *thermal_;}
    const RunQualification& Qualification()const{return qualification_;}
    const UrcaCoefficientEntry& ChannelEntry(UrcaProcess p)const{RequireCheapCurrent();return reaction_.Entry(p);}
    SecularEvaluation Evaluate(double t,const Evolution::StateVector& state,const Evolution::DriverContext& ctx,const PrescribedSpinHistory* owner)const
    {
        RequireOwners(ctx,owner);Need(state.GetThermal().Size()==1&&state.GetChem().Size()==2,"wrong state layout");
        const double T=state.GetThermal().Tinf();Need(T>0&&std::isfinite(T),"invalid temperature");
        SecularEvaluation out;const auto eta=ChemicalImbalanceState::Read(state.GetChem());out.spin=spin_->Sample(t);Need(std::isfinite(out.spin.omega_rad_s)&&std::isfinite(out.spin.omega_dot_rad_s2),"invalid spin sample");
        out.reaction=reaction_.Evaluate(T,eta);out.beta=RotochemicalThermalPower::From(out.reaction);
        for(auto row:{BetaChannel::Npe,BetaChannel::NpMu}){const auto i=ChannelIndex(row);out.eta_dot_MeV_s[i]=2*w_values_[i]*out.spin.omega_rad_s*out.spin.omega_dot_rad_s2;for(auto col:{BetaChannel::Npe,BetaChannel::NpMu})out.eta_dot_MeV_s[i]-=z_values_[i][ChannelIndex(col)]*out.reaction.Rate(col);}
        auto ph=Driver::Thermal::Detail::ComputeDerived(*photon_,state,ctx);auto other=Driver::Thermal::Detail::NeutrinoCooling_Details::ComputeDerived(*other_,state,ctx);
        Need(ph.ok&&other.ok&&ph.C_star_erg_K>0,"invalid thermal authorities");out.Cstar_erg_K=ph.C_star_erg_K;out.Lgamma_erg_s=ph.L_gamma_inf_erg_s;out.Tsurface_infinity_K=ph.Tsurf_K*std::sqrt(ph.exp2nu_surf);out.Lother_neutrino_erg_s=other.L_nu_PBF_inf_erg_s;
        const double eq=out.reaction.EquilibriumErgPerSecond();out.Pnet_erg_s=out.beta.incremental_beta_erg_s-eq-out.Lgamma_erg_s-out.Lother_neutrino_erg_s;
        out.equilibrium_x_dot_s=(-eq-out.Lgamma_erg_s-out.Lother_neutrino_erg_s)/(T*out.Cstar_erg_K);out.x_dot_s=out.Pnet_erg_s/(T*out.Cstar_erg_K);
        for(double v:{out.x_dot_s,out.equilibrium_x_dot_s,out.Cstar_erg_K,out.Lgamma_erg_s,out.Lother_neutrino_erg_s,out.Tsurface_infinity_K,out.Pnet_erg_s,out.eta_dot_MeV_s[0],out.eta_dot_MeV_s[1],out.beta.heating_erg_s,out.beta.neutrino_increment_erg_s,out.beta.incremental_beta_erg_s,out.beta.full_beta_erg_s})Need(std::isfinite(v),"nonfinite secular RHS");RequireCheapCurrent();return out;
    }
  private:
    static void ValidateSemanticZ(const std::vector<BetaChannel>& channels,const Analysis::ChemicalMatrix& matrix){Need(channels==std::vector<BetaChannel>{BetaChannel::Npe,BetaChannel::NpMu},"wrong semantic Z order");Need(matrix.size()==2&&matrix[0].size()==2&&matrix[1].size()==2,"wrong semantic Z shape");}
    static void Need(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
    static std::shared_ptr<const Analysis::FixedBaryonNumberResponse> CopyFixed(const std::shared_ptr<const Analysis::ChemicalImbalanceResponse>& z,const std::shared_ptr<const Analysis::FixedBaryonNumberResponse>& p){Need(bool(z)&&bool(p),"missing chemical/fixed authority");z->RequireCurrent();const auto& owners=z->Global()->Lifetime()->stars;for(const auto& source:p->metadata.sources)Need(source.star&&std::any_of(owners.begin(),owners.end(),[&](const auto& star){return star&&star.get()==source.star;}),"foreign fixed response lifetime owner");p->RequireCurrent();return std::make_shared<const Analysis::FixedBaryonNumberResponse>(*p);}
    static bool SameEOS(const Analysis::NumberEosSource&a,const Analysis::NumberEosSource&b){return a.identity==b.identity&&a.revision==b.revision&&a.physical_domain==b.physical_domain&&a.table_path==b.table_path;}
    static std::string ProviderIdentity(const CompactStar::ILocalThermodynamicProvider& p){const auto& m=p.Metadata();std::ostringstream s;for(const auto* x:{&m.model_id,&m.model_revision,&m.particle_content,&m.coordinate_chart,&m.temperature_scope,&m.rest_mass_convention,&m.lepton_ownership,&m.smooth_domain})s<<x->size()<<':'<<*x;return s.str();}
    const std::shared_ptr<const Analysis::ChemicalImbalanceResponse> z_;
    const std::shared_ptr<const Analysis::FixedBaryonNumberResponse> fixed_;
    const std::shared_ptr<const GlobalUrcaChannelCoefficient> channels_;
    const std::shared_ptr<const FrozenThermalSource> thermal_;
    const std::shared_ptr<const PrescribedSpinHistory> spin_;
    const RunQualification qualification_;
    const std::shared_ptr<const RunDependencyToken> source_token_;
    const FrozenReactionEvaluator reaction_;
    std::shared_ptr<const Analysis::GlobalChemicalNumberResponse> global_;
    std::unique_ptr<const Analysis::RotochemicalSpinDrive> w_;
    std::shared_ptr<Analysis::ChemicalLifetime> lifetime_;
    std::shared_ptr<const CompactStar::ILocalThermodynamicProvider> provider_;
    std::shared_ptr<Analysis::ChemicalRevision> revision_;
    std::vector<std::shared_ptr<Core::NStar>> stars_;
    std::vector<uint64_t> profile_versions_;
    std::vector<std::string> source_paths_;
    std::string revision_snapshot_,provider_snapshot_,domain_,metric_,spin_identity_;
    std::vector<double> partition_;
    std::vector<FrozenSource> qualification_sources_;
    uint64_t source_generation_=0;
    std::array<std::array<double,2>,2> z_values_{};
    std::array<double,2> w_values_{},i_values_{};
    UrcaProcessSelection selection_{};
    std::unique_ptr<Evolution::StarContext> star_;
    std::unique_ptr<Evolution::GeometryCache> geo_;
    Evolution::Config cfg_;
    std::unique_ptr<const Driver::Thermal::PhotonCooling> photon_;
    std::unique_ptr<const Driver::Thermal::NeutrinoCooling> other_;
};
}
