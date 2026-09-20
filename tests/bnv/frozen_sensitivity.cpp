#define main phase5b_original_validation_main
#include "../analysis/phase5b_freegas_validation.cpp"
#undef main

#include <CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp>
#include <CompactStar/Analysis/ChemicalResponse.hpp>
#include <CompactStar/Physics/BNV/FrozenBnvValidity.hpp>
#include <CompactStar/Physics/Rotochemical/GlobalUrcaChannelCoefficient.hpp>
#include <CompactStar/Physics/Rotochemical/FrozenThermalSource.hpp>
#include <CompactStar/Physics/Evolution/StarContext.hpp>
#include <CompactStar/Physics/Evolution/GeometryCache.hpp>
#include <CompactStar/Physics/Evolution/EvolutionConfig.hpp>
#include <CompactStar/Physics/State/ThermalState.hpp>
#include <CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp>
#include <CompactStar/Physics/Driver/Thermal/PhotonCooling_Details.hpp>
#include <CompactStar/Physics/Driver/Thermal/Boundary/SurfaceGravity.hpp>

#include <fstream>
#include <iomanip>

namespace BNV=CompactStar::Physics::BNV;
namespace RC=CompactStar::Physics::Rotochemical;
namespace EV=CompactStar::Physics::Evolution;
namespace TH=CompactStar::Physics::Driver::Thermal;

namespace
{
struct TargetStar
{
    double rho_c=0,B=0,target=0,residual=0,bracket_width=0;
    std::shared_ptr<CompactStar::Core::NStar> star;
    CompactStar::Analysis::ParticleNumbers numbers;
};

std::size_t Axis(const std::vector<CompactStar::Analysis::Species>& axes,const std::string& label)
{
    for(std::size_t i=0;i<axes.size();++i)if(axes[i].label==label)return i;
    throw std::runtime_error("missing particle-number axis "+label);
}

TargetStar Evaluate(
    const std::shared_ptr<CompactStar::Analysis::NumberEosSource>& source,
    double rho,const std::filesystem::path& directory,double target)
{
    TargetStar out;out.rho_c=rho;out.target=target;
    out.star=solve(source->table_path,rho,80000,directory);
    out.numbers=CompactStar::Analysis::ParticleNumbers::Compute(input(*out.star,source));
    out.numbers.RequireCurrent();
    const auto& values=out.numbers.Values();
    out.B=values[Axis(species,"10")]+values[Axis(species,"11")];
    out.residual=out.B-target;
    return out;
}

TargetStar SolveTarget(
    const std::shared_ptr<CompactStar::Analysis::NumberEosSource>& source,
    double target,double rho_hi,const TargetStar& high,
    double initial_rho_step,const std::filesystem::path& root,std::size_t target_index)
{
    std::size_t evaluation=0;
    auto at=[&](double rho){return Evaluate(source,rho,root/("target-"+std::to_string(target_index)+"-eval-"+std::to_string(evaluation++)),target);};
    TargetStar hi=high;hi.target=target;hi.residual=hi.B-target;
    if(!(hi.residual>=0))throw std::runtime_error("target upper bracket is below target");
    double step=initial_rho_step;
    TargetStar lo=at(rho_hi-step);
    while(lo.residual>0)
    {
        step*=2;
        if(!(rho_hi-step>0))throw std::runtime_error("target-B lower bracket exhausted");
        lo=at(rho_hi-step);
    }
    const double B0=high.B;
    const double tau=5e-11*B0;
    for(unsigned iteration=0;iteration<32&&hi.B-lo.B>tau;++iteration)
    {
        const double denominator=hi.residual-lo.residual;
        if(!(denominator>0))throw std::runtime_error("nonmonotone target-B bracket");
        double rho=lo.rho_c+(hi.rho_c-lo.rho_c)*(-lo.residual)/denominator;
        const double margin=64*std::numeric_limits<double>::epsilon()*std::max(std::abs(lo.rho_c),std::abs(hi.rho_c));
        if(!(rho>lo.rho_c+margin&&rho<hi.rho_c-margin)||iteration%6==5)rho=lo.rho_c+(hi.rho_c-lo.rho_c)/2;
        auto mid=at(rho);
        if(mid.residual<=0)lo=std::move(mid);else hi=std::move(mid);
    }
    const double width=hi.B-lo.B;
    if(width>tau)throw std::runtime_error("target-B bracket did not reach predeclared tolerance");
    TargetStar answer=std::abs(lo.residual)<=std::abs(hi.residual)?std::move(lo):std::move(hi);
    answer.bracket_width=width;
    if(std::abs(answer.residual)>tau)throw std::runtime_error("target-B residual did not reach predeclared tolerance");
    return answer;
}

void WriteChemicalProfile(const CompactStar::Core::NStar& star,
    const CompactStar::TrackRFreeGasThermodynamicProvider& provider,
    const std::filesystem::path& authenticated_profile,const std::filesystem::path& directory)
{
    std::filesystem::create_directories(directory);
    for(const auto& name:{"freegas.tsv","model.txt","windows.tsv"})
        std::filesystem::copy_file(authenticated_profile/name,directory/name);
    const auto& p=star.Profile();auto r=p.GetRadius(),mass=p.GetMass(),nu=p.GetMetricNu(),nb=p.GetBaryonDensity(),pressure=p.GetPressure();
    std::ofstream out(directory/"profile.tsv");
    out<<std::setprecision(17)<<"r\tm\tnu\tnB\tPgeom\tdim\tnn\tnp\tne\tnmu\tepsMeVfm3\tPMeVfm3"
      <<"\tH00\tH01\tH02\tH10\tH11\tH12\tH20\tH21\tH22"
      <<"\teps_profile_geom\tnn_profile\tnp_profile\tne_profile\tnmu_profile\n";
    for(std::size_t i=0;i<r->Size();++i)
    {
        std::array<std::array<double,3>,3> h{};int dimension=-1;const auto value=provider.BarotropeAt((*nb)[i]);
        try{std::visit([&](const auto& v){using V=std::decay_t<decltype(v)>;dimension=V::response_dimension;
          if constexpr(V::response_dimension>0)for(int a=0;a<dimension;++a)for(int b=0;b<dimension;++b)h[a][b]=v.hessian(a,b);},provider.EquilibriumAt((*nb)[i]));}
        catch(const CompactStar::EquilibriumResolutionError&){}
        out<<(*r)[i]<<'\t'<<(*mass)[i]<<'\t'<<(*nu)[i]<<'\t'<<(*nb)[i]<<'\t'<<(*pressure)[i]<<'\t'<<dimension;
        for(double x:value.number_densities_fm3)out<<'\t'<<x;
        out<<'\t'<<value.energy_density_MeV_fm3<<'\t'<<value.pressure_MeV_fm3;
        for(const auto& row:h)for(double x:row)out<<'\t'<<x;
        out<<'\t'<<(*p.GetEnergyDensity())[i];
        for(const char* label:{"10","11","0","1"}){const auto* column=p.GetSpeciesPtr(label);require(column,"missing profile species");out<<'\t'<<(*column)[i]*(*nb)[i];}
        out<<'\n';
    }
    require(bool(out),"chemical profile write failed");
}

CompactStar::Analysis::ChemicalMatrix ReadMatrix(std::istream& in,std::size_t n)
{
    CompactStar::Analysis::ChemicalMatrix a(n,std::vector<double>(n));
    for(auto& row:a)for(double& value:row)require(bool(in>>value),"incomplete chemical certificate matrix");
    return a;
}

struct ChemicalOwners
{
    std::shared_ptr<CompactStar::TrackRFreeGasThermodynamicProvider> provider;
    std::shared_ptr<CompactStar::Analysis::GlobalChemicalNumberResponse> global;
    std::shared_ptr<CompactStar::Analysis::ChemicalImbalanceResponse> z;
    std::function<std::array<double,3>(double)> interpolate;
    std::function<RC::UrcaMetric(double)> metric;
};

ChemicalOwners BuildChemicalOwners(const std::shared_ptr<CompactStar::Core::NStar>& star,
    const std::filesystem::path& profile,const std::filesystem::path& certificate)
{
    using namespace CompactStar::Analysis;
    auto provider=std::make_shared<CompactStar::TrackRFreeGasThermodynamicProvider>();
    auto lifetime=std::make_shared<ChemicalLifetime>();lifetime->provider=provider;lifetime->stars.push_back(star);
    lifetime->revision=std::make_shared<ChemicalRevision>();
    *lifetime->revision={provider->Metadata().model_id,provider->Metadata().model_revision,
      "authenticated target-B source files retained by bytes","Zaki cold fermions and AngularVelocity unit owner",
      "Phase-6A-1 achieved target-B radial80000 star","nu=Phi, one inverse lapse, proper volume once",
      "Neutron,Electron,Muon","WholeStar P=0, finite cut plus certified tail","all profile nodes, onsets, refusal edges",
      "source/ULP-derived continuous onset certificates","positive pe mass-upper radius construction",
      "Phase-6A-1 frozen certificate; Phase-5C goals and owner","current target-B star owned",true};
    lifetime->source_files={profile/"freegas.tsv",profile/"profile.tsv",profile/"model.txt",certificate,
      "CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp","CompactStar/EOS/src/LocalThermodynamics.cpp",
      "CompactStar/EOS/TrackRFreeGasThermodynamics.hpp","CompactStar/EOS/LocalThermodynamics.hpp"};
    auto r=star->Profile().GetRadius(),mass=star->Profile().GetMass(),nu=star->Profile().GetMetricNu(),nb=star->Profile().GetBaryonDensity();
    auto interpolate=[star,r,mass,nu,nb](double x){std::size_t k=std::upper_bound(r->Values().begin(),r->Values().end(),x)-r->Values().begin();
      require(x>=0&&x<=(*r)[-1],"chemical node outside target profile");if(k==0){double t=x/(*r)[0];return std::array<double,3>{(*mass)[0]*t*t*t,(*nu)[0],(*nb)[0]};}
      if(k==r->Size())return std::array<double,3>{(*mass)[-1],(*nu)[-1],(*nb)[-1]};double t=(x-(*r)[k-1])/((*r)[k]-(*r)[k-1]);return std::array<double,3>{(*mass)[k-1]+t*((*mass)[k]-(*mass)[k-1]),(*nu)[k-1]+t*((*nu)[k]-(*nu)[k-1]),(*nb)[k-1]+t*((*nb)[k]-(*nb)[k-1])};};
    ChemicalIntegrationRequest request;request.lifetime=lifetime;
    request.background=[interpolate](double x){auto v=interpolate(x);return ChemicalMetric{v[0],v[1]};};
    request.local=[provider,interpolate](double x){const auto value=provider->EquilibriumAt(interpolate(x)[2]);const std::size_t n=std::visit([](const auto& v){return std::decay_t<decltype(v)>::response_dimension;},value);require(n>0,"target certificate threshold query");auto c=ChargeNeutralNumberSusceptibility::Compute(value,CompactStar::Analysis::ChemicalMatrix(n,std::vector<double>(n)));return ChemicalLocalResponse{c.Values(),c.NumericalError(),c.Support()};};
    std::ifstream input_certificate(certificate);require(bool(input_certificate),"target chemical certificate missing");std::size_t count=0;require(bool(input_certificate>>count),"target partition count missing");
    for(std::size_t i=0;i<count;++i){ChemicalInterval interval;require(bool(input_certificate>>interval.left>>interval.right>>interval.left_continuous_onset>>interval.right_continuous_onset),"target interval incomplete");request.intervals.push_back(interval);}
    require(bool(input_certificate>>count),"target refusal count missing");
    for(std::size_t i=0;i<count;++i){ChemicalRefusalCertificate f;require(bool(input_certificate>>f.left>>f.right>>f.containing_left>>f.containing_right>>f.first_cell>>f.last_cell>>f.radius_upper>>f.mass_upper>>f.nu_lower),"target refusal incomplete");f.branch=i==0?"radial left npe; radial right pe":"radial left npemu; radial right npe";f.availability_authority="Track-R source guards and current-profile interpolation";f.geometry_extrema_authority="entire current containing profile cells";f.numerical_error=ReadMatrix(input_certificate,3);request.refusals.push_back(std::move(f));}
    ChemicalPeTailCertificate tail;require(bool(input_certificate>>tail.cut_radius>>tail.cut_mass>>tail.cut_nu>>tail.h_upper>>tail.epsilon_upper>>tail.bootstrap_radius_upper>>tail.total_mass_upper>>tail.radius_upper>>tail.susceptibility_upper>>tail.lapse_error_upper>>tail.response_upper),"target tail certificate incomplete");tail.source_hypotheses="authenticated Track-R pe positive monotone energy/pressure, current whole-star cut";request.pe_tail=tail;
    request.center_error=ReadMatrix(input_certificate,3);request.tail_error=ReadMatrix(input_certificate,3);request.background_error=ReadMatrix(input_certificate,3);request.absolute_goal=ReadMatrix(input_certificate,3);
    auto qgoal=ReadMatrix(input_certificate,2),zgoal=ReadMatrix(input_certificate,2);
    request.center_authority="Phase-5C current-profile regular-center construction";request.tail_authority="Phase-5C current-profile pe comparison";request.background_error_authority="governed Phase-5C numerical-error components on current target profile";
    request.structural_zeros={{NumberAxis::Neutron,NumberAxis::Electron},{NumberAxis::Electron,NumberAxis::Neutron},{NumberAxis::Neutron,NumberAxis::Muon},{NumberAxis::Muon,NumberAxis::Neutron}};request.structural_zero_authority="independent free-gas neutral y-chart polynomial";
    auto global=std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(request));
    auto z=std::make_shared<ChemicalImbalanceResponse>(ChemicalImbalanceResponse::Compute(global,qgoal,zgoal));
    auto metric=[interpolate](double x){auto v=interpolate(x);return RC::UrcaMetric{v[1],x==0?0:-.5*std::log1p(-2*v[0]/x)};};
    return {provider,global,z,interpolate,metric};
}

std::shared_ptr<const RC::GlobalUrcaChannelCoefficient> BuildChannels(const ChemicalOwners& owner,unsigned order)
{
    RC::UrcaIntegrationRequest request;request.chemical_domain=owner.global;request.metric=owner.metric;request.radial_partition_km=owner.global->Partition();request.domain_identity=owner.global->Lifetime()->revision->domain;request.metric_identity="qualified Structure-1 radial80000 canonical nu/lambda";request.order=order;request.selection=RC::UrcaProcessSelection::ControlledModifiedOnly();
    auto edge=[&](double onset){double lo=0,hi=request.radial_partition_km.back();for(int i=0;i<80;++i){double mid=(lo+hi)/2;if(owner.interpolate(mid)[2]>onset)lo=mid;else hi=mid;}return (lo+hi)/2;};
    request.normalizations.emplace_back(RC::UrcaProcess::Me,[](double){return 1e-51;},std::vector<RC::UrcaSupportInterval>{{0,edge(owner.provider->NeutronOnsetBaryonDensityFm3())}},"predeclared mathematical benchmark SMe=1e-51 erg cm^-3 s^-1 K^-8");
    request.normalizations.emplace_back(RC::UrcaProcess::Mmu,[](double){return 2e-51;},std::vector<RC::UrcaSupportInterval>{{0,edge(owner.provider->MuonOnsetBaryonDensityFm3())}},"predeclared mathematical benchmark SMmu=2e-51 erg cm^-3 s^-1 K^-8");
    return std::make_shared<const RC::GlobalUrcaChannelCoefficient>(RC::GlobalUrcaChannelCoefficient::Compute(request));
}

double MuNMeV(const CompactStar::ActiveLocalThermodynamicEvaluation& value)
{
    return std::visit([](const auto& v)->double{using V=std::decay_t<decltype(v)>;
      if constexpr(std::is_same_v<V,CompactStar::LocalThermodynamicEvaluation>||std::is_same_v<V,CompactStar::NpeThermodynamicEvaluation>)return v.conjugates.MuNMeV();
      else if constexpr(std::is_same_v<V,CompactStar::MuonThresholdEvaluation>)return v.limiting_npe_conjugates.MuNMeV();
      else if constexpr(std::is_same_v<V,CompactStar::NeutronThresholdEvaluation>)return v.limiting_pe_conjugates.HPeMeV();
      else if constexpr(std::is_same_v<V,CompactStar::PeThermodynamicEvaluation>)return v.conjugates.HPeMeV();
      else return std::numeric_limits<double>::quiet_NaN();},value);
}

struct DirectReference
{
    double mu_B_inf_MeV=0,p0_inf_MeV=0,p1_inf_MeV=0,p2_inf_MeV=0;
    double neutron_measure=0;
};

double RelativisticUniformSeaAverage(double pf,double mass)
{
    const double x=pf/mass;
    if(x<1e-3){const double x2=x*x;return mass*(1+3*x2/10-3*x2*x2/56+x2*x2*x2/48);}
    const double ef=std::hypot(pf,mass);
    return 3/(8*std::pow(pf,3))*(pf*ef*(2*pf*pf+mass*mass)-std::pow(mass,4)*std::asinh(pf/mass));
}

DirectReference DirectReferenceAverages(const CompactStar::Core::NStar& star,
    const CompactStar::TrackRFreeGasThermodynamicProvider& provider,unsigned order)
{
    const auto& p=star.Profile();auto r=p.GetRadius(),mass=p.GetMass(),nu=p.GetMetricNu(),nb=p.GetBaryonDensity();
    auto rule=gsl_integration_glfixed_table_alloc(order);require(rule,"direct-reference quadrature allocation");long double den=0,p0=0,p1=0;
    try{for(std::size_t i=0;i<r->Size();++i){const double a=i?(*r)[i-1]:0,b=(*r)[i];for(unsigned j=0;j<order;++j){double x,w;gsl_integration_glfixed_point(a,b,j,&x,&w,rule);const double u=(x-a)/(b-a);
      const double mm=i?(*mass)[i-1]+u*((*mass)[i]-(*mass)[i-1]):(*mass)[0]*u*u*u;
      const double nn=i?(*nu)[i-1]+u*((*nu)[i]-(*nu)[i-1]):(*nu)[0];
      const double density=i?(*nb)[i-1]+u*((*nb)[i]-(*nb)[i-1]):(*nb)[0];const auto bar=provider.BarotropeAt(density);const double neutron=bar.number_densities_fm3[0];if(!(neutron>0))continue;
      const double lapse=std::exp(nn),measure=w*4*M_PI*x*x/std::sqrt(1-2*mm/x)*lapse*neutron*1e54;
      const double mu=MuNMeV(provider.EquilibriumAt(density));require(std::isfinite(mu),"missing neutron potential on source support");
      const double pf=197.3269804*std::cbrt(3*M_PI*M_PI*neutron),mn=939.56542052,average=RelativisticUniformSeaAverage(pf,mn);
      den+=measure;p0+=measure*lapse*mu;p1+=measure*lapse*average;}}
    }catch(...){gsl_integration_glfixed_table_free(rule);throw;}gsl_integration_glfixed_table_free(rule);
    require(den>0,"empty neutron direct-reference support");DirectReference out;out.neutron_measure=double(den);out.p0_inf_MeV=double(p0/den);out.p1_inf_MeV=double(p1/den);out.mu_B_inf_MeV=out.p0_inf_MeV;return out;
}

void WriteProfileMetrics(std::size_t index,const CompactStar::Core::NStar& star,
    const ChemicalOwners& owners,std::ofstream& output)
{
    const auto& p=star.Profile();const auto r=p.GetRadius(),mass=p.GetMass(),nu=p.GetMetricNu(),nb=p.GetBaryonDensity();
    double support_left=0,support_right=(*r)[-1];
    const double neutron_onset=owners.provider->NeutronOnsetBaryonDensityFm3();
    for(int iteration=0;iteration<80;++iteration){const double mid=(support_left+support_right)/2;
      if(owners.interpolate(mid)[2]>neutron_onset)support_left=mid;else support_right=mid;}
    const double source_radius=(support_left+support_right)/2;
    for(unsigned k=0;k<=2048;++k)
    {
        const double q=static_cast<double>(k)/2049.0,x=q*source_radius;
        const auto background=owners.interpolate(x);
        const double lambda=x==0?0:-.5*std::log1p(-2*background[0]/x);
        const auto local=owners.provider->EquilibriumAt(background[2]);
        const int dimension=std::visit([](const auto& value){return std::decay_t<decltype(value)>::response_dimension;},local);
        const double mu_n_inf=std::exp(background[1])*MuNMeV(local);
        require(std::isfinite(lambda)&&std::isfinite(mu_n_inf),"invalid frozen profile metric");
        output<<index<<'\t'<<q<<'\t'<<x<<'\t'<<background[0]<<'\t'<<background[1]<<'\t'<<lambda
          <<'\t'<<background[2]<<'\t'<<dimension<<'\t'<<mu_n_inf<<'\n';
    }
}

int RunCoefficients(const std::filesystem::path& authenticated_profile,
    const std::filesystem::path& work,int selected=-1)
{
    auto source_owner=source(authenticated_profile/"freegas.tsv");
    CompactStar::EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    std::ifstream targets(work/"target_stars.tsv");require(bool(targets),"target-star table missing");std::string line;std::getline(targets,line);
    const std::string suffix=selected<0?"":"-shard-"+std::to_string(selected);
    std::ofstream coefficients(work/("coefficients"+suffix+".tsv")),thermals(work/("thermal_surface"+suffix+".tsv")),profiles(work/("profile_metrics"+suffix+".tsv"));
    coefficients<<std::setprecision(17)<<"index\tB_solved\tZ00\tZ01\tZ10\tZ11\tuZ00\tuZ01\tuZ10\tuZ11\tLtilde_Me\tLtilde_Mmu\tuLtilde_Me\tuLtilde_Mmu\tmu_B_inf\tP0_inf\tP1_inf\tP2_inf\tu_mu_B\tu_P0\tu_P1\tu_P2\tneutron_measure\n";
    thermals<<std::setprecision(17)<<"index\tTinf_K\tCstar_erg_K\tu_Cstar_erg_K\tTsurface_local_K\tTsurface_inf_K\tu_Tsurface_inf_K\tsurface_g14\tu_surface_g14\n";
    profiles<<std::setprecision(17)<<"index\tq\tr_km\tm_km\tnu\tlambda\tnB_fm3\tactive_dimension\tmu_n_inf_MeV\n";
    for(std::size_t expected=0;std::getline(targets,line);++expected)
    {
        std::istringstream row(line);std::size_t index;double fraction,rho,B,target,residual,bracket,Nn,Ne,Nmu,uNn,uNe,uNmu,R,M,nu;
        require(bool(row>>index>>fraction>>rho>>B>>target>>residual>>bracket>>Nn>>Ne>>Nmu>>uNn>>uNe>>uNmu>>R>>M>>nu)&&index==expected,"malformed target-star table");
        if(selected>=0&&static_cast<int>(index)!=selected)continue;
        auto current=Evaluate(source_owner,rho,work/(selected<0?"coefficient-stars":"coefficient-stars-shard")/("target-"+std::to_string(index)),target);
        require(current.B==B,"coefficient star differs from achieved target-B star");
        const auto profile=work/"chemical-profiles"/("target-"+std::to_string(index));
        RC::FrozenThermalSource thermal(profile/"thermal",options,
          "current target-B mathematical free-gas entropy adapter; radial80000; Phase-5D formula unchanged");
        auto owners=BuildChemicalOwners(current.star,profile,profile/"certificate.txt");auto channels=BuildChannels(owners,16),fine_channels=BuildChannels(owners,32);
        const auto& z=owners.z->Values();const auto& uz=owners.z->NumericalError();
        const auto direct=DirectReferenceAverages(*current.star,*owners.provider,32),coarse_direct=DirectReferenceAverages(*current.star,*owners.provider,16);
        const double l_me=channels->LuminosityCoefficient(RC::UrcaProcess::Me),l_mmu=channels->LuminosityCoefficient(RC::UrcaProcess::Mmu);
        const double ul_me=std::abs(l_me-fine_channels->LuminosityCoefficient(RC::UrcaProcess::Me));
        const double ul_mmu=std::abs(l_mmu-fine_channels->LuminosityCoefficient(RC::UrcaProcess::Mmu));
        coefficients<<index<<'\t'<<B<<'\t'<<z[0][0]<<'\t'<<z[0][1]<<'\t'<<z[1][0]<<'\t'<<z[1][1]
          <<'\t'<<uz[0][0]<<'\t'<<uz[0][1]<<'\t'<<uz[1][0]<<'\t'<<uz[1][1]
          <<'\t'<<l_me<<'\t'<<l_mmu<<'\t'<<ul_me<<'\t'<<ul_mmu
          <<'\t'<<direct.mu_B_inf_MeV<<'\t'<<direct.p0_inf_MeV<<'\t'<<direct.p1_inf_MeV<<'\t'<<direct.p2_inf_MeV
          <<'\t'<<std::abs(direct.mu_B_inf_MeV-coarse_direct.mu_B_inf_MeV)
          <<'\t'<<std::abs(direct.p0_inf_MeV-coarse_direct.p0_inf_MeV)
          <<'\t'<<std::abs(direct.p1_inf_MeV-coarse_direct.p1_inf_MeV)
          <<'\t'<<std::abs(direct.p2_inf_MeV-coarse_direct.p2_inf_MeV)<<'\t'<<direct.neutron_measure<<'\n';
        WriteProfileMetrics(index,*current.star,owners,profiles);
        EV::StarContext star_context(current.star->Profile());EV::GeometryCache geometry(star_context);EV::Config config;
        TH::PhotonCooling::Options po;po.surface_model=TH::PhotonCooling::Options::SurfaceModel::EnvelopeTbTs;TH::PhotonCooling photon(po);
        const double g14=TH::Boundary::SurfaceGravity_g14(star_context,&geometry);
        for(int k=0;k<=12;++k){const double T=std::pow(10.,6+.25*k);CompactStar::Physics::State::ThermalState thermal_state;thermal_state.Resize(1);thermal_state.SetTinf(T);EV::StateVector state;state.Register(CompactStar::Physics::State::StateTag::Thermal,thermal_state);EV::DriverContext context;context.star=&star_context;context.geo=&geometry;context.cfg=&config;context.thermo=&thermal.Table();const auto d=TH::Detail::ComputeDerived(photon,state,context);require(d.ok&&d.C_star_erg_K>0,"current thermal/surface evaluation failed");const double tsinf=d.Tsurf_K*std::sqrt(d.exp2nu_surf);thermals<<index<<'\t'<<T<<'\t'<<d.C_star_erg_K<<'\t'<<64*std::numeric_limits<double>::epsilon()*d.C_star_erg_K<<'\t'<<d.Tsurf_K<<'\t'<<tsinf<<'\t'<<64*std::numeric_limits<double>::epsilon()*tsinf<<'\t'<<g14<<'\t'<<64*std::numeric_limits<double>::epsilon()*g14<<'\n';}
        std::cout<<"FROZEN_COEFFICIENTS "<<index<<" Z "<<z[0][0]<<' '<<z[0][1]<<' '<<z[1][1]<<" Ltilde "<<channels->LuminosityCoefficient(RC::UrcaProcess::Me)<<' '<<channels->LuminosityCoefficient(RC::UrcaProcess::Mmu)<<'\n';
    }
    require(bool(coefficients)&&bool(thermals)&&bool(profiles),"frozen coefficient evidence write failed");
    std::cout<<"BA13_COEFFICIENT_SAMPLES PASS count "<<(selected<0?21:1)<<" temperature_knots 13\n";return 0;
}
}

int main(int argc,char** argv)
{
 try
 {
    gsl_set_error_handler_off();
    if(argc==4&&std::string(argv[1])=="coefficients")return RunCoefficients(argv[2],argv[3]);
    if(argc==5&&std::string(argv[1])=="coefficient")return RunCoefficients(argv[2],argv[3],std::stoi(argv[4]));
    require(argc==4&&std::string(argv[1])=="stars","stars profile-dir fresh-work-dir OR coefficients profile-dir thermal-dir work-dir");
    std::cout<<std::setprecision(17)<<std::unitbuf;
    const auto profile=std::filesystem::path(argv[2]);
    const auto work=std::filesystem::path(argv[3]);
    require(!std::filesystem::exists(work),"fresh frozen-sensitivity work directory required");
    std::filesystem::create_directories(work);
    auto source_owner=source(profile/"freegas.tsv");
    CompactStar::TrackRFreeGasThermodynamicProvider thermodynamic_provider;
    const double rho0=1.10e15;
    auto reference=Evaluate(source_owner,rho0,work/"target-0-eval-0",0);
    reference.target=reference.B;reference.residual=0;reference.bracket_width=0;
    const double B0=reference.B,tau=5e-11*B0;
    std::vector<TargetStar> stars;stars.push_back(reference);
    double rho_step=rho0*1e-6;
    for(std::size_t i=1;i<=20;++i)
    {
        const double target=B0*(1-5e-8*static_cast<double>(i));
        auto solved=SolveTarget(source_owner,target,stars.back().rho_c,stars.back(),rho_step,work/"target-solves",i);
        rho_step=std::max(stars.back().rho_c-solved.rho_c,rho0*1e-12);
        std::cout<<"TARGET_B "<<i<<" fraction "<<(solved.B-B0)/B0<<" rho_c "<<solved.rho_c
                 <<" residual "<<solved.residual<<" bracket "<<solved.bracket_width<<" tau "<<tau<<'\n';
        stars.push_back(std::move(solved));
    }
    std::ofstream out(work/"target_stars.tsv");
    out<<std::setprecision(17)<<"index\ttarget_fraction\trho_c_g_cm3\tB_solved\tB_target\tresidual\tbracket_width\tN_n\tN_e\tN_mu\tu_N_n\tu_N_e\tu_N_mu\tradius_km\tmass_km\tnu_surface\n";
    const auto n=Axis(species,"10"),e=Axis(species,"0"),m=Axis(species,"1");
    for(std::size_t i=0;i<stars.size();++i)
    {
        const auto& s=stars[i];const auto& values=s.numbers.Values();const auto& errors=s.numbers.Errors();const auto& p=s.star->Profile();
        out<<i<<'\t'<<-5e-8*static_cast<double>(i)<<'\t'<<s.rho_c<<'\t'<<s.B<<'\t'<<s.target<<'\t'<<s.residual<<'\t'<<s.bracket_width
           <<'\t'<<values[n]<<'\t'<<values[e]<<'\t'<<values[m]<<'\t'<<errors[n]<<'\t'<<errors[e]<<'\t'<<errors[m]<<'\t'<<(*p.GetRadius())[-1]<<'\t'<<(*p.GetMass())[-1]<<'\t'<<(*p.GetMetricNu())[-1]<<'\n';
    }
    require(bool(out),"target-star evidence write failed");
    for(std::size_t i=0;i<stars.size();++i)WriteChemicalProfile(*stars[i].star,thermodynamic_provider,profile,work/"chemical-profiles"/("target-"+std::to_string(i)));
    std::cout<<"BA13_TARGET_STARS PASS count 21 tau_B_target "<<tau<<'\n';

    std::ofstream tangent_out(work/"tangent.tsv");
    tangent_out<<std::setprecision(17)<<"index\tB_solved\tt_n\tt_e\tt_mu\tu_t_n\tu_t_e\tu_t_mu\traw_closure\tclosure_budget\n";
    const auto domain=CompactStar::Physics::BNV::EquilibriumBaryonTangent::SerializeDomainIdentity(whole);
    for(std::size_t i=0;i<stars.size();++i)
    {
        const auto rho=stars[i].rho_c;
        NumberSequenceRecipe recipe{source_owner,species,whole,rho,rho-5e12,rho+5e12,
          "Structure-1 monotone midpoint smooth central branch",{.001,.0005,.00025},80000,
          (work/"tangent-sequences"/("target-"+std::to_string(i))).string(),tail_bound};
        recipe.tail_policy_identity="positive-source pe comparison inequalities";
        recipe.tail_policy_revision="PB13-comparison-v1";
        auto sequence=std::make_shared<const EquilibriumSequenceNumberDerivative>(EquilibriumSequenceNumberDerivative::Compute(recipe));
        auto tangent=CompactStar::Physics::BNV::EquilibriumBaryonTangent::Compute(sequence,stars[i].B,
          "Structure-1 target-B frozen certificate radial80000 EOS8192 index="+std::to_string(i),domain);
        const auto& t=tangent.ClosedValues();const auto& u=tangent.NumericalErrors();
        tangent_out<<i<<'\t'<<stars[i].B<<'\t'<<t[0]<<'\t'<<t[1]<<'\t'<<t[2]
          <<'\t'<<u[0]<<'\t'<<u[1]<<'\t'<<u[2]<<'\t'<<tangent.RawClosureResidual()<<'\t'<<tangent.ClosureBudget()<<'\n';
        std::cout<<"FROZEN_T "<<i<<' '<<t[0]<<' '<<t[1]<<' '<<t[2]<<" errors "<<u[0]<<' '<<u[1]<<' '<<u[2]<<'\n';
    }
    require(bool(tangent_out),"tangent evidence write failed");
    std::cout<<"BA13_TANGENT_SAMPLES PASS count 21\n";
    return 0;
 }
 catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
