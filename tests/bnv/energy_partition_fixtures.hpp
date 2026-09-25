#pragma once

#include <CompactStar/Physics/BNV/DirectEnergyLedger.hpp>
#include <CompactStar/Core/NStar.hpp>
#include <CompactStar/EOS/TrackRFreeGasThermodynamics.hpp>
#include <gsl/gsl_integration.h>
#include <cmath>
#include <memory>

namespace Phase6A1Test
{
namespace BNV=CompactStar::Physics::BNV;
namespace RC=CompactStar::Physics::Rotochemical;

class PartitionBase:public BNV::IDirectBnvEnergyPartition
{
  public:
    PartitionBase(std::string id,std::string fate,double efkin,std::shared_ptr<const BNV::BnvDependencyToken> token)
      :id_(std::move(id)),fate_(std::move(fate)),efkin_(efkin),token_(std::move(token)),generation_(token_?token_->generation:0){}
    void RequireCurrent() const override{if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale partition");}
    const std::string& Identity() const override{RequireCurrent();return id_;}
  protected:
    double Floor(const BNV::OrdinaryMatterBnvHistorySample& s,double T,double factor) const
    {double rate=0;for(const auto&e:s.events)rate+=e.rate_count_s;return rate*factor*std::pow(RC::BoltzmannMeVPerK*T,2)/efkin_*RC::MeVToErg;}
    std::string id_,fate_;double efkin_;std::shared_ptr<const BNV::BnvDependencyToken> token_;std::uint64_t generation_;
};

class P0Partition final:public PartitionBase
{
  public:using PartitionBase::PartitionBase;
    std::vector<BNV::DirectEventEnergy> Evaluate(const BNV::OrdinaryMatterBnvHistorySample&s,const BNV::ActualOrdinaryPotential&p)const override
    {RequireCurrent();std::vector<BNV::DirectEventEnergy> out;for(const auto&e:s.events)out.push_back({e.event_id,id_,fate_,p.mu_n_actual_inf_MeV,p.mu_n_actual_inf_MeV,0,0,false,1,0});return out;}
    const std::string& FiniteTemperatureWeightingClass()const override{static const std::string s="P0 smooth cold moving-potential cancellation";return s;}
    double OmittedFiniteTemperaturePowerErgPerSecond(const BNV::OrdinaryMatterBnvHistorySample&s,double T)const override{return Floor(s,T,M_PI*M_PI/6);}
};

class P2Partition final:public PartitionBase
{
  public:using PartitionBase::PartitionBase;
    std::vector<BNV::DirectEventEnergy> Evaluate(const BNV::OrdinaryMatterBnvHistorySample&s,const BNV::ActualOrdinaryPotential&)const override
    {RequireCurrent();std::vector<BNV::DirectEventEnergy> out;for(const auto&e:s.events)out.push_back({e.event_id,id_,fate_,0,0,0,0,false,1,0});return out;}
    const std::string& FiniteTemperatureWeightingClass()const override{static const std::string s="P2 cold full retention with conservative smooth-weight floor";return s;}
    double OmittedFiniteTemperaturePowerErgPerSecond(const BNV::OrdinaryMatterBnvHistorySample&s,double T)const override{return Floor(s,T,M_PI*M_PI/3);}
};

class P1UniformSeaPartition final:public PartitionBase
{
  public:
    P1UniformSeaPartition(std::string id,std::string fate,double pF,double mass,double redshift,
      std::shared_ptr<const BNV::BnvDependencyToken> token)
      :PartitionBase(std::move(id),std::move(fate),std::hypot(pF,mass)-mass,std::move(token)),pF_(pF),mass_(mass),redshift_(redshift)
    {if(!(pF_>0&&mass_>0&&redshift_>0&&redshift_<=1))throw std::runtime_error("invalid P1 fixture");}
    double AverageLocalEnergyMeV()const
    {const double ef=std::hypot(pF_,mass_);return 3/(8*std::pow(pF_,3))*(pF_*ef*(2*pF_*pF_+mass_*mass_)-std::pow(mass_,4)*std::asinh(pF_/mass_));}
    std::vector<BNV::DirectEventEnergy> Evaluate(const BNV::OrdinaryMatterBnvHistorySample&s,const BNV::ActualOrdinaryPotential&)const override
    {RequireCurrent();std::vector<BNV::DirectEventEnergy> out;for(const auto&v:s.events){const double pf=v.local_neutron_density_fm3>0?197.3269804*std::cbrt(3*M_PI*M_PI*v.local_neutron_density_fm3):pF_;const double ef=std::hypot(pf,mass_);const double avg=3/(8*std::pow(pf,3))*(pf*ef*(2*pf*pf+mass_*mass_)-std::pow(mass_,4)*std::asinh(pf/mass_));const double redshift=v.local_neutron_density_fm3>0?v.exp_phi:redshift_;const double e=redshift*avg;out.push_back({v.event_id,id_,fate_,e,e,0,0,false,redshift,0});}return out;}
    const std::string& FiniteTemperatureWeightingClass()const override{static const std::string s="P1 relativistic uniform occupied Fermi sea w(p)=3p^2/pF^3";return s;}
    double OmittedFiniteTemperaturePowerErgPerSecond(const BNV::OrdinaryMatterBnvHistorySample&s,double T)const override
    {double value=0;bool spatial=false;for(const auto&e:s.events)if(e.local_neutron_density_fm3>0){spatial=true;const double pf=hbarc_*std::cbrt(3*M_PI*M_PI*e.local_neutron_density_fm3);const double efkin=std::hypot(pf,mass_)-mass_;value+=e.rate_count_s*(M_PI*M_PI/3)*std::pow(RC::BoltzmannMeVPerK*T,2)/efkin*RC::MeVToErg;}return spatial?value:Floor(s,T,M_PI*M_PI/3);}
  private:double pF_,mass_,redshift_;const double hbarc_=197.3269804;
};

// The campaign P1 owner integrates the R10 relativistic occupied-sea oracle
// over the same frozen proper-neutron event measure as the uniform sink.  A
// one-event history may therefore stay compact without replacing a nonlinear
// whole-star average by the value at a mean density.
class P1IntegratedUniformSeaPartition final:public PartitionBase
{
  public:
    P1IntegratedUniformSeaPartition(std::string id,std::string fate,
      std::shared_ptr<const CompactStar::Core::NStar> star,
      std::shared_ptr<const CompactStar::TrackRFreeGasThermodynamicProvider> provider,
      std::shared_ptr<const BNV::BnvDependencyToken> token)
      :PartitionBase(std::move(id),std::move(fate),1,std::move(token)),star_(std::move(star)),provider_(std::move(provider)),
       profile_version_(star_?star_->Profile().Version():0)
    {
        if(!star_||!provider_)throw std::runtime_error("missing P1 integrated fixture dependency");
        const auto& p=star_->Profile();auto radius=p.GetRadius(),mass=p.GetMass(),nu=p.GetMetricNu(),nb=p.GetBaryonDensity();
        auto rule=gsl_integration_glfixed_table_alloc(32);if(!rule)throw std::bad_alloc();
        long double measure=0,energy=0,inverse_kinetic=0;
        try
        {
            for(std::size_t i=0;i<radius->Size();++i)
            {
                const double left=i?(*radius)[i-1]:0,right=(*radius)[i];
                for(int j=0;j<32;++j)
                {
                    double r,w;gsl_integration_glfixed_point(left,right,j,&r,&w,rule);const double u=(r-left)/(right-left);
                    const double m=i?(*mass)[i-1]+u*((*mass)[i]-(*mass)[i-1]):(*mass)[0]*u*u*u;
                    const double phi=i?(*nu)[i-1]+u*((*nu)[i]-(*nu)[i-1]):(*nu)[0];
                    const double density=i?(*nb)[i-1]+u*((*nb)[i]-(*nb)[i-1]):(*nb)[0];
                    const double neutron=provider_->BarotropeAt(density).number_densities_fm3[0];
                    if(!(neutron>0))continue;
                    const double lapse=std::exp(phi),event_measure=w*4*M_PI*r*r/std::sqrt(1-2*m/r)*lapse*neutron*1e54;
                    const double pf=hbarc_*std::cbrt(3*M_PI*M_PI*neutron),ef=std::hypot(pf,mass_MeV_);
                    const double average=3/(8*std::pow(pf,3))*(pf*ef*(2*pf*pf+mass_MeV_*mass_MeV_)
                      -std::pow(mass_MeV_,4)*std::asinh(pf/mass_MeV_));
                    measure+=event_measure;energy+=event_measure*lapse*average;
                    inverse_kinetic+=event_measure/(ef-mass_MeV_);
                }
            }
        }
        catch(...){gsl_integration_glfixed_table_free(rule);throw;}
        gsl_integration_glfixed_table_free(rule);
        if(!(measure>0)||!std::isfinite(double(energy/measure))||!std::isfinite(double(inverse_kinetic/measure)))
            throw std::runtime_error("invalid P1 integrated fixture");
        escape_inf_MeV_=double(energy/measure);average_inverse_efkin_MeV_inv_=double(inverse_kinetic/measure);
        RequireCurrent();
    }
    void RequireCurrent()const override
    {PartitionBase::RequireCurrent();if(!star_||star_->Profile().Version()!=profile_version_)throw std::runtime_error("stale P1 integrated fixture");}
    std::vector<BNV::DirectEventEnergy> Evaluate(const BNV::OrdinaryMatterBnvHistorySample&s,const BNV::ActualOrdinaryPotential&)const override
    {RequireCurrent();std::vector<BNV::DirectEventEnergy> out;for(const auto&e:s.events)out.push_back({e.event_id,id_,fate_,escape_inf_MeV_,escape_inf_MeV_,0,0,false,1,0});return out;}
    const std::string& FiniteTemperatureWeightingClass()const override
    {static const std::string s="P1 whole-star relativistic uniform occupied Fermi sea w(p)=3p^2/pF^3";return s;}
    double OmittedFiniteTemperaturePowerErgPerSecond(const BNV::OrdinaryMatterBnvHistorySample&s,double T)const override
    {RequireCurrent();double rate=0;for(const auto&e:s.events)rate+=e.rate_count_s;return rate*(M_PI*M_PI/3)*std::pow(RC::BoltzmannMeVPerK*T,2)*average_inverse_efkin_MeV_inv_*RC::MeVToErg;}
    double EscapeInfinityMeV()const{RequireCurrent();return escape_inf_MeV_;}
  private:
    std::shared_ptr<const CompactStar::Core::NStar> star_;
    std::shared_ptr<const CompactStar::TrackRFreeGasThermodynamicProvider> provider_;
    std::uint64_t profile_version_=0;
    double escape_inf_MeV_=0,average_inverse_efkin_MeV_inv_=0;
    const double hbarc_=197.3269804,mass_MeV_=939.56542052;
};

} // namespace Phase6A1Test

