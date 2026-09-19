#pragma once

#include <CompactStar/Physics/BNV/DirectEnergyLedger.hpp>
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

} // namespace Phase6A1Test
