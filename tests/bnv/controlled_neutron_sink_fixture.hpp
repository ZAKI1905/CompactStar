#pragma once

#include <CompactStar/Physics/BNV/OrdinaryMatterSource.hpp>
#include <CompactStar/Core/NStar.hpp>
#include <gsl/gsl_integration.h>

namespace Phase6A1Test
{
namespace BNV=CompactStar::Physics::BNV;

class ConstantNeutronSinkHistory final:public BNV::OrdinaryMatterBnvHistory
{
  public:
    ConstantNeutronSinkHistory(double B0,double Bdot,std::string domain,std::string star,
      std::string sequence,std::string partition,std::string fate,
      std::shared_ptr<const BNV::BnvDependencyToken> token)
      :B0_(B0),Bdot_(Bdot),domain_(std::move(domain)),star_(std::move(star)),sequence_(std::move(sequence)),
       partition_(std::move(partition)),fate_(std::move(fate)),token_(std::move(token)),generation_(token_?token_->generation:0)
    {if(!(B0_>0)||!(Bdot_<0)||!std::isfinite(Bdot_))throw std::runtime_error("invalid controlled neutron drive");RequireCurrent();}
    BNV::OrdinaryMatterBnvHistorySample Sample(double t) const override
    {
        RequireCurrent();
        BNV::OrdinaryMatterBnvHistorySample s;
        s.epoch_s=t;s.B_count=B0_+Bdot_*t;s.Bdot_count_s=Bdot_;s.source_count_s={Bdot_,0,0};s.proton_source_count_s=0;
        s.source_identity=identity_;s.channel_provenance="controlled abstract uniform proper neutron sink proportional to n_n";
        s.revision_identity="phase6a1-controlled-neutron-sink-v1";s.domain_identity=domain_;s.star_identity=star_;s.sequence_state_identity=sequence_;
        s.events.push_back({"controlled-neutron-event","abstract-neutron-disappearance",fate_,partition_,{-1,0,0},-Bdot_});
        s.Validate();return s;
    }
    void RequireCurrent() const override
    {if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale controlled neutron history");}
    const std::string& Identity() const override{RequireCurrent();return identity_;}
  private:
    double B0_,Bdot_;std::string domain_,star_,sequence_,partition_,fate_;
    std::shared_ptr<const BNV::BnvDependencyToken> token_;std::uint64_t generation_;
    const std::string identity_="Phase-6A-1 controlled abstract uniform proper neutron sink v1";
};

// The controlled campaign specialization: a mathematical proper neutron sink
// proportional to n_n.  gamma is only the normalization of this test measure.
class UniformProperNeutronSinkHistory final:public BNV::OrdinaryMatterBnvHistory
{
  public:
    UniformProperNeutronSinkHistory(std::shared_ptr<const CompactStar::Core::NStar> star,double B0,double Bdot,
      std::string domain,std::string star_identity,std::string sequence,std::string partition,std::string fate,
      std::shared_ptr<const BNV::BnvDependencyToken> token)
      :star_(std::move(star)),B0_(B0),Bdot_(Bdot),domain_(std::move(domain)),star_identity_(std::move(star_identity)),
       sequence_(std::move(sequence)),partition_(std::move(partition)),fate_(std::move(fate)),token_(std::move(token)),
       generation_(token_?token_->generation:0),profile_version_(star_?star_->Profile().Version():0)
    {
        if(!star_||!(B0_>0)||!(Bdot_<0)||!std::isfinite(Bdot_))throw std::runtime_error("invalid uniform proper neutron sink");
        const auto& p=star_->Profile();auto R=p.GetRadius(),M=p.GetMass(),Nu=p.GetMetricNu(),Nb=p.GetBaryonDensity(),Yn=p.GetSpeciesPtr("10");
        auto rule=gsl_integration_glfixed_table_alloc(32);if(!rule)throw std::bad_alloc();
        std::vector<double> raw;double total=0;
        for(std::size_t i=0;i<R->Size();++i)
        {
            const double a=i?(*R)[i-1]:0,b=(*R)[i];long double wsum=0,nsum=0,rsum=0,phisum=0;
            for(int j=0;j<32;++j){double r,w;gsl_integration_glfixed_point(a,b,j,&r,&w,rule);const double u=(r-a)/(b-a);
              const double mass=i?(*M)[i-1]+u*((*M)[i]-(*M)[i-1]):(*M)[0]*u*u*u;
              const double nu=i?(*Nu)[i-1]+u*((*Nu)[i]-(*Nu)[i-1]):(*Nu)[0];
              const double nb=i?(*Nb)[i-1]+u*((*Nb)[i]-(*Nb)[i-1]):(*Nb)[0];
              const double yn=i?(*Yn)[i-1]+u*((*Yn)[i]-(*Yn)[i-1]):(*Yn)[0];
              const double measure=w*4*M_PI*r*r/std::sqrt(1-2*mass/r)*std::exp(nu)*nb*yn*1e54;
              wsum+=measure;nsum+=measure*nb*yn;rsum+=measure*r;phisum+=measure*std::exp(nu);}
            const double q=double(wsum);if(q>0){raw.push_back(q);radii_.push_back(double(rsum/wsum));exp_phi_.push_back(double(phisum/wsum));density_.push_back(double(nsum/wsum));total+=q;}
        }
        gsl_integration_glfixed_table_free(rule);if(!(total>0)||raw.empty())throw std::runtime_error("empty uniform proper neutron measure");
        gamma_s_=std::abs(Bdot_)/total;
        for(std::size_t i=0;i<raw.size();++i){mean_radius_km_+=raw[i]*radii_[i]/total;mean_exp_phi_+=raw[i]*exp_phi_[i]/total;mean_density_fm3_+=raw[i]*density_[i]/total;}
        radii_.clear();exp_phi_.clear();density_.clear();
        RequireCurrent();
    }
    BNV::OrdinaryMatterBnvHistorySample Sample(double t)const override
    {
        RequireCurrent();BNV::OrdinaryMatterBnvHistorySample s;s.epoch_s=t;s.B_count=B0_+Bdot_*t;s.Bdot_count_s=Bdot_;
        s.source_count_s={Bdot_,0,0};s.source_identity=identity_;s.channel_provenance="controlled abstract uniform proper neutron sink proportional to n_n; gamma mathematical normalization only";
        s.revision_identity="phase6a1-uniform-proper-neutron-sink-v1";s.domain_identity=domain_;s.star_identity=star_identity_;s.sequence_state_identity=sequence_;
        s.events.push_back({"controlled-neutron-event","abstract-neutron-disappearance",fate_,partition_,{-1,0,0},std::abs(Bdot_),mean_radius_km_,mean_exp_phi_,mean_density_fm3_});
        s.Validate();return s;
    }
    void RequireCurrent()const override{if(!token_||!token_->alive||token_->generation!=generation_||!star_||star_->Profile().Version()!=profile_version_)throw std::runtime_error("stale uniform proper neutron history");}
    const std::string& Identity()const override{RequireCurrent();return identity_;}
    double GammaPerSecond()const{RequireCurrent();return gamma_s_;}
  private:
    std::shared_ptr<const CompactStar::Core::NStar> star_;double B0_,Bdot_,gamma_s_=0;
    std::string domain_,star_identity_,sequence_,partition_,fate_;std::shared_ptr<const BNV::BnvDependencyToken> token_;
    std::uint64_t generation_,profile_version_;double mean_radius_km_=0,mean_exp_phi_=0,mean_density_fm3_=0;std::vector<double> radii_,exp_phi_,density_;
    const std::string identity_="Phase-6A-1 controlled abstract uniform proper neutron sink v1";
};

inline BNV::OrdinaryMatterBnvHistorySample GenericSample(
    double B,const std::array<double,3>& source,double proton_source,
    const std::string& domain,const std::string& star,const std::string& sequence,
    const std::string& partition="synthetic-partition",const std::string& fate="synthetic-fate")
{
    BNV::OrdinaryMatterBnvHistorySample s;
    s.epoch_s=0;s.B_count=B;s.source_count_s=source;s.Bdot_count_s=source[0]+source[1]+source[2];s.proton_source_count_s=proton_source;
    s.source_identity="generic-source";s.channel_provenance="generic charge-balanced mathematical event";s.revision_identity="generic-v1";
    s.domain_identity=domain;s.star_identity=star;s.sequence_state_identity=sequence;
    s.events.push_back({"generic-event","generic-channel",fate,partition,source,1});s.Validate();return s;
}

} // namespace Phase6A1Test
