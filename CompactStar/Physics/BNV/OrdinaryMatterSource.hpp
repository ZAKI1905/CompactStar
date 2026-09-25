#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace CompactStar::Physics::BNV
{

struct BnvDependencyToken { std::uint64_t generation=0; bool alive=true; };

struct OrdinaryMatterEventMeasure
{
    std::string event_id, channel_id, product_fate_identity, partition_identity;
    std::array<double,3> delta_y{};
    double rate_count_s = 0;
    // Optional immutable quadrature metadata for a spatially resolved event
    // measure.  Global sources may leave these exact zeros except exp_phi.
    double radius_km=0, exp_phi=1, local_neutron_density_fm3=0;
};

struct OrdinaryMatterBnvHistorySample
{
    double epoch_s=0, B_count=0, Bdot_count_s=0;
    std::array<double,3> source_count_s{}; // n,e,mu
    double proton_source_count_s=0;
    std::string source_identity, channel_provenance, revision_identity;
    std::string domain_identity, star_identity, sequence_state_identity;
    std::vector<OrdinaryMatterEventMeasure> events;

    double SourceTolerance() const
    {
        const double sum=std::abs(source_count_s[0])+std::abs(source_count_s[1])+std::abs(source_count_s[2]);
        return 32*std::numeric_limits<double>::epsilon()*std::max({1.0,std::abs(Bdot_count_s),sum});
    }
    void Validate() const
    {
        for(double v:{epoch_s,B_count,Bdot_count_s,source_count_s[0],source_count_s[1],source_count_s[2],proton_source_count_s})
            if(!std::isfinite(v)) throw std::runtime_error("nonfinite ordinary-matter source");
        if (!(epoch_s>=0) || !(B_count>0)) throw std::runtime_error("invalid ordinary-matter history state");
        if (source_identity.empty()||channel_provenance.empty()||revision_identity.empty()||
            domain_identity.empty()||star_identity.empty()||sequence_state_identity.empty())
            throw std::runtime_error("incomplete ordinary-matter source identity");
        const double tau=SourceTolerance();
        const double sum=source_count_s[0]+source_count_s[1]+source_count_s[2];
        if(std::abs(Bdot_count_s-sum)>tau) throw std::runtime_error("Bdot differs from b^T S_y");
        if(std::abs(proton_source_count_s-source_count_s[1]-source_count_s[2])>tau)
            throw std::runtime_error("ordinary charge closure failed");
        std::array<double,3> event_sum{};
        std::set<std::string> ids;
        for(const auto& e:events)
        {
            if(e.event_id.empty()||e.channel_id.empty()||e.product_fate_identity.empty()||e.partition_identity.empty()||
               !ids.insert(e.event_id).second) throw std::runtime_error("malformed or duplicate source event");
            if(!(e.rate_count_s>=0)||!std::isfinite(e.rate_count_s)||!(e.radius_km>=0)||
               !(e.exp_phi>0&&e.exp_phi<=1)||!(e.local_neutron_density_fm3>=0)||
               !std::isfinite(e.radius_km)||!std::isfinite(e.exp_phi)||!std::isfinite(e.local_neutron_density_fm3))
                throw std::runtime_error("invalid event measure metadata");
            for(std::size_t i=0;i<3;++i){if(!std::isfinite(e.delta_y[i]))throw std::runtime_error("nonfinite event stoichiometry");event_sum[i]+=e.delta_y[i]*e.rate_count_s;}
        }
        for(std::size_t i=0;i<3;++i)if(std::abs(event_sum[i]-source_count_s[i])>tau)
            throw std::runtime_error("event measure does not reproduce ordinary source");
    }
};

class OrdinaryMatterBnvHistory
{
  public:
    virtual ~OrdinaryMatterBnvHistory()=default;
    virtual OrdinaryMatterBnvHistorySample Sample(double epoch_s) const=0;
    virtual void RequireCurrent() const=0;
    virtual const std::string& Identity() const=0;
};

} // namespace CompactStar::Physics::BNV

