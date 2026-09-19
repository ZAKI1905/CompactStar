#include <CompactStar/Physics/BNV/MovingReferenceSource.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CompactStar::Physics::BNV
{

MovingReferenceSample MovingReferenceSource::Project(
    const OrdinaryMatterBnvHistorySample& source,
    const Analysis::EquilibriumBaryonTangent& tangent)
{
    return Project(source,tangent.Snapshot());
}

MovingReferenceSample MovingReferenceSource::Project(
    const OrdinaryMatterBnvHistorySample& source,
    const Analysis::ValidatedTangentSnapshot& tangent)
{
    source.Validate();
    if(source.domain_identity!=tangent.domain_identity) throw std::runtime_error("source/t domain mismatch");
    if(source.star_identity!=tangent.star_identity) throw std::runtime_error("source/t star mismatch");
    if(source.sequence_state_identity!=tangent.sequence_state_identity) throw std::runtime_error("source/t sequence mismatch");
    const auto& t=tangent.closed;
    const auto& te=tangent.numerical_error;
    MovingReferenceSample out;
    out.source=source;
    double max_lift=0;
    for(std::size_t i=0;i<3;++i)
    {
        out.Sigma_count_s[i]=source.source_count_s[i]-t[i]*source.Bdot_count_s;
        const double tau=64*std::numeric_limits<double>::epsilon()*
          std::max({1.0,std::abs(source.source_count_s[i]),std::abs(t[i]*source.Bdot_count_s),std::abs(out.Sigma_count_s[i])})+
          std::abs(source.Bdot_count_s)*te[i];
        max_lift=std::max(max_lift,tau);
    }
    out.sigma_count_s={out.Sigma_count_s[1],out.Sigma_count_s[2]};
    const std::array<double,3> lift{-out.sigma_count_s[0]-out.sigma_count_s[1],out.sigma_count_s[0],out.sigma_count_s[1]};
    for(std::size_t i=0;i<3;++i)out.lift_residual_count_s=std::max(out.lift_residual_count_s,std::abs(out.Sigma_count_s[i]-lift[i]));
    out.baryon_residual_count_s=out.Sigma_count_s[0]+out.Sigma_count_s[1]+out.Sigma_count_s[2];
    out.baryon_budget_count_s=source.SourceTolerance()+std::abs(source.Bdot_count_s)*tangent.closure_budget+
      64*std::numeric_limits<double>::epsilon()*std::max({1.0,
      std::abs(source.source_count_s[0])+std::abs(source.source_count_s[1])+std::abs(source.source_count_s[2]),
      std::abs(source.Bdot_count_s)});
    out.lift_budget_count_s=max_lift;
    if(std::abs(out.baryon_residual_count_s)>out.baryon_budget_count_s)
        throw std::runtime_error("moving-reference baryon neutrality failed");
    if(out.lift_residual_count_s>out.lift_budget_count_s)
        throw std::runtime_error("moving-reference exact lift failed");
    return out;
}

} // namespace CompactStar::Physics::BNV
