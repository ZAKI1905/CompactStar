#pragma once

#include <CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp>
#include <CompactStar/Physics/BNV/OrdinaryMatterSource.hpp>
#include <array>

namespace CompactStar::Physics::BNV
{

struct MovingReferenceSample
{
    OrdinaryMatterBnvHistorySample source;
    std::array<double,3> Sigma_count_s{};
    std::array<double,2> sigma_count_s{};
    double baryon_residual_count_s=0, baryon_budget_count_s=0;
    double lift_residual_count_s=0, lift_budget_count_s=0;
};

class MovingReferenceSource final
{
  public:
    static MovingReferenceSample Project(const OrdinaryMatterBnvHistorySample&,
                                         const EquilibriumBaryonTangent&);
    static MovingReferenceSample Project(const OrdinaryMatterBnvHistorySample&,
                                         const ValidatedTangentSnapshot&);
};

} // namespace CompactStar::Physics::BNV

