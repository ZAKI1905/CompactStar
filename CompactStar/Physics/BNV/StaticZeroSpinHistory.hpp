#pragma once

#include <CompactStar/Physics/Rotochemical/PrescribedSpinHistory.hpp>

namespace CompactStar::Physics::BNV
{

class StaticZeroSpinHistory final:public Rotochemical::PrescribedSpinHistory
{
  public:
    explicit StaticZeroSpinHistory(std::shared_ptr<const Rotochemical::RunDependencyToken>);
    Rotochemical::SpinHistorySample Sample(double epoch_s) const override;
    void RequireCurrent() const override;
    const std::string& Identity() const override;
  private:
    const std::shared_ptr<const Rotochemical::RunDependencyToken> token_;
    const std::uint64_t generation_;
    const std::string identity_=
      "static zero spin Omega=0 rad s^-1 OmegaDot=0 rad s^-2; Phase-6A-1 controlled abstract BNV";
};

} // namespace CompactStar::Physics::BNV

