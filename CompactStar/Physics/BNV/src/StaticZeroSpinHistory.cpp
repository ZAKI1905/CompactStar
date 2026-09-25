#include <CompactStar/Physics/BNV/StaticZeroSpinHistory.hpp>
#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::BNV
{
StaticZeroSpinHistory::StaticZeroSpinHistory(std::shared_ptr<const Rotochemical::RunDependencyToken> token)
  :token_(std::move(token)),generation_(token_?token_->generation:0){RequireCurrent();}
void StaticZeroSpinHistory::RequireCurrent() const
{if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale static zero-spin owner");}
Rotochemical::SpinHistorySample StaticZeroSpinHistory::Sample(double t) const
{RequireCurrent();if(!(t>=0)||!std::isfinite(t))throw std::runtime_error("invalid static zero-spin epoch");return {0,0};}
const std::string& StaticZeroSpinHistory::Identity() const {RequireCurrent();return identity_;}
} // namespace CompactStar::Physics::BNV

