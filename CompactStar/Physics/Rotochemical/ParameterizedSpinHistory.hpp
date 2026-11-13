#pragma once
#include <CompactStar/Physics/Rotochemical/PrescribedSpinHistory.hpp>
#include <iomanip>
#include <sstream>

namespace CompactStar::Physics::Rotochemical
{
// ADR-0025 H-new/S-obs: imposed timing, no separate BNV torque. Candidate.
// P(t)^2=P0^2+2 P0 Pdot0 t. The run context checks the star-dependent limit.
class ParameterizedDipoleHistory final : public PrescribedSpinHistory
{
  public:
    ParameterizedDipoleHistory(double period_s, double period_dot,
                              std::shared_ptr<const RunDependencyToken> token)
      :period_s_(period_s),period_dot_(period_dot),token_(std::move(token)),
       generation_(token_?token_->generation:0)
    {
        if(!(period_s_>0)||!(period_dot_>=0)||!std::isfinite(period_s_)||!std::isfinite(period_dot_))
            throw std::runtime_error("invalid prescribed period or period derivative");
        std::ostringstream s;s<<std::setprecision(17)
          <<"ADR-0025 prescribed PPdot; P0_s="<<period_s_<<" Pdot0="<<period_dot_
          <<"; common timing S-obs; no additional BNV torque";
        identity_=s.str();RequireCurrent();
    }
    void RequireCurrent() const override
    {
        if(!token_||!token_->alive||token_->generation!=generation_)
            throw std::runtime_error("stale parameterized spin owner");
    }
    SpinHistorySample Sample(double t) const override
    {
        RequireCurrent();
        if(!(t>=0)||!std::isfinite(t))throw std::runtime_error("invalid prescribed spin epoch");
        const double pp=period_s_*period_dot_;
        const double period=std::sqrt(period_s_*period_s_+2*pp*t);
        if(!(period>0)||!std::isfinite(period))throw std::runtime_error("nonfinite prescribed spin period");
        const double two_pi=2*std::acos(-1.);
        return {two_pi/period,-two_pi*pp/(period*period*period)};
    }
    const std::string& Identity() const override {RequireCurrent();return identity_;}
  private:
    const double period_s_,period_dot_;
    const std::shared_ptr<const RunDependencyToken> token_;
    const std::uint64_t generation_;
    std::string identity_;
};
}
