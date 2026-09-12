#pragma once
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
namespace CompactStar::Physics::Rotochemical
{
struct SpinHistorySample {double omega_rad_s,omega_dot_rad_s2;};
// Shared revocation token: changing it invalidates every consumer before sampling.
struct RunDependencyToken {std::uint64_t generation=0;bool alive=true;};
class PrescribedSpinHistory
{
  public:
    virtual ~PrescribedSpinHistory()=default;
    virtual SpinHistorySample Sample(double t)const=0;
    virtual void RequireCurrent()const=0;
    virtual const std::string& Identity()const=0;
};
// External timing-law adapter. The chemical evaluator consumes its single sample.
class PrescribedDipoleHistory final:public PrescribedSpinHistory
{
  public:
    explicit PrescribedDipoleHistory(std::shared_ptr<const RunDependencyToken> token):token_(std::move(token)),generation_(token_?token_->generation:0){RequireCurrent();}
    void RequireCurrent()const override {if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale spin owner");}
    SpinHistorySample Sample(double t)const override {RequireCurrent();if(!(t>=0)||!std::isfinite(t))throw std::runtime_error("invalid spin epoch");const double pp=std::pow(1e8/3.2e19,2),p=std::sqrt(1e-6+2*pp*t),two_pi=2*std::acos(-1.);return {two_pi/p,-two_pi*pp/(p*p*p)};}
    const std::string& Identity()const override {RequireCurrent();return identity_;}
  private:
    const std::shared_ptr<const RunDependencyToken> token_;
    const std::uint64_t generation_;
    const std::string identity_="prescribed dipole B=1e8 G P0=1e-3 s PPdot=(B/3.2e19)^2";
};
}
