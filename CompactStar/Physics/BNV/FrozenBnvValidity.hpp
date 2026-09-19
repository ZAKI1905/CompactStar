#pragma once

#include <CompactStar/Physics/BNV/OrdinaryMatterSource.hpp>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace CompactStar::Physics::BNV
{

struct FrozenValiditySample
{
    double fractional_depletion=0;
    double B_solved_count=0,B_target_residual_count=0,final_bracket_width_count=0;
    std::map<std::string,double> drift_bound,threshold,utilization;
};

class FrozenSensitivityCertificate final
{
  public:
    FrozenSensitivityCertificate(double B0_count,std::string star_identity,
      std::string domain_identity,std::string provenance,
      std::shared_ptr<const BnvDependencyToken>,std::vector<FrozenValiditySample>);
    void RequireCurrent() const;
    double B0Count() const {RequireCurrent();return B0_count_;}
    const std::string& StarIdentity() const {RequireCurrent();return star_identity_;}
    const std::string& DomainIdentity() const {RequireCurrent();return domain_identity_;}
    const std::string& Provenance() const {RequireCurrent();return provenance_;}
    const std::vector<FrozenValiditySample>& Samples() const {RequireCurrent();return samples_;}
  private:
    double B0_count_=0;
    std::string star_identity_,domain_identity_,provenance_;
    std::shared_ptr<const BnvDependencyToken> token_;
    std::uint64_t generation_=0;
    std::vector<FrozenValiditySample> samples_;
};

struct FrozenValidityResult
{
    double fractional_depletion=0,max_utilization=0;
    std::string limiting_quantity;
    std::map<std::string,double> drift_bound,threshold,utilization;
    bool valid=false;
};

class FrozenBnvValidityMonitor final
{
  public:
    explicit FrozenBnvValidityMonitor(std::shared_ptr<const FrozenSensitivityCertificate>);
    FrozenValidityResult Evaluate(double B_count) const;
    FrozenValidityResult RequireValid(double B_count) const;
    const std::shared_ptr<const FrozenSensitivityCertificate>& Certificate() const{return certificate_;}
  private:
    std::shared_ptr<const FrozenSensitivityCertificate> certificate_;
};

} // namespace CompactStar::Physics::BNV
