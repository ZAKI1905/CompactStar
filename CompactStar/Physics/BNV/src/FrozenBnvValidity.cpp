#include <CompactStar/Physics/BNV/FrozenBnvValidity.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::BNV
{
namespace
{
const std::vector<std::string> Required{
 "t_n","t_e","t_mu","Z_row_npe","Z_row_npmu","Cstar","Ltilde_Me","Ltilde_Mmu",
 "mu_B","mu_n_profile","P0_average","P1_average","P2_average","N_n","N_e","N_mu",
 "species_support","metric_structure","radius","surface_gravity","envelope","Tsurface_inf"};
}

FrozenSensitivityCertificate::FrozenSensitivityCertificate(double B0,std::string star,std::string domain,
    std::string provenance,std::shared_ptr<const BnvDependencyToken> token,std::vector<FrozenValiditySample> samples)
  :B0_count_(B0),star_identity_(std::move(star)),domain_identity_(std::move(domain)),provenance_(std::move(provenance)),
   token_(std::move(token)),generation_(token_?token_->generation:0),samples_(std::move(samples))
{
    if(!(B0_count_>0)||!std::isfinite(B0_count_)||star_identity_.empty()||domain_identity_.empty()||provenance_.empty()||!token_)
        throw std::runtime_error("invalid frozen certificate identity");
    if(samples_.size()!=21)throw std::runtime_error("frozen certificate requires 21 stars");
    const double tau=5e-11*B0_count_;
    std::map<std::string,double> previous;
    for(std::size_t i=0;i<samples_.size();++i)
    {
        const double expected=-5e-8*static_cast<double>(i);
        auto& s=samples_[i];
        if(std::abs(s.fractional_depletion-expected)>32*std::numeric_limits<double>::epsilon())
            throw std::runtime_error("frozen certificate grid changed");
        if(std::abs(s.B_target_residual_count)>tau||s.final_bracket_width_count>tau||s.final_bracket_width_count<0)
            throw std::runtime_error("frozen certificate target tolerance failed");
        const double target=B0_count_*(1+expected);
        const double consistency=64*std::numeric_limits<double>::epsilon()*B0_count_;
        if(!(s.B_solved_count>0)||!std::isfinite(s.B_solved_count)||
           std::abs((s.B_solved_count-target)-s.B_target_residual_count)>consistency)
            throw std::runtime_error("frozen certificate achieved-B identity failed");
        for(const auto& name:Required){const auto it=s.utilization.find(name),d=s.drift_bound.find(name),t=s.threshold.find(name);
          if(it==s.utilization.end()||d==s.drift_bound.end()||t==s.threshold.end()||!(d->second>=0)||!(t->second>0)||
             !std::isfinite(d->second)||!std::isfinite(t->second)||!(it->second>=0)||!std::isfinite(it->second)||it->second>1||
             std::abs(it->second-d->second/t->second)>64*std::numeric_limits<double>::epsilon()*std::max(1.0,it->second)||
             (i&&it->second<previous.at(name)))
              throw std::runtime_error("frozen certificate quantity failed: "+name);
          previous[name]=it->second;}
        if(i&&!(samples_[i].B_solved_count<samples_[i-1].B_solved_count))
            throw std::runtime_error("frozen certificate achieved-B abscissae not monotone");
    }
    RequireCurrent();
}
void FrozenSensitivityCertificate::RequireCurrent() const
{if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale frozen certificate");}

FrozenBnvValidityMonitor::FrozenBnvValidityMonitor(std::shared_ptr<const FrozenSensitivityCertificate> c):certificate_(std::move(c))
{if(!certificate_)throw std::runtime_error("missing frozen certificate");certificate_->RequireCurrent();}

FrozenValidityResult FrozenBnvValidityMonitor::Evaluate(double B) const
{
    certificate_->RequireCurrent();
    if(!(B>0)||!std::isfinite(B))throw std::runtime_error("invalid runtime baryon count");
    FrozenValidityResult out;out.fractional_depletion=(B-certificate_->B0Count())/certificate_->B0Count();
    if(out.fractional_depletion>0||out.fractional_depletion < -1e-6){out.max_utilization=std::numeric_limits<double>::infinity();out.limiting_quantity="DeltaB_over_B0";return out;}
    const auto& v=certificate_->Samples();
    const double depletion=-out.fractional_depletion;
    std::size_t hi=0;
    while(hi<v.size() && -(v[hi].B_solved_count-certificate_->B0Count())/certificate_->B0Count()<depletion)++hi;
    hi=std::min(hi,v.size()-1);const std::size_t lo=hi?hi-1:0;
    const double dlo=-(v[lo].B_solved_count-certificate_->B0Count())/certificate_->B0Count();
    const double dhi=-(v[hi].B_solved_count-certificate_->B0Count())/certificate_->B0Count();
    const double a=hi==lo?0:std::clamp((depletion-dlo)/(dhi-dlo),0.0,1.0);
    for(const auto& name:Required)
    {
        const double x=(1-a)*v[lo].utilization.at(name)+a*v[hi].utilization.at(name);
        out.drift_bound[name]=(1-a)*v[lo].drift_bound.at(name)+a*v[hi].drift_bound.at(name);
        out.threshold[name]=(1-a)*v[lo].threshold.at(name)+a*v[hi].threshold.at(name);
        out.utilization[name]=x;
        if(x>=out.max_utilization){out.max_utilization=x;out.limiting_quantity=name;}
    }
    out.valid=out.max_utilization<=1;
    return out;
}
FrozenValidityResult FrozenBnvValidityMonitor::RequireValid(double B) const
{auto r=Evaluate(B);if(!r.valid)throw std::runtime_error("frozen BNV validity failed: "+r.limiting_quantity);return r;}

} // namespace CompactStar::Physics::BNV
