#include <CompactStar/Physics/Rotochemical/GlobalUrcaChannelCoefficient.hpp>
#include <CompactStar/Physics/Rotochemical/UrcaImbalanceFunctions.hpp>
#include <gsl/gsl_integration.h>
#include <algorithm>
#include <memory>

namespace CompactStar::Physics::Rotochemical
{
void GlobalUrcaChannelCoefficient::RequireCurrent() const
{
    if(!domain_) throw std::runtime_error("missing chemical domain");
    domain_->RequireCurrent();
    if(domain_->Lifetime()->revision->domain!=domain_identity_)
        throw std::runtime_error("Urca chemical domain changed");
}
GlobalUrcaChannelCoefficient GlobalUrcaChannelCoefficient::Compute(const UrcaIntegrationRequest& r)
{
    if(!r.chemical_domain || !r.metric || r.metric_identity.empty()) throw std::runtime_error("Urca integration provenance missing");
    r.chemical_domain->RequireCurrent();
    if(r.domain_identity!=r.chemical_domain->Lifetime()->revision->domain) throw std::runtime_error("Urca/G_y domain mismatch");
    if(r.order<2 || r.order>128 || r.radial_partition_km.size()<2) throw std::runtime_error("invalid Urca quadrature");
    double last=-1;
    for(double x:r.radial_partition_km) {
        if(!std::isfinite(x)||x<0||x<=last) throw std::runtime_error("radial profile must be innermost-first");
        last=x;
    }
    const auto& gp=r.chemical_domain->Partition();
    if(r.radial_partition_km.front()!=gp.front() || r.radial_partition_km.back()!=gp.back())
        throw std::runtime_error("Urca/G_y radial domain mismatch");
    GlobalUrcaChannelCoefficient result;
    result.domain_=r.chemical_domain;result.domain_identity_=r.domain_identity;result.metric_identity_=r.metric_identity;
    result.partition_=r.radial_partition_km;result.order_=r.order;result.selection_=r.selection;
    std::unique_ptr<gsl_integration_glfixed_table,decltype(&gsl_integration_glfixed_table_free)>
        gl(gsl_integration_glfixed_table_alloc(r.order),gsl_integration_glfixed_table_free);
    if(!gl) throw std::runtime_error("Urca quadrature allocation failed");
    std::array<bool,4> seen{};
    for(const auto& local:r.normalizations) {
        auto p=local.Process();auto index=ProcessIndex(p);
        if(seen[index]) throw std::runtime_error("duplicate Urca normalization authority");
        seen[index]=true;
        // Disabled processes remain empty even if positive source and kinematics exist.
        if(!r.selection.Enabled(p)) continue;
        auto& out=result.entries_[index];out.support=local.Support();out.normalization_identity=local.Identity();
        long double sum=0;
        for(auto support:local.Support()) {
            if(support.left_km<gp.front() || support.right_km>gp.back()) throw std::runtime_error("Urca support outside chemical domain");
            std::vector<double> knots{support.left_km};
            for(double x:r.radial_partition_km) if(x>support.left_km&&x<support.right_km) knots.push_back(x);
            knots.push_back(support.right_km);
            for(std::size_t cell=1;cell<knots.size();++cell) for(unsigned j=0;j<r.order;++j) {
                double x,weight;gsl_integration_glfixed_point(knots[cell-1],knots[cell],j,&x,&weight,gl.get());
                const auto metric=r.metric(x);
                if(!std::isfinite(metric.nu)||!std::isfinite(metric.lambda)||metric.lambda<0) throw std::runtime_error("invalid Urca metric");
                // FR2005 (43). Proper volume once; lapse 2-q, not G_y's -nu.
                double term=4*UrcaImbalanceFunctions::Pi()*x*x*std::exp(metric.lambda+(2-static_cast<int>(TemperatureExponent(p)))*metric.nu)*local.LocalS(x)*weight*Units::KM3_TO_CM3;
                if(!(term>=0)||!std::isfinite(term)) throw std::runtime_error("nonfinite Urca integrand");
                sum+=term;
            }
        }
        out.luminosity_erg_s_Kq=static_cast<double>(sum);
    }
    for(auto p:UrcaProcesses) if(r.selection.Enabled(p)&&!seen[ProcessIndex(p)]) throw std::runtime_error("enabled Urca process has no normalization owner");
    result.RequireCurrent();return result;
}
}
