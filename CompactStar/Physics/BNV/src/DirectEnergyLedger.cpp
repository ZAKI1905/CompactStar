#include <CompactStar/Physics/BNV/DirectEnergyLedger.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace CompactStar::Physics::BNV
{

ActualOrdinaryPotential BnvDirectEnergyLedger::ActualPotential(
    double mu_B,const Rotochemical::ChemicalImbalanceState& eta,
    const EquilibriumBaryonTangent& tangent,std::string provenance)
{
    return ActualPotential(mu_B,eta,tangent.Snapshot(),std::move(provenance));
}

ActualOrdinaryPotential BnvDirectEnergyLedger::ActualPotential(
    double mu_B,const Rotochemical::ChemicalImbalanceState& eta,
    const ValidatedTangentSnapshot& tangent,std::string provenance)
{
    if(!std::isfinite(mu_B)||provenance.empty()) throw std::runtime_error("invalid actual-potential authority");
    const auto& t=tangent.closed;
    const double ee=eta.InfinityMeV(Rotochemical::BetaChannel::Npe);
    const double em=eta.InfinityMeV(Rotochemical::BetaChannel::NpMu);
    const double scalar=t[1]*ee+t[2]*em;
    ActualOrdinaryPotential out;
    out.mu_B_inf_MeV=mu_B;
    out.equilibrium_inf_MeV={mu_B,mu_B,mu_B};
    // g=mu_B b-(I-b t^T)P^T eta, with P selecting the e and mu
    // components of baryon-neutral departures.
    out.actual_inf_MeV={mu_B+scalar,mu_B-ee+scalar,mu_B-em+scalar};
    out.mu_n_actual_inf_MeV=out.actual_inf_MeV[0];
    out.reconstruction_residual_MeV=out.mu_n_actual_inf_MeV-(mu_B+t[1]*ee+t[2]*em);
    out.provenance=std::move(provenance);
    for(double v:out.actual_inf_MeV)if(!std::isfinite(v))throw std::runtime_error("nonfinite actual potential");
    return out;
}

DirectEnergyResult BnvDirectEnergyLedger::Evaluate(
    const OrdinaryMatterBnvHistorySample& source,const ActualOrdinaryPotential& potential,
    const std::vector<DirectEventEnergy>& energies,const ProductFateLedger& fate)
{
    source.Validate();
    if(energies.size()!=source.events.size())throw std::runtime_error("direct-energy event cardinality mismatch");
    std::map<std::string,const DirectEventEnergy*> by_id;
    for(const auto& e:energies)
    {
        if(e.event_id.empty()||!by_id.emplace(e.event_id,&e).second)throw std::runtime_error("duplicate direct-energy event");
        for(double v:{e.Eesc_fluid_inf_MeV,e.Eesc_star_inf_MeV,e.EX_inf_MeV,e.Ein_explicit_inf_MeV,e.exp_phi,e.Qlocal_MeV})
            if(!std::isfinite(v))throw std::runtime_error("nonfinite direct-energy input");
        if(e.Eesc_fluid_inf_MeV<0||e.Eesc_star_inf_MeV<0||e.EX_inf_MeV<0||e.Ein_explicit_inf_MeV<0)
            throw std::runtime_error("negative direct-energy partition");
        if(e.partition_identity.empty()||e.product_fate_identity!=fate.Identity())
            throw std::runtime_error("direct-energy fate identity mismatch");
    }
    DirectEnergyResult out;out.potential=potential;
    for(const auto& event:source.events)
    {
        fate.RequireChannel(event.channel_id);
        const auto it=by_id.find(event.event_id);if(it==by_id.end())throw std::runtime_error("missing direct-energy event");
        const auto& e=*it->second;
        if(e.partition_identity!=event.partition_identity||e.product_fate_identity!=event.product_fate_identity)
            throw std::runtime_error("source/direct partition identity mismatch");
        const double fate_residual=e.Eesc_fluid_inf_MeV-e.Eesc_star_inf_MeV-e.EX_inf_MeV;
        out.escape_fate_residual_MeV=std::max(out.escape_fate_residual_MeV,std::abs(fate_residual));
        const double scale=std::max({1.0,std::abs(e.Eesc_fluid_inf_MeV),std::abs(e.Eesc_star_inf_MeV),std::abs(e.EX_inf_MeV)});
        if(std::abs(fate_residual)>32*std::numeric_limits<double>::epsilon()*scale)
            throw std::runtime_error("fluid/star/product energy closure failed");
        double qeq=e.Ein_explicit_inf_MeV-e.Eesc_fluid_inf_MeV;
        double qactual=qeq;
        for(std::size_t i=0;i<3;++i){qeq-=event.delta_y[i]*potential.equilibrium_inf_MeV[i];qactual-=event.delta_y[i]*potential.actual_inf_MeV[i];}
        if(e.local_frame_check)
        {
            if(!(e.exp_phi>0&&e.exp_phi<=1))throw std::runtime_error("invalid local/infinity redshift");
            out.local_infinity_residual_MeV=std::max(out.local_infinity_residual_MeV,std::abs(qactual-e.exp_phi*e.Qlocal_MeV));
            const double b=64*std::numeric_limits<double>::epsilon()*std::max({1.0,std::abs(qactual),std::abs(e.exp_phi*e.Qlocal_MeV)});
            if(std::abs(qactual-e.exp_phi*e.Qlocal_MeV)>b)throw std::runtime_error("local/infinity direct-energy mismatch");
        }
        out.Qeq_event_inf_MeV.push_back(qeq);out.Qactual_event_inf_MeV.push_back(qactual);
        out.power_eq_MeV_s+=qeq*event.rate_count_s;out.power_actual_MeV_s+=qactual*event.rate_count_s;
    }
    for(double v:{out.power_eq_MeV_s,out.power_actual_MeV_s})if(!std::isfinite(v))throw std::runtime_error("nonfinite direct power");
    // The sole direct-BNV MeV/s -> erg/s conversion boundary.
    out.power_eq_erg_s=Rotochemical::MeVToErg*out.power_eq_MeV_s;
    out.power_actual_erg_s=Rotochemical::MeVToErg*out.power_actual_MeV_s;
    return out;
}

} // namespace CompactStar::Physics::BNV
