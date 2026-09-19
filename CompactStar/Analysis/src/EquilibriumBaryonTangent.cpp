#include <CompactStar/Analysis/EquilibriumBaryonTangent.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace CompactStar::Analysis
{
namespace
{
std::size_t Find(const NumberMetadata& m, const std::string& name)
{
    for (std::size_t i=0;i<m.species.size();++i)
        if (m.species[i].label==name) return i;
    throw std::runtime_error("equilibrium tangent missing species " + name);
}

std::string SerializeDomain(const NumberDomain& d)
{
    std::ostringstream s;
    s.precision(17);
    s << static_cast<int>(d.type) << ':' << d.inner_pressure_km_minus2 << ':'
      << d.outer_pressure_km_minus2 << ':' << d.boundary_definition;
    return s.str();
}

std::string SequenceIdentity(const NumberMetadata& m)
{
    std::ostringstream s;
    s.precision(17);
    s << m.central_state_definition << ':' << m.central_energy_km_minus2 << ':'
      << m.sequence_branch_policy << ':' << m.sequence_radial_resolution << ':'
      << m.differentiation_coordinate << ':' << m.tail_callback_identity << ':'
      << m.tail_callback_revision;
    if (!m.sources.empty() && m.sources.front().eos)
        s << ':' << m.sources.front().eos->identity << ':'
          << m.sources.front().eos->revision << ':'
          << m.sources.front().eos->table_path;
    return s.str();
}
}

EquilibriumBaryonTangent EquilibriumBaryonTangent::Compute(
    std::shared_ptr<const EquilibriumSequenceNumberDerivative> sequence,
    double B0_count, std::string star_identity, std::string source_domain_identity)
{
    if (!sequence) throw std::runtime_error("equilibrium tangent requires Phase-5B sequence derivative");
    sequence->RequireCurrent();
    if (!(B0_count>0) || !std::isfinite(B0_count)) throw std::runtime_error("invalid tangent B0");
    if (star_identity.empty() || source_domain_identity.empty()) throw std::runtime_error("missing tangent identity");
    const auto& m=sequence->metadata;
    if (m.domain.type!=DomainType::WholeStar) throw std::runtime_error("tangent requires whole-star domain");
    if (SerializeDomain(m.domain)!=source_domain_identity) throw std::runtime_error("tangent/source domain mismatch");
    const auto& v=sequence->Values();
    const auto& e=sequence->Errors();
    if (v.size()!=m.species.size() || e.size()!=v.size()) throw std::runtime_error("malformed Phase-5B derivative");
    for (double x:v) if (!std::isfinite(x)) throw std::runtime_error("nonfinite Phase-5B derivative");
    for (double x:e) if (!(x>=0) || !std::isfinite(x)) throw std::runtime_error("invalid Phase-5B derivative error");

    // These are the exact authenticated Phase-5B StarProfile axes.  Their
    // particle semantics are additionally checked below; aliases are refused.
    const auto n=Find(m,"10"), p=Find(m,"11"), el=Find(m,"0"), mu=Find(m,"1");
    if (m.species[n].baryon_number!=1 || m.species[p].baryon_number!=1 ||
        m.species[el].baryon_number!=0 || m.species[mu].baryon_number!=0 ||
        m.species[n].charge!=0 || m.species[p].charge!=1 ||
        m.species[el].charge!=-1 || m.species[mu].charge!=-1)
        throw std::runtime_error("unexpected Phase-5B species semantics");

    EquilibriumBaryonTangent out;
    out.sequence_=std::move(sequence);
    out.B0_count_=B0_count;
    out.star_identity_=std::move(star_identity);
    out.domain_identity_=std::move(source_domain_identity);
    out.sequence_identity_=SequenceIdentity(m);
    out.B_B_=v[n]+v[p];
    out.B_B_error_=e[n]+e[p];
    if (!(out.B_B_>out.B_B_error_) || !std::isfinite(out.B_B_))
        throw std::runtime_error("ill-conditioned canonical baryon derivative");
    const double reduced=v[n]+v[el]+v[mu];
    const double reduced_budget=out.B_B_error_+e[n]+e[el]+e[mu]
      +32*std::numeric_limits<double>::epsilon()*(std::abs(out.B_B_)+std::abs(reduced));
    if (std::abs(reduced-out.B_B_)>reduced_budget)
        throw std::runtime_error("Phase-5B reduced baryon closure failed");

    const std::array<std::size_t,3> index{n,el,mu};
    for (std::size_t i=0;i<3;++i)
    {
        const auto k=index[i];
        const double t=v[k]/out.B_B_;
        const double error=(e[k]+std::abs(t)*out.B_B_error_)/(std::abs(out.B_B_)-out.B_B_error_);
        if (!std::isfinite(t) || !(error>=0) || !std::isfinite(error))
            throw std::runtime_error("invalid equilibrium tangent component");
        out.components_[i]={t,t,error};
    }
    out.raw_closure_residual_=out.components_[0].raw+out.components_[1].raw+out.components_[2].raw-1;
    out.closure_budget_=out.components_[0].numerical_error+out.components_[1].numerical_error+
      out.components_[2].numerical_error+32*std::numeric_limits<double>::epsilon()*
      (std::abs(out.components_[0].raw)+std::abs(out.components_[1].raw)+std::abs(out.components_[2].raw));
    if (std::abs(out.raw_closure_residual_)>out.closure_budget_)
        throw std::runtime_error("equilibrium tangent closure outside certification");
    out.components_[0].closed=1-out.components_[1].closed-out.components_[2].closed;
    if (std::abs(out.components_[0].closed-out.components_[0].raw)>out.components_[0].numerical_error)
        throw std::runtime_error("tangent closure adjustment exceeds numerical budget");
    out.RequireCurrent();
    return out;
}

std::string EquilibriumBaryonTangent::SerializeDomainIdentity(const NumberDomain& d)
{ return SerializeDomain(d); }

void EquilibriumBaryonTangent::RequireCurrent() const
{
    if (!sequence_) throw std::runtime_error("missing equilibrium tangent source");
    sequence_->RequireCurrent();
    if (SerializeDomain(sequence_->metadata.domain)!=domain_identity_ ||
        SequenceIdentity(sequence_->metadata)!=sequence_identity_)
        throw std::runtime_error("stale equilibrium tangent identity");
}

const TangentComponent& EquilibriumBaryonTangent::Component(OrdinaryMatterAxis a) const
{
    RequireCurrent();
    switch(a) {
      case OrdinaryMatterAxis::Neutron: return components_[0];
      case OrdinaryMatterAxis::Electron: return components_[1];
      case OrdinaryMatterAxis::Muon: return components_[2];
    }
    throw std::runtime_error("invalid ordinary-matter tangent axis");
}

std::array<double,3> EquilibriumBaryonTangent::ClosedValues() const
{ RequireCurrent(); return {components_[0].closed,components_[1].closed,components_[2].closed}; }
std::array<double,3> EquilibriumBaryonTangent::RawValues() const
{ RequireCurrent(); return {components_[0].raw,components_[1].raw,components_[2].raw}; }
std::array<double,3> EquilibriumBaryonTangent::NumericalErrors() const
{ RequireCurrent(); return {components_[0].numerical_error,components_[1].numerical_error,components_[2].numerical_error}; }

ValidatedTangentSnapshot EquilibriumBaryonTangent::Snapshot() const
{
    RequireCurrent();ValidatedTangentSnapshot out;
    out.closed={components_[0].closed,components_[1].closed,components_[2].closed};
    out.numerical_error={components_[0].numerical_error,components_[1].numerical_error,components_[2].numerical_error};
    out.B0_count=B0_count_;out.closure_budget=closure_budget_;out.star_identity=star_identity_;
    out.domain_identity=domain_identity_;out.sequence_state_identity=sequence_identity_;return out;
}

} // namespace CompactStar::Analysis
