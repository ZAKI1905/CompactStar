#pragma once

#include <CompactStar/Analysis/ParticleNumberResponse.hpp>
#include <array>
#include <memory>
#include <string>

namespace CompactStar::Physics::BNV
{

enum class OrdinaryMatterAxis { Neutron, Electron, Muon };

struct TangentComponent
{
    double raw = 0;
    double closed = 0;
    double numerical_error = 0;
};

struct ValidatedTangentSnapshot
{
    std::array<double,3> closed{},numerical_error{};
    double B0_count=0,closure_budget=0;
    std::string star_identity,domain_identity,sequence_state_identity;
};

// Phase-6A-1 typed view of the Phase-5B equilibrium-sequence derivative.
// It never computes a sequence and never consumes G_y.
class EquilibriumBaryonTangent final
{
  public:
    static EquilibriumBaryonTangent Compute(
        std::shared_ptr<const Analysis::EquilibriumSequenceNumberDerivative> sequence,
        double B0_count,
        std::string star_identity,
        std::string source_domain_identity);
    static std::string SerializeDomainIdentity(const Analysis::NumberDomain&);

    void RequireCurrent() const;
    void RequireCheapCurrent() const;
    const TangentComponent& Component(OrdinaryMatterAxis) const;
    std::array<double,3> ClosedValues() const;
    std::array<double,3> RawValues() const;
    std::array<double,3> NumericalErrors() const;
    ValidatedTangentSnapshot Snapshot() const;
    ValidatedTangentSnapshot SnapshotCheap() const;
    double B0Count() const { RequireCurrent(); return B0_count_; }
    double CanonicalBaryonDerivative() const { RequireCurrent(); return B_B_; }
    double CanonicalBaryonDerivativeError() const { RequireCurrent(); return B_B_error_; }
    double RawClosureResidual() const { RequireCurrent(); return raw_closure_residual_; }
    double ClosureBudget() const { RequireCurrent(); return closure_budget_; }
    const std::string& StarIdentity() const { RequireCurrent(); return star_identity_; }
    const std::string& DomainIdentity() const { RequireCurrent(); return domain_identity_; }
    const std::string& SequenceStateIdentity() const { RequireCurrent(); return sequence_identity_; }
    const std::shared_ptr<const Analysis::EquilibriumSequenceNumberDerivative>& Source() const
    { RequireCurrent(); return sequence_; }

  private:
    struct CheapSourceSnapshot
    {
        const Core::StarProfile* profile=nullptr;
        const Core::NStar* star=nullptr;
        const Core::HartleFirstOrderResponse* first_order=nullptr;
        const Core::HartleMonopoleResponse* monopole=nullptr;
        std::shared_ptr<const Analysis::NumberEosSource> eos;
        Analysis::NumberEosSource eos_snapshot;
        std::uint64_t profile_version=0;
    };
    std::shared_ptr<const Analysis::EquilibriumSequenceNumberDerivative> sequence_;
    std::array<TangentComponent,3> components_{};
    double B0_count_ = 0, B_B_ = 0, B_B_error_ = 0;
    double raw_closure_residual_ = 0, closure_budget_ = 0;
    std::string star_identity_, domain_identity_, sequence_identity_;
    std::vector<CheapSourceSnapshot> cheap_sources_;
};

} // namespace CompactStar::Physics::BNV
