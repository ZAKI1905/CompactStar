// ADR-0011 structural layer. Scientific candidate; no chemical susceptibility or evolution.
#pragma once
#include <CompactStar/Core/NStar.hpp>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace CompactStar::Analysis
{
struct Species
{
    std::string label; // exact StarProfile column label; ordering is never inferred
    double baryon_number = 0, charge = 0;
};
enum class DomainType { WholeStar, ExplicitFixedIsobar };
struct NumberDomain
{
    DomainType type = DomainType::WholeStar;
    double outer_pressure_km_minus2 = 0;
    // Zero means the regular centre. A positive value declares a moving inner isobar.
    double inner_pressure_km_minus2 = 0;
    std::string boundary_definition;
    bool operator==(const NumberDomain &) const;
};
// This descriptor is owned by the EOS adapter. A change of revision/domain invalidates results.
// If a table is used, every result also retains and checks its exact file contents.
struct NumberEosSource
{
    std::string identity, revision, physical_domain, table_path;
};
enum class NumberSurface { ExactVacuum, FiniteCutWithBoundedRemainder };
struct NumberTail
{
    NumberSurface semantics = NumberSurface::ExactVacuum;
    std::string authority;
    // Net correction from the cut-isobar quantity (including its terminal atom) to P=0.
    // A zero estimate with a nonzero proven bound is permitted; R_cut is never called R_0.
    std::vector<double> count_correction, count_error;
    std::vector<double> response_correction, response_error; // count km^2
    std::string policy_revision;
};
struct NumberInput
{
    const Core::NStar *star = nullptr; // non-owning: results must not outlive this star
    std::shared_ptr<const NumberEosSource> eos;
    std::vector<Species> species;
    NumberDomain domain;
    NumberTail tail;
    std::string central_state_definition;
};
struct NumberProvenance
{
    const Core::StarProfile *profile = nullptr;
    std::uint64_t profile_version = 0;
    const Core::NStar *star = nullptr;
    const Core::HartleFirstOrderResponse *first_order = nullptr;
    const Core::HartleMonopoleResponse *monopole = nullptr;
    std::shared_ptr<const NumberEosSource> eos;
    NumberEosSource eos_snapshot;
    std::string table_contents;
    void RequireCurrent() const;
};
struct NumberMetadata
{
    std::vector<Species> species;
    NumberDomain domain;
    NumberSurface surface;
    std::string surface_authority, central_state_definition;
    // Eager snapshots expose the actual tail inputs and realized policy for every
    // contributing star, including a sequence callback's returned corrections/bounds.
    std::vector<NumberTail> tail_inputs;
    std::vector<std::size_t> profile_node_counts;
    std::size_t sequence_radial_resolution = 0;
    std::string sequence_branch_policy, tail_callback_identity, tail_callback_revision;
    std::string refinement_policy = "original profile nodes; no immutable refinement knots";
    std::string q_normalization = "q=Omega_geom^2; Omega_geom=Omega_phys/c; q in km^-2";
    std::string profile_units = "r,m:km; epsilon,p:km^-2; n_i=Y_i*n_B:fm^-3; nu:dimensionless";
    std::string differentiation_coordinate;
    double central_energy_km_minus2 = 0;
    std::vector<double> steps, achieved_central_energy_km_minus2;
    std::vector<NumberProvenance> sources;
};
// Eager values. All supported numerical access validates EVERY source, including sequence stars.
class NumberResult
{
  public:
    bool valid = false;
    std::string refusal_reason;
    NumberMetadata metadata;
    const std::vector<double> &Values() const;
    const std::vector<double> &Errors() const;
    void RequireCurrent() const;
  protected:
    std::vector<double> values_, errors_;
};
class ParticleNumbers : public NumberResult
{
  public:
    static ParticleNumbers Compute(const NumberInput &); // count
};
class FixedCentralEnergyNumberResponse : public NumberResult
{
  public:
    static FixedCentralEnergyNumberResponse Compute(const NumberInput &); // count km^2
    // Density-measure, radial metric, velocity, outer-minus-inner boundary terms, before tail.
    std::vector<std::vector<double>> contributions;
};
struct NumberSequenceRecipe
{
    std::shared_ptr<const NumberEosSource> eos;
    std::vector<Species> species;
    NumberDomain domain;
    double rho_reference_g_cm3 = 0;
    // One smooth physical sequence branch, authenticated by the EOS/validation adapter.
    double minimum_rho_g_cm3 = 0, maximum_rho_g_cm3 = 0;
    std::string branch_definition;
    std::vector<double> log_steps{0.004, 0.002, 0.001};
    std::size_t radial_resolution = 40000;
    std::string work_directory;
    // Computes this newly solved star's own surface remainder, never the central star's tail.
    std::function<NumberTail(const Core::NStar &)> tail_for_star;
    // Callback must be pure for its owning star and immutable declared policy.
    // Results retain its identity/revision and every realized value, not a lazy callback.
    std::string tail_policy_identity, tail_policy_revision;
};
class EquilibriumSequenceNumberDerivative : public NumberResult
{
  public:
    static EquilibriumSequenceNumberDerivative Compute(const NumberSequenceRecipe &);
    std::vector<std::vector<double>> three_point, five_point; // dN/d epsilon_c, count km^2
    // Complete independently solved stars; shared ownership is solely lifetime management.
    std::vector<std::shared_ptr<Core::NStar>> contributing_stars;
};
class FixedBaryonNumberResponse : public NumberResult
{
  public:
    static FixedBaryonNumberResponse Compute(const FixedCentralEnergyNumberResponse &,
                                            const EquilibriumSequenceNumberDerivative &);
    double A_B = 0, B_B = 0, B_B_error = 0, conditioning = 0;
    double central_energy_per_q = 0, central_energy_per_q_error = 0;
    double baryon_residual = 0, baryon_budget = 0, charge_residual = 0, charge_budget = 0;
    // PN7/PN8: for a core/shell, reduce using the WHOLE STAR constraint supplied here.
    static FixedBaryonNumberResponse ReduceDomain(const FixedCentralEnergyNumberResponse &,
          const EquilibriumSequenceNumberDerivative &, const FixedBaryonNumberResponse &whole);
    std::vector<double> WholeStarIPhysical() const; // count s^2, conditional P=0 mapping
    std::vector<double> WholeStarEquilibriumNumberRate(double omega_rad_per_second,
                                                     double omega_dot_rad_per_second2) const; // count/s
    // Positive outward flux Q_B at each moving isobar; shell inner term has opposite sign.
    std::vector<double> FixedIsobarIGeometric(const std::vector<double> &outer_Y,
          double outer_Q_B, const std::vector<double> &inner_Y = {}, double inner_Q_B = 0) const;
};

// Generic numerical representation, also independently testable with analytic fixtures.
// A zero-length segment is allowed ONLY when explicitly declared as a true jump.
// Continuous onsets are ordinary segments. There is no gradient/steepness classification.
struct NumberMeasureNode
{
    double r = 0, m = 0, nu = 0, density = 0, xi_hat = 0, m_hat = 0, s = 0;
};
struct NumberMeasureSegment
{
    NumberMeasureNode left, right;
    bool declared_jump = false;
};
struct NumberMeasureIntegral
{
    double count = 0, count_error = 0;
    double density_measure = 0, metric = 0, velocity = 0, boundary = 0, response_error = 0;
};
// Units returned already include fm^-3 km^3 -> count (1e54).
NumberMeasureIntegral IntegrateNumberMeasure(const std::vector<NumberMeasureSegment> &,
                                             bool moving_outer, bool moving_inner);
} // namespace CompactStar::Analysis
