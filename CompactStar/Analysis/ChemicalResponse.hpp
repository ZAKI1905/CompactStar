// ADR-0013 generic cold chemical response. AI-authored candidate; no evolution.
#pragma once
#include <CompactStar/Analysis/ParticleNumberResponse.hpp>
#include <CompactStar/EOS/LocalThermodynamics.hpp>
#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace CompactStar::Analysis
{
using ChemicalMatrix = std::vector<std::vector<double>>;
enum class NumberAxis
{
	Neutron,
	Electron,
	Muon
};
enum class ImbalanceChannel
{
	Npe,
	NpMu
};

// A genuinely active Hessian is solved before its response is embedded in y.
// Boundaries have no susceptibility result; no padded-H or projector API exists.
class ChargeNeutralNumberSusceptibility
{
  public:
	static ChargeNeutralNumberSusceptibility Compute(
		const ActiveLocalThermodynamicEvaluation &, const ChemicalMatrix &h_error);
	const ChemicalMatrix &Values() const { return c_; }
	const ChemicalMatrix &NumericalError() const { return error_; }
	const std::vector<NumberAxis> &Support() const { return support_; }
	double Rho() const { return rho_; }
	double Condition() const { return condition_; }

  private:
	ChemicalMatrix c_, error_;
	std::vector<NumberAxis> support_;
	double rho_ = 0, condition_ = 0;
};

struct ChemicalRevision
{
	std::string model, provider_revision, provider_bytes, particle_constants;
	std::string background, metric, basis, domain, partition, onset, tail, accuracy;
	std::string lifetime_token;
	bool alive = true;
	std::string Serialize() const;
};

// Every raw structural source must be covered by these shared owners BEFORE
// RequireCurrent can dereference it. Destruction of caller handles is safe.
struct ChemicalLifetime
{
	std::vector<std::shared_ptr<Core::NStar>> stars;
	std::shared_ptr<const ILocalThermodynamicProvider> provider;
	std::shared_ptr<ChemicalRevision> revision;
	std::vector<std::string> source_files;
};

struct ChemicalMetric
{
	double mass_km, nu;
};
struct ChemicalLocalResponse
{
	ChemicalMatrix c, numerical_error;
	std::vector<NumberAxis> support;
};
struct ChemicalInterval
{
	double left = 0, right = 0;
	bool left_continuous_onset = false, right_continuous_onset = false;
	bool first_order_discontinuity = false;
};
struct ChemicalRefusalCertificate
{
	double left = 0, right = 0, containing_left = 0, containing_right = 0;
	std::size_t first_cell = 0, last_cell = 0;
	std::string branch, availability_authority, geometry_extrema_authority;
	double radius_upper = 0, mass_upper = 0, nu_lower = 0;
	ChemicalMatrix numerical_error;
};
// Source-qualified pe comparison hypotheses are supplied by the provider.
// The integration owner verifies the analytic geometric inequalities itself.
struct ChemicalPeTailCertificate
{
	double cut_radius, cut_mass, cut_nu, h_upper, epsilon_upper;
	double bootstrap_radius_upper, total_mass_upper, radius_upper;
	double susceptibility_upper, lapse_error_upper, response_upper;
	std::string source_hypotheses;
};
struct ChemicalIntegrationRequest
{
	std::shared_ptr<ChemicalLifetime> lifetime;
	std::function<ChemicalMetric(double)> background;
	std::function<ChemicalLocalResponse(double)> local;
	std::vector<ChemicalInterval> intervals;
	std::vector<ChemicalRefusalCertificate> refusals;
	std::optional<ChemicalPeTailCertificate> pe_tail;
	ChemicalMatrix center_error, tail_error, background_error, absolute_goal;
	std::string center_authority, tail_authority, background_error_authority;
	unsigned order = 16;
	std::vector<std::pair<NumberAxis, NumberAxis>> structural_zeros;
	std::string structural_zero_authority;
};

struct ChemicalIntegrationDiagnostics
{
	unsigned order = 16;
	std::vector<unsigned> validation_orders;
	std::vector<ChemicalMatrix> validation_integrals;
	ChemicalMatrix local_arithmetic_error, center_error, tail_error, background_error;
	ChemicalMatrix refusal_error, absolute_goal;
	std::vector<ChemicalInterval> source_intervals;
	std::vector<ChemicalRefusalCertificate> refusal_certificates;
	std::optional<ChemicalPeTailCertificate> pe_tail;
	std::string accumulation = "Neumaier compensated sums in declared segment/node order";
	std::string center_authority, tail_authority, background_error_authority;
	std::size_t validation_node_count = 0;
	std::vector<std::pair<NumberAxis, NumberAxis>> structural_zeros;
	std::string structural_zero_authority;
};

class GlobalChemicalNumberResponse
{
  public:
	GlobalChemicalNumberResponse(const GlobalChemicalNumberResponse &) = default;
	GlobalChemicalNumberResponse(GlobalChemicalNumberResponse &&) = default;
	GlobalChemicalNumberResponse &operator=(const GlobalChemicalNumberResponse &) = delete;
	GlobalChemicalNumberResponse &operator=(GlobalChemicalNumberResponse &&) = delete;
	static GlobalChemicalNumberResponse Compute(const ChemicalIntegrationRequest &);
	const ChemicalMatrix &Values() const;
	const ChemicalMatrix &NumericalError() const;
	const std::vector<NumberAxis> &Support() const;
	const std::vector<double> &Eigenvalues() const;
	void RequireCurrent() const;
	std::size_t NodeCount() const
	{
		RequireCurrent();
		return nodes_;
	}
	const std::vector<double> &Partition() const
	{
		RequireCurrent();
		return partition_;
	}
	const ChemicalMatrix &QuadratureError() const
	{
		RequireCurrent();
		return quadrature_error_;
	}
	const ChemicalIntegrationDiagnostics &Diagnostics() const
	{
		RequireCurrent();
		return diagnostics_;
	}
	double Condition() const
	{
		RequireCurrent();
		return condition_;
	}
	const std::shared_ptr<ChemicalLifetime> &Lifetime() const
	{
		RequireCurrent();
		return lifetime_;
	}

  private:
	GlobalChemicalNumberResponse() = default;
	ChemicalMatrix g_, error_, quadrature_error_;
	std::vector<NumberAxis> support_;
	std::vector<double> eigenvalues_, partition_;
	std::size_t nodes_ = 0;
	double condition_ = 0;
	std::shared_ptr<ChemicalLifetime> lifetime_;
	std::string snapshot_;
	std::vector<NumberProvenance> profiles_;
	std::vector<std::string> file_contents_;
	ChemicalIntegrationDiagnostics diagnostics_;
};

struct ChemicalReductionDiagnostic
{
	ChemicalMatrix q, numerical_error, schur_arithmetic;
	double denominator = 0, denominator_error = 0;
};

class ChemicalImbalanceResponse
{
  public:
	ChemicalImbalanceResponse(const ChemicalImbalanceResponse &) = default;
	ChemicalImbalanceResponse(ChemicalImbalanceResponse &&) = default;
	ChemicalImbalanceResponse &operator=(const ChemicalImbalanceResponse &) = delete;
	ChemicalImbalanceResponse &operator=(ChemicalImbalanceResponse &&) = delete;
	static ChemicalImbalanceResponse Compute(std::shared_ptr<const GlobalChemicalNumberResponse>,
											 const ChemicalMatrix &q_goal, const ChemicalMatrix &z_goal);
	void RequireCurrent() const;
	const ChemicalMatrix &Values() const
	{
		RequireCurrent();
		return z_;
	}
	const ChemicalMatrix &NumericalError() const
	{
		RequireCurrent();
		return error_;
	}
	const ChemicalMatrix &GlobalSolveError() const
	{
		RequireCurrent();
		return global_solve_error_;
	}
	const ChemicalMatrix &ArithmeticError() const
	{
		RequireCurrent();
		return arithmetic_error_;
	}
	ChemicalReductionDiagnostic ReductionDiagnostic() const;
	const std::vector<ImbalanceChannel> &Channels() const
	{
		RequireCurrent();
		return channels_;
	}
	double Rho() const
	{
		RequireCurrent();
		return rho_;
	}
	double Condition() const
	{
		RequireCurrent();
		return condition_;
	}
	double PaperZ(ImbalanceChannel row, ImbalanceChannel column) const;
	const std::shared_ptr<const GlobalChemicalNumberResponse> &Global() const
	{
		RequireCurrent();
		return global_;
	}

  private:
	ChemicalImbalanceResponse() = default;
	std::shared_ptr<const GlobalChemicalNumberResponse> global_;
	ChemicalMatrix z_, error_, global_solve_error_, arithmetic_error_, q_goal_, z_goal_;
	std::vector<ImbalanceChannel> channels_;
	double rho_ = 0, condition_ = 0;
};

class RotochemicalSpinDrive
{
  public:
	RotochemicalSpinDrive(const RotochemicalSpinDrive &) = default;
	RotochemicalSpinDrive(RotochemicalSpinDrive &&) = default;
	RotochemicalSpinDrive &operator=(const RotochemicalSpinDrive &) = delete;
	RotochemicalSpinDrive &operator=(RotochemicalSpinDrive &&) = delete;
	static RotochemicalSpinDrive Compute(std::shared_ptr<const ChemicalImbalanceResponse>,
										 std::shared_ptr<const FixedBaryonNumberResponse>, std::vector<double> validation_I,
										 std::vector<double> numerical_goal, std::vector<double> validation_goal);
	void RequireCurrent() const;
	const std::vector<double> &Values() const
	{
		RequireCurrent();
		return w_;
	}
	const std::vector<double> &NumericalError() const
	{
		RequireCurrent();
		return error_;
	}
	const std::vector<double> &ValidationEnvelope() const
	{
		RequireCurrent();
		return validation_;
	}
	const std::vector<double> &IPhysical() const
	{
		RequireCurrent();
		return i_;
	}
	const std::vector<double> &INumericalError() const
	{
		RequireCurrent();
		return i_error_;
	}
	const std::vector<double> &IValidationEnvelope() const
	{
		RequireCurrent();
		return i_validation_;
	}
	const std::vector<double> &ArithmeticError() const
	{
		RequireCurrent();
		return arithmetic_;
	}
	std::vector<double> Evaluate(double omega, double omega_dot) const;

  private:
	RotochemicalSpinDrive() = default;
	std::shared_ptr<const ChemicalImbalanceResponse> chemical_;
	std::shared_ptr<const FixedBaryonNumberResponse> structural_;
	std::vector<double> w_, error_, validation_, arithmetic_, i_, i_error_, i_validation_;
	std::vector<double> numerical_goal_, validation_goal_;
	std::string structural_snapshot_;
	std::vector<NumberProvenance> structural_sources_snapshot_;
};
} // namespace CompactStar::Analysis
