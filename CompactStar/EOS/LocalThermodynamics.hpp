// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file LocalThermodynamics.hpp
 * @brief Cold charge-neutral local thermodynamic contract from ADR-0010.
 *
 * This boundary is independent of stars, rotation, and secular evolution.  Its
 * canonical coordinates are x=(n_B,n_e,n_mu), in fm^-3.  The provider returns
 * the conjugates g=(mu_n,-eta_npe,-eta_npmu), in MeV, and their constrained
 * chemical Hessian d g / d x, in MeV fm^3.
 *
 * The constrained Hessian already represents the physical charge-neutral
 * response.  This API deliberately exposes neither an electrostatic projection,
 * a full charged-species susceptibility, nor any paper-B matrix or inverse.
 */

#ifndef CompactStar_LocalThermodynamics_H
#define CompactStar_LocalThermodynamics_H

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <variant>

namespace CompactStar
{

/// Canonical cold local coordinates, all in fm^-3.
struct ChargeNeutralCoordinates
{
	double n_B_fm3 = 0.0;
	double n_e_fm3 = 0.0;
	double n_mu_fm3 = 0.0;
};

/**
 * @brief Validated npe-mu state on the local charge-neutral manifold.
 *
 * Instances can only be created through MakeChargeNeutralCompositionState(),
 * which fails closed on non-finite or unphysical input.  Charge neutrality is
 * imposed by construction, never by clipping an input.
 */
class ChargeNeutralCompositionState
{
  public:
	[[nodiscard]] double BaryonDensityFm3() const noexcept { return n_B_fm3_; }
	[[nodiscard]] double NeutronDensityFm3() const noexcept { return n_n_fm3_; }
	[[nodiscard]] double ProtonDensityFm3() const noexcept { return n_p_fm3_; }
	[[nodiscard]] double ElectronDensityFm3() const noexcept { return n_e_fm3_; }
	[[nodiscard]] double MuonDensityFm3() const noexcept { return n_mu_fm3_; }

	[[nodiscard]] ChargeNeutralCoordinates Coordinates() const noexcept
	{
		return {n_B_fm3_, n_e_fm3_, n_mu_fm3_};
	}

  private:
	friend ChargeNeutralCompositionState
	MakeChargeNeutralCompositionState(const ChargeNeutralCoordinates &coordinates);

	ChargeNeutralCompositionState(double n_B_fm3, double n_n_fm3,
								  double n_p_fm3, double n_e_fm3,
								  double n_mu_fm3) noexcept
		: n_B_fm3_(n_B_fm3), n_n_fm3_(n_n_fm3),
		  n_p_fm3_(n_p_fm3), n_e_fm3_(n_e_fm3), n_mu_fm3_(n_mu_fm3)
	{
	}

	double n_B_fm3_;
	double n_n_fm3_;
	double n_p_fm3_;
	double n_e_fm3_;
	double n_mu_fm3_;
};

/// Validate canonical coordinates and reconstruct n_p=n_e+n_mu and n_n=n_B-n_p.
[[nodiscard]] ChargeNeutralCompositionState
MakeChargeNeutralCompositionState(const ChargeNeutralCoordinates &coordinates);

/// g=(mu_n,-eta_npe,-eta_npmu), all in MeV.
struct NeutralConjugates
{
	std::array<double, 3> value_MeV{};

	[[nodiscard]] double MuNMeV() const noexcept { return value_MeV[0]; }
	[[nodiscard]] double EtaNpeMeV() const noexcept { return -value_MeV[1]; }
	[[nodiscard]] double EtaNpmuMeV() const noexcept { return -value_MeV[2]; }
};

/// H_ab=partial g_a/partial x_b at fixed remaining x coordinates, in MeV fm^3.
struct ChargeNeutralChemicalHessian
{
	using Matrix = std::array<std::array<double, 3>, 3>;
	Matrix value_MeV_fm3{};

	[[nodiscard]] double operator()(std::size_t row, std::size_t column) const
	{
		return value_MeV_fm3.at(row).at(column);
	}
};

/// Species active in the returned chart (not a mask over a larger matrix).
/// A zero-density species at its appearance threshold is not yet active.
struct ActiveParticleContent
{
	bool neutron;
	bool proton;
	bool electron;
	bool muon;
};

enum class LocalResponseDomainStatus
{
	SmoothInterior,
	SpeciesThreshold,
	VacuumBoundary
};

/// One complete smooth three-coordinate local evaluation.
/// Energy density includes the provider's declared rest masses.
struct LocalThermodynamicEvaluation
{
	ChargeNeutralCompositionState state;
	double energy_density_MeV_fm3;
	NeutralConjugates conjugates;
	ChargeNeutralChemicalHessian hessian;

	static constexpr ActiveParticleContent active_particles{true, true, true, true};
	static constexpr std::size_t response_dimension = 3;
	static constexpr auto domain_status = LocalResponseDomainStatus::SmoothInterior;
};

/// Explicit active chart z=(n_B,n_e) in fm^-3; n_mu is identically zero.
struct NpeCoordinates
{
	double n_B_fm3;
	double n_e_fm3;
};

/// h=(mu_n,-eta_npe) in MeV. No active muon conjugate exists in this chart.
struct NpeConjugates
{
	std::array<double, 2> value_MeV{};
	[[nodiscard]] double MuNMeV() const noexcept { return value_MeV[0]; }
	[[nodiscard]] double EtaNpeMeV() const noexcept { return -value_MeV[1]; }
};

/// d h / d z, holding the other active coordinate fixed; units MeV fm^3.
struct NpeChemicalHessian
{
	using Matrix = std::array<std::array<double, 2>, 2>;
	Matrix value_MeV_fm3{};
	[[nodiscard]] double operator()(std::size_t row, std::size_t column) const
	{
		return value_MeV_fm3.at(row).at(column);
	}
};

struct NpeThermodynamicEvaluation
{
	ChargeNeutralCompositionState state;
	double energy_density_MeV_fm3;
	NpeConjugates conjugates;
	NpeChemicalHessian hessian;
	// Value diagnostic only: mu_n-mu_p-m_mu at zero muon density.
	// It is deliberately outside the active conjugates and Hessian.
	double eta_npmu_threshold_diagnostic_MeV;

	static constexpr ActiveParticleContent active_particles{true, true, true, false};
	static constexpr std::size_t response_dimension = 2;
	static constexpr auto domain_status = LocalResponseDomainStatus::SmoothInterior;
};

/// Explicit active chart z_pe=(n_B) in fm^-3; n_p=n_e=n_B and n_n=n_mu=0.
struct PeCoordinates
{
	double n_B_fm3;
};

/// h_pe=d epsilon_pe/d n_B=mu_p+mu_e in MeV.
struct PeConjugates
{
	std::array<double, 1> value_MeV{};
	[[nodiscard]] double HPeMeV() const noexcept { return value_MeV[0]; }
};

/// d h_pe/d n_B in MeV fm^3. This is a genuine 1x1 active response.
struct PeChemicalHessian
{
	using Matrix = std::array<std::array<double, 1>, 1>;
	Matrix value_MeV_fm3{};
	[[nodiscard]] double operator()(std::size_t row, std::size_t column) const
	{
		return value_MeV_fm3.at(row).at(column);
	}
};

struct PeThermodynamicEvaluation
{
	ChargeNeutralCompositionState state;
	double energy_density_MeV_fm3;
	PeConjugates conjugates;
	PeChemicalHessian hessian;
	// Value diagnostic only: m_n-mu_p-mu_e while the neutron is inactive.
	double eta_npe_threshold_diagnostic_MeV;

	static constexpr ActiveParticleContent active_particles{false, true, true, false};
	static constexpr std::size_t response_dimension = 1;
	static constexpr auto domain_status = LocalResponseDomainStatus::SmoothInterior;
};

/// Value-only neutron appearance: p/e values exist, but the active chart changes 1D -> 2D.
struct NeutronThresholdEvaluation
{
	ChargeNeutralCompositionState state;
	double energy_density_MeV_fm3;
	PeConjugates limiting_pe_conjugates;
	double eta_npe_threshold_diagnostic_MeV;

	static constexpr ActiveParticleContent active_particles{false, true, true, false};
	// Zero means no smooth response returned, not a zero-dimensional matrix.
	static constexpr std::size_t response_dimension = 0;
	static constexpr auto domain_status = LocalResponseDomainStatus::SpeciesThreshold;
};

/// All-zero physical composition at the n_B=0 boundary.
struct VacuumCompositionState
{
	[[nodiscard]] constexpr double BaryonDensityFm3() const noexcept { return 0.0; }
	[[nodiscard]] constexpr double NeutronDensityFm3() const noexcept { return 0.0; }
	[[nodiscard]] constexpr double ProtonDensityFm3() const noexcept { return 0.0; }
	[[nodiscard]] constexpr double ElectronDensityFm3() const noexcept { return 0.0; }
	[[nodiscard]] constexpr double MuonDensityFm3() const noexcept { return 0.0; }
	[[nodiscard]] constexpr ChargeNeutralCoordinates Coordinates() const noexcept
	{
		return {0.0, 0.0, 0.0};
	}
};

/// Value-only vacuum boundary. The one-sided pe conjugate limit is finite; its Hessian is not.
struct VacuumBoundaryEvaluation
{
	VacuumCompositionState state;
	double energy_density_MeV_fm3;
	PeConjugates limiting_pe_conjugates;
	double eta_npe_threshold_diagnostic_MeV;

	static constexpr ActiveParticleContent active_particles{false, false, false, false};
	static constexpr std::size_t response_dimension = 0;
	static constexpr auto domain_status = LocalResponseDomainStatus::VacuumBoundary;
};

/// Value-only muon-appearance result: no smooth Hessian of either size.
struct MuonThresholdEvaluation
{
	ChargeNeutralCompositionState state;
	double energy_density_MeV_fm3;
	NpeConjugates limiting_npe_conjugates;
	double eta_npmu_threshold_diagnostic_MeV;

	static constexpr ActiveParticleContent active_particles{true, true, true, false};
	// Zero means no smooth response returned, not a zero-dimensional matrix.
	static constexpr std::size_t response_dimension = 0;
	static constexpr auto domain_status = LocalResponseDomainStatus::SpeciesThreshold;
};

/// A consumer must inspect the alternative before accessing a Hessian.
/// There is no conversion from a lower-dimensional or boundary alternative to the full one.
using ActiveLocalThermodynamicEvaluation = std::variant<
	LocalThermodynamicEvaluation, NpeThermodynamicEvaluation, MuonThresholdEvaluation,
	PeThermodynamicEvaluation, NeutronThresholdEvaluation, VacuumBoundaryEvaluation>;

/// Minimum provenance and domain description required of a local provider.
struct LocalThermodynamicProviderMetadata
{
	std::string model_id;
	std::string model_revision;
	std::string particle_content;
	std::string coordinate_chart;
	std::string temperature_scope;
	std::string rest_mass_convention;
	std::string lepton_ownership;
	std::string smooth_domain;
};

/**
 * @brief Generic ADR-0010 cold local thermodynamic provider.
 *
 * Individual charged intrinsic chemical potentials are intentionally not part
 * of this minimum contract.  A model may expose them through its own optional
 * API, but consumers of this interface use only the neutral conjugates.
 */
class ILocalThermodynamicProvider
{
  public:
	virtual ~ILocalThermodynamicProvider() = default;

	[[nodiscard]] virtual const LocalThermodynamicProviderMetadata &Metadata() const noexcept = 0;
	[[nodiscard]] virtual LocalThermodynamicEvaluation
	Evaluate(const ChargeNeutralCoordinates &coordinates) const = 0;
	[[nodiscard]] virtual ChargeNeutralCompositionState
	EquilibriumStateAt(double n_B_fm3) const = 0;

	/// Generic active-response entry point; existing smooth providers keep their semantics.
	[[nodiscard]] virtual ActiveLocalThermodynamicEvaluation
	EvaluateActive(const ChargeNeutralCoordinates &coordinates) const
	{
		return Evaluate(coordinates);
	}
	/// Composition plus per-result active species, response dimension and domain status.
	[[nodiscard]] virtual ActiveLocalThermodynamicEvaluation
	EquilibriumAt(double n_B_fm3) const
	{
		return EvaluateActive(EquilibriumStateAt(n_B_fm3).Coordinates());
	}
};

enum class ColdFreeLeptonKind
{
	Electron,
	Muon
};

enum class ColdIdealFermionKind
{
	Neutron,
	Proton,
	Electron,
	Muon
};

/**
 * @brief Cold relativistic ideal spin-1/2 fermion result.
 *
 * The total energy density includes rest mass.  The derivative is deliberately
 * unavailable at zero density, where the active-species Hessian is singular.
 */
struct ColdIdealFermionEvaluation
{
	double number_density_fm3 = 0.0;
	double rest_mass_energy_MeV = 0.0;
	double fermi_momentum_MeV = 0.0;
	double chemical_potential_MeV = 0.0;
	double energy_density_MeV_fm3 = 0.0;
	std::optional<double> dchemical_potential_dn_MeV_fm3;
};

/**
 * @brief Analytic T=0, spin-1/2, noninteracting relativistic fermion.
 *
 * The factories bind species masses and hbar-c to the repository's Zaki
 * constant authority.  This primitive has no equilibrium or stellar policy.
 */
class ColdRelativisticIdealFermion
{
  public:
	[[nodiscard]] static ColdRelativisticIdealFermion Neutron();
	[[nodiscard]] static ColdRelativisticIdealFermion Proton();
	[[nodiscard]] static ColdRelativisticIdealFermion Electron();
	[[nodiscard]] static ColdRelativisticIdealFermion Muon();

	[[nodiscard]] ColdIdealFermionKind Kind() const noexcept { return kind_; }
	[[nodiscard]] const char *Name() const noexcept;
	[[nodiscard]] double RestMassEnergyMeV() const noexcept { return rest_mass_energy_MeV_; }
	[[nodiscard]] double HbarCMeVFm() const noexcept { return hbar_c_MeV_fm_; }

	[[nodiscard]] ColdIdealFermionEvaluation Evaluate(double number_density_fm3) const;
	[[nodiscard]] double ChemicalPotentialMeV(double number_density_fm3) const;
	[[nodiscard]] double EnergyDensityMeVFm3(double number_density_fm3) const;
	[[nodiscard]] double ChemicalPotentialDerivativeMeVFm3(double number_density_fm3) const;

	/// Inverse of the positive-density chemical-potential relation; zero at mu=m.
	[[nodiscard]] double NumberDensityForChemicalPotentialFm3(
		double chemical_potential_MeV) const;

  private:
	ColdRelativisticIdealFermion(ColdIdealFermionKind kind,
								 double rest_mass_energy_MeV,
								 double hbar_c_MeV_fm) noexcept;

	ColdIdealFermionKind kind_;
	double rest_mass_energy_MeV_;
	double hbar_c_MeV_fm_;
};

/**
 * @brief Cold relativistic free-Fermi lepton result.
 *
 * Momentum and chemical potential include c and are reported in MeV.  Energy
 * density includes the rest mass and is in MeV fm^-3.  dmu/dn is in MeV fm^3;
 * it is explicitly unavailable at n=0 rather than silently regularized.
 */
struct ColdFreeLeptonEvaluation
{
	double number_density_fm3 = 0.0;
	double rest_mass_energy_MeV = 0.0;
	double fermi_momentum_MeV = 0.0;
	double chemical_potential_MeV = 0.0;
	double energy_density_MeV_fm3 = 0.0;
	std::optional<double> dchemical_potential_dn_MeV_fm3;
};

/**
 * @brief Analytic T=0, spin-1/2, noninteracting electron or muon component.
 */
class ColdRelativisticFreeLepton
{
  public:
	[[nodiscard]] static ColdRelativisticFreeLepton Electron();
	[[nodiscard]] static ColdRelativisticFreeLepton Muon();

	[[nodiscard]] ColdFreeLeptonKind Kind() const noexcept { return kind_; }
	[[nodiscard]] const char *Name() const noexcept;
	[[nodiscard]] double RestMassEnergyMeV() const noexcept { return rest_mass_energy_MeV_; }
	[[nodiscard]] double HbarCMeVFm() const noexcept { return hbar_c_MeV_fm_; }

	[[nodiscard]] ColdFreeLeptonEvaluation Evaluate(double number_density_fm3) const;
	[[nodiscard]] double ChemicalPotentialMeV(double number_density_fm3) const;
	[[nodiscard]] double EnergyDensityMeVFm3(double number_density_fm3) const;

	/// @throws std::runtime_error at n<=0 or for any non-finite/invalid state.
	[[nodiscard]] double ChemicalPotentialDerivativeMeVFm3(double number_density_fm3) const;

  private:
	ColdRelativisticFreeLepton(ColdFreeLeptonKind kind,
							   double rest_mass_energy_MeV,
							   double hbar_c_MeV_fm) noexcept;

	ColdFreeLeptonKind kind_;
	double rest_mass_energy_MeV_;
	double hbar_c_MeV_fm_;
};

} // namespace CompactStar

#endif /* CompactStar_LocalThermodynamics_H */
