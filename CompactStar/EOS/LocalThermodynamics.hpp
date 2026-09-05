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

/// One complete local evaluation.  Energy density includes the provider's declared rest masses.
struct LocalThermodynamicEvaluation
{
	ChargeNeutralCompositionState state;
	double energy_density_MeV_fm3;
	NeutralConjugates conjugates;
	ChargeNeutralChemicalHessian hessian;
};

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
