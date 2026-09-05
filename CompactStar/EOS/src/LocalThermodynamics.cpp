#include "CompactStar/EOS/LocalThermodynamics.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <Zaki/Physics/Constants.hpp>

namespace CompactStar
{
namespace
{

void RequireFiniteNonnegative(double value, const char *name)
{
	if (!std::isfinite(value) || value < 0.0)
		throw std::runtime_error(std::string(name) + " must be finite and nonnegative.");
}

long double DimensionlessEnergyBracket(long double x)
{
	// The direct closed form loses its leading x^3 term when x is very small.
	// This analytic series is the same T=0 integral, not a density clip.
	if (x < 1.0e-2L)
	{
		const long double x2 = x * x;
		return x * x2 * (8.0L / 3.0L + x2 * (4.0L / 5.0L + x2 * (-1.0L / 7.0L + x2 * (1.0L / 18.0L - x2 * 5.0L / 176.0L))));
	}

	return x * std::sqrt(1.0L + x * x) * (1.0L + 2.0L * x * x) -
		   std::asinh(x);
}

} // namespace

ChargeNeutralCompositionState
MakeChargeNeutralCompositionState(const ChargeNeutralCoordinates &coordinates)
{
	if (!std::isfinite(coordinates.n_B_fm3) || coordinates.n_B_fm3 <= 0.0)
		throw std::runtime_error("n_B must be finite and strictly positive.");
	RequireFiniteNonnegative(coordinates.n_e_fm3, "n_e");
	RequireFiniteNonnegative(coordinates.n_mu_fm3, "n_mu");

	const double n_p_fm3 = coordinates.n_e_fm3 + coordinates.n_mu_fm3;
	if (!std::isfinite(n_p_fm3) || n_p_fm3 < 0.0)
		throw std::runtime_error("Reconstructed n_p is outside the finite physical domain.");

	const double n_n_fm3 = coordinates.n_B_fm3 - n_p_fm3;
	if (!std::isfinite(n_n_fm3) || n_n_fm3 < 0.0)
		throw std::runtime_error("Reconstructed n_n is outside the finite physical domain.");

	return ChargeNeutralCompositionState(coordinates.n_B_fm3, n_n_fm3,
										 n_p_fm3, coordinates.n_e_fm3,
										 coordinates.n_mu_fm3);
}

ColdRelativisticFreeLepton::ColdRelativisticFreeLepton(
	ColdFreeLeptonKind kind, double rest_mass_energy_MeV,
	double hbar_c_MeV_fm) noexcept
	: kind_(kind), rest_mass_energy_MeV_(rest_mass_energy_MeV),
	  hbar_c_MeV_fm_(hbar_c_MeV_fm)
{
}

ColdRelativisticFreeLepton ColdRelativisticFreeLepton::Electron()
{
	return {ColdFreeLeptonKind::Electron, Zaki::Physics::ELECTRON_M_MEV,
			1.0 / Zaki::Physics::MEV_2_INV_FM};
}

ColdRelativisticFreeLepton ColdRelativisticFreeLepton::Muon()
{
	return {ColdFreeLeptonKind::Muon, Zaki::Physics::MUON_M_MEV,
			1.0 / Zaki::Physics::MEV_2_INV_FM};
}

const char *ColdRelativisticFreeLepton::Name() const noexcept
{
	return kind_ == ColdFreeLeptonKind::Electron ? "electron" : "muon";
}

ColdFreeLeptonEvaluation
ColdRelativisticFreeLepton::Evaluate(double number_density_fm3) const
{
	RequireFiniteNonnegative(number_density_fm3, "Free-lepton number density");
	if (!std::isfinite(rest_mass_energy_MeV_) || rest_mass_energy_MeV_ <= 0.0 ||
		!std::isfinite(hbar_c_MeV_fm_) || hbar_c_MeV_fm_ <= 0.0)
		throw std::runtime_error("Free-lepton constant authority is invalid.");

	ColdFreeLeptonEvaluation result;
	result.number_density_fm3 = number_density_fm3;
	result.rest_mass_energy_MeV = rest_mass_energy_MeV_;
	result.chemical_potential_MeV = rest_mass_energy_MeV_;

	if (number_density_fm3 == 0.0)
	{
		// The value limit is physical; the derivative is deliberately unavailable.
		return result;
	}

	const long double pi = std::acos(-1.0L);
	const long double n = number_density_fm3;
	const long double hbar_c = hbar_c_MeV_fm_;
	const long double mass = rest_mass_energy_MeV_;
	const long double p_F = hbar_c * std::cbrt(3.0L * pi * pi * n);
	const long double mu = std::hypot(mass, p_F);
	const long double x = p_F / mass;
	const long double bracket = DimensionlessEnergyBracket(x);
	const long double energy =
		mass * mass * mass * mass * bracket /
		(8.0L * pi * pi * hbar_c * hbar_c * hbar_c);
	const long double derivative = p_F * p_F / (3.0L * mu * n);

	result.fermi_momentum_MeV = static_cast<double>(p_F);
	result.chemical_potential_MeV = static_cast<double>(mu);
	result.energy_density_MeV_fm3 = static_cast<double>(energy);
	result.dchemical_potential_dn_MeV_fm3 = static_cast<double>(derivative);

	if (!std::isfinite(result.fermi_momentum_MeV) ||
		!std::isfinite(result.chemical_potential_MeV) ||
		!std::isfinite(result.energy_density_MeV_fm3) ||
		!std::isfinite(*result.dchemical_potential_dn_MeV_fm3))
		throw std::runtime_error("Free-lepton evaluation left the finite analytic domain.");

	return result;
}

double ColdRelativisticFreeLepton::ChemicalPotentialMeV(double number_density_fm3) const
{
	return Evaluate(number_density_fm3).chemical_potential_MeV;
}

double ColdRelativisticFreeLepton::EnergyDensityMeVFm3(double number_density_fm3) const
{
	return Evaluate(number_density_fm3).energy_density_MeV_fm3;
}

double ColdRelativisticFreeLepton::ChemicalPotentialDerivativeMeVFm3(
	double number_density_fm3) const
{
	const auto evaluation = Evaluate(number_density_fm3);
	if (!evaluation.dchemical_potential_dn_MeV_fm3)
		throw std::runtime_error(
			"Free-lepton dmu/dn is unavailable at zero density; the active-species "
			"smooth-domain derivative must not be regularized.");
	return *evaluation.dchemical_potential_dn_MeV_fm3;
}

} // namespace CompactStar
