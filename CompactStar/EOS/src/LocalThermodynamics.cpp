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

ColdRelativisticIdealFermion::ColdRelativisticIdealFermion(
	ColdIdealFermionKind kind, double rest_mass_energy_MeV,
	double hbar_c_MeV_fm) noexcept
	: kind_(kind), rest_mass_energy_MeV_(rest_mass_energy_MeV),
	  hbar_c_MeV_fm_(hbar_c_MeV_fm)
{
}

ColdRelativisticIdealFermion ColdRelativisticIdealFermion::Neutron()
{
	return {ColdIdealFermionKind::Neutron, Zaki::Physics::NEUTRON_M_MEV,
			1.0 / Zaki::Physics::MEV_2_INV_FM};
}

ColdRelativisticIdealFermion ColdRelativisticIdealFermion::Proton()
{
	return {ColdIdealFermionKind::Proton, Zaki::Physics::PROTON_M_MEV,
			1.0 / Zaki::Physics::MEV_2_INV_FM};
}

ColdRelativisticIdealFermion ColdRelativisticIdealFermion::Electron()
{
	return {ColdIdealFermionKind::Electron, Zaki::Physics::ELECTRON_M_MEV,
			1.0 / Zaki::Physics::MEV_2_INV_FM};
}

ColdRelativisticIdealFermion ColdRelativisticIdealFermion::Muon()
{
	return {ColdIdealFermionKind::Muon, Zaki::Physics::MUON_M_MEV,
			1.0 / Zaki::Physics::MEV_2_INV_FM};
}

const char *ColdRelativisticIdealFermion::Name() const noexcept
{
	switch (kind_)
	{
	case ColdIdealFermionKind::Neutron:
		return "neutron";
	case ColdIdealFermionKind::Proton:
		return "proton";
	case ColdIdealFermionKind::Electron:
		return "electron";
	case ColdIdealFermionKind::Muon:
		return "muon";
	}
	return "unknown";
}

ColdIdealFermionEvaluation
ColdRelativisticIdealFermion::Evaluate(double number_density_fm3) const
{
	RequireFiniteNonnegative(number_density_fm3, "Ideal-fermion number density");
	if (!std::isfinite(rest_mass_energy_MeV_) || rest_mass_energy_MeV_ <= 0.0 ||
		!std::isfinite(hbar_c_MeV_fm_) || hbar_c_MeV_fm_ <= 0.0)
		throw std::runtime_error("Ideal-fermion constant authority is invalid.");

	ColdIdealFermionEvaluation result;
	result.number_density_fm3 = number_density_fm3;
	result.rest_mass_energy_MeV = rest_mass_energy_MeV_;
	result.chemical_potential_MeV = rest_mass_energy_MeV_;
	if (number_density_fm3 == 0.0)
		return result;

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
		throw std::runtime_error("Ideal-fermion evaluation left the finite analytic domain.");

	return result;
}

double ColdRelativisticIdealFermion::ChemicalPotentialMeV(double number_density_fm3) const
{
	return Evaluate(number_density_fm3).chemical_potential_MeV;
}

ColdIdealFermionValues ColdRelativisticIdealFermion::Values(double n) const
{
	RequireFiniteNonnegative(n, "Ideal-fermion number density");
	const double m = rest_mass_energy_MeV_, b = hbar_c_MeV_fm_;
	if (n == 0.0) return {0.0, 0.0, m};
	const double pi = std::acos(-1.0), k = b * std::cbrt(3.0*pi*pi*n);
	const double x = k/m, x2 = x*x, mu = std::hypot(m,k);
	const double energy = m*m*m*m * DimensionlessEnergyBracket(x)/(8*pi*pi*b*b*b);
	double pressure;
	if (x < 0.1)
	{
		// Integrate q^4/sqrt(1+q^2) term by term. The convergent
		// binomial series avoids cancellation of the leading x^5 term.
		double sum = 1.0/5.0, coefficient = 1.0, power = 1.0;
		for (int j=1; j<=12; ++j) {
			coefficient *= -(2.0*j-1.0)/(2.0*j);
			power *= x2;
			sum += coefficient*power/(5.0+2.0*j);
		}
		pressure = n*k*k/m*sum;
	}
	else
		pressure = m*m*m*m*(x*std::sqrt(1+x2)*(2*x2-3)+3*std::asinh(x)) /
			(24*pi*pi*b*b*b);
	if (!std::isfinite(energy) || !std::isfinite(pressure) || pressure <= 0)
		throw std::runtime_error("Ideal-fermion values exceed floating-point resolution.");
	return {energy, pressure, mu};
}

double ColdRelativisticIdealFermion::EnergyDensityMeVFm3(double number_density_fm3) const
{
	return Evaluate(number_density_fm3).energy_density_MeV_fm3;
}

double ColdRelativisticIdealFermion::ChemicalPotentialDerivativeMeVFm3(
	double number_density_fm3) const
{
	const auto evaluation = Evaluate(number_density_fm3);
	if (!evaluation.dchemical_potential_dn_MeV_fm3)
		throw std::runtime_error(
			"Ideal-fermion dmu/dn is unavailable at zero density; the active-species "
			"smooth-domain derivative must not be regularized.");
	return *evaluation.dchemical_potential_dn_MeV_fm3;
}

double ColdRelativisticIdealFermion::NumberDensityForChemicalPotentialFm3(
	double chemical_potential_MeV) const
{
	if (!std::isfinite(chemical_potential_MeV) ||
		chemical_potential_MeV < rest_mass_energy_MeV_)
		throw std::runtime_error(
			"Ideal-fermion chemical potential must be finite and at least its rest mass.");
	if (chemical_potential_MeV == rest_mass_energy_MeV_)
		return 0.0;

	const long double chemical_potential = chemical_potential_MeV;
	const long double mass = rest_mass_energy_MeV_;
	const long double hbar_c = hbar_c_MeV_fm_;
	const long double p_F = std::sqrt(
		(chemical_potential - mass) * (chemical_potential + mass));
	const long double pi = std::acos(-1.0L);
	const long double ratio = p_F / hbar_c;
	const double density = static_cast<double>(ratio * ratio * ratio /
											 (3.0L * pi * pi));
	if (!std::isfinite(density) || density < 0.0)
		throw std::runtime_error(
			"Ideal-fermion inverse chemical-potential relation left its finite domain.");
	return density;
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
	const auto fermion = kind_ == ColdFreeLeptonKind::Electron
		? ColdRelativisticIdealFermion::Electron()
		: ColdRelativisticIdealFermion::Muon();
	const auto generic = fermion.Evaluate(number_density_fm3);
	return {generic.number_density_fm3,
			generic.rest_mass_energy_MeV,
			generic.fermi_momentum_MeV,
			generic.chemical_potential_MeV,
			generic.energy_density_MeV_fm3,
			generic.dchemical_potential_dn_MeV_fm3};
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
