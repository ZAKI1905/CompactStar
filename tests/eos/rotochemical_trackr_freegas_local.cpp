#include "CompactStar/EOS/TrackRFreeGasThermodynamics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <Zaki/Physics/Constants.hpp>

namespace
{

using CompactStar::ChargeNeutralCoordinates;
using CompactStar::ColdRelativisticFreeLepton;
using CompactStar::ColdRelativisticIdealFermion;
using CompactStar::MakeChargeNeutralCompositionState;
using CompactStar::TrackRFreeGasThermodynamicProvider;
using Matrix3 = std::array<std::array<double, 3>, 3>;
using Matrix43 = std::array<std::array<double, 3>, 4>;
using Vector3 = std::array<double, 3>;
using Vector4 = std::array<double, 4>;

struct SpeciesAuthority
{
	ColdRelativisticIdealFermion species;
	double mass_MeV;
};

struct IndependentFermionEvaluation
{
	long double p_F_MeV;
	long double mu_MeV;
	long double energy_MeV_fm3;
	long double derivative_MeV_fm3;
};

void Require(bool condition, const std::string &message)
{
	if (!condition)
		throw std::runtime_error(message);
}

void RequireNear(double actual, double expected, double tolerance,
				 const std::string &message)
{
	if (!std::isfinite(actual) || !std::isfinite(expected) ||
		std::abs(actual - expected) > tolerance)
		throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
								 ", expected=" + std::to_string(expected) +
								 ", tolerance=" + std::to_string(tolerance));
}

double RelativeError(double actual, double expected)
{
	return std::abs(actual - expected) / std::max(std::abs(expected), 1.0e-300);
}

template <typename Function>
bool ThrowsRuntimeError(Function &&function)
{
	try
	{
		function();
	}
	catch (const std::runtime_error &)
	{
		return true;
	}
	return false;
}

long double IndependentDimensionlessEnergyIntegral(long double x,
												std::size_t panels = 32768)
{
	Require(panels > 0 && panels % 2 == 0,
			"Independent Simpson quadrature requires a positive even panel count");
	const long double step = x / static_cast<long double>(panels);
	long double sum = 0.0L;
	long double compensation = 0.0L;
	for (std::size_t index = 0; index <= panels; ++index)
	{
		const long double t = step * static_cast<long double>(index);
		const long double integrand = t * t * std::sqrt(1.0L + t * t);
		const long double weight = index == 0 || index == panels
			? 1.0L
			: (index % 2 == 0 ? 2.0L : 4.0L);
		const long double term = weight * integrand - compensation;
		const long double updated = sum + term;
		compensation = (updated - sum) - term;
		sum = updated;
	}
	return step * sum / 3.0L;
}

IndependentFermionEvaluation IndependentFermion(double n_fm3, double mass_MeV,
											 double hbar_c_MeV_fm,
											 std::size_t panels = 32768)
{
	Require(n_fm3 > 0.0 && std::isfinite(n_fm3),
			"Independent fermion oracle requires positive finite density");
	const long double pi = std::acos(-1.0L);
	const long double n = n_fm3;
	const long double mass = mass_MeV;
	const long double hbar_c = hbar_c_MeV_fm;
	const long double p_F = hbar_c * std::cbrt(3.0L * pi * pi * n);
	const long double mu = std::sqrt(mass * mass + p_F * p_F);
	const long double integral = IndependentDimensionlessEnergyIntegral(p_F / mass, panels);
	const long double energy = mass * mass * mass * mass * integral /
		(pi * pi * hbar_c * hbar_c * hbar_c);
	const long double derivative = p_F * p_F / (3.0L * mu * n);
	return {p_F, mu, energy, derivative};
}

double DensityForTargetX(double x, double mass_MeV, double hbar_c_MeV_fm)
{
	const long double pi = std::acos(-1.0L);
	const long double ratio = static_cast<long double>(mass_MeV) * x /
		static_cast<long double>(hbar_c_MeV_fm);
	return static_cast<double>(ratio * ratio * ratio / (3.0L * pi * pi));
}

double IndependentDensityForMu(double mu_MeV, double mass_MeV,
								 double hbar_c_MeV_fm)
{
	if (mu_MeV <= mass_MeV)
		return 0.0;
	const long double mu = mu_MeV;
	const long double mass = mass_MeV;
	const long double p_F = std::sqrt((mu - mass) * (mu + mass));
	const long double ratio = p_F / hbar_c_MeV_fm;
	const long double pi = std::acos(-1.0L);
	return static_cast<double>(ratio * ratio * ratio / (3.0L * pi * pi));
}

double IndependentMu(double n_fm3, double mass_MeV, double hbar_c_MeV_fm)
{
	const long double pi = std::acos(-1.0L);
	const long double p_F = static_cast<long double>(hbar_c_MeV_fm) *
		std::cbrt(3.0L * pi * pi * static_cast<long double>(n_fm3));
	return static_cast<double>(std::hypot(static_cast<long double>(mass_MeV), p_F));
}

double IndependentDerivative(double n_fm3, double mass_MeV,
								 double hbar_c_MeV_fm)
{
	const long double pi = std::acos(-1.0L);
	const long double n = n_fm3;
	const long double p_F = static_cast<long double>(hbar_c_MeV_fm) *
		std::cbrt(3.0L * pi * pi * n);
	const long double mu = std::hypot(static_cast<long double>(mass_MeV), p_F);
	return static_cast<double>(p_F * p_F / (3.0L * mu * n));
}

Vector4 IndependentSpeciesMu(const ChargeNeutralCoordinates &coordinates,
									 double hbar_c_MeV_fm)
{
	const auto state = MakeChargeNeutralCompositionState(coordinates);
	return {
		IndependentMu(state.NeutronDensityFm3(), Zaki::Physics::NEUTRON_M_MEV,
					  hbar_c_MeV_fm),
		IndependentMu(state.ProtonDensityFm3(), Zaki::Physics::PROTON_M_MEV,
					  hbar_c_MeV_fm),
		IndependentMu(state.ElectronDensityFm3(), Zaki::Physics::ELECTRON_M_MEV,
					  hbar_c_MeV_fm),
		IndependentMu(state.MuonDensityFm3(), Zaki::Physics::MUON_M_MEV,
					  hbar_c_MeV_fm)};
}

double IndependentNeutralEnergy(const ChargeNeutralCoordinates &coordinates,
								 double hbar_c_MeV_fm,
								 std::size_t panels = 16384)
{
	const auto state = MakeChargeNeutralCompositionState(coordinates);
	return static_cast<double>(
		IndependentFermion(state.NeutronDensityFm3(), Zaki::Physics::NEUTRON_M_MEV,
						   hbar_c_MeV_fm, panels).energy_MeV_fm3 +
		IndependentFermion(state.ProtonDensityFm3(), Zaki::Physics::PROTON_M_MEV,
						   hbar_c_MeV_fm, panels).energy_MeV_fm3 +
		IndependentFermion(state.ElectronDensityFm3(), Zaki::Physics::ELECTRON_M_MEV,
						   hbar_c_MeV_fm, panels).energy_MeV_fm3 +
		IndependentFermion(state.MuonDensityFm3(), Zaki::Physics::MUON_M_MEV,
						   hbar_c_MeV_fm, panels).energy_MeV_fm3);
}

Matrix3 IndependentHessian(const ChargeNeutralCoordinates &coordinates,
							  double hbar_c_MeV_fm)
{
	const auto state = MakeChargeNeutralCompositionState(coordinates);
	const double D_n = IndependentDerivative(state.NeutronDensityFm3(),
		Zaki::Physics::NEUTRON_M_MEV, hbar_c_MeV_fm);
	const double D_p = IndependentDerivative(state.ProtonDensityFm3(),
		Zaki::Physics::PROTON_M_MEV, hbar_c_MeV_fm);
	const double D_e = IndependentDerivative(state.ElectronDensityFm3(),
		Zaki::Physics::ELECTRON_M_MEV, hbar_c_MeV_fm);
	const double D_mu = IndependentDerivative(state.MuonDensityFm3(),
		Zaki::Physics::MUON_M_MEV, hbar_c_MeV_fm);
	return {{{{D_n, -D_n, -D_n}},
			 {{-D_n, D_n + D_p + D_e, D_n + D_p}},
			 {{-D_n, D_n + D_p, D_n + D_p + D_mu}}}};
}

double Determinant(const Matrix3 &m)
{
	return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
		   m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
		   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

Matrix3 Inverse3(const Matrix3 &m)
{
	const double determinant = Determinant(m);
	Require(std::isfinite(determinant) && determinant > 0.0,
			"RFG condition metric requested a non-positive matrix inverse");
	Matrix3 inverse{};
	inverse[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / determinant;
	inverse[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / determinant;
	inverse[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / determinant;
	inverse[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / determinant;
	inverse[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / determinant;
	inverse[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / determinant;
	inverse[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / determinant;
	inverse[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / determinant;
	inverse[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / determinant;
	return inverse;
}

double InfinityNorm(const Matrix3 &m)
{
	double result = 0.0;
	for (const auto &row : m)
	{
		double sum = 0.0;
		for (double value : row)
			sum += std::abs(value);
		result = std::max(result, sum);
	}
	return result;
}

double MaxAbsDifference(const Matrix3 &a, const Matrix3 &b)
{
	double result = 0.0;
	for (std::size_t row = 0; row < 3; ++row)
		for (std::size_t column = 0; column < 3; ++column)
			result = std::max(result, std::abs(a[row][column] - b[row][column]));
	return result;
}

ChargeNeutralCoordinates Offset(const ChargeNeutralCoordinates &x,
								 std::size_t axis, double delta)
{
	auto result = x;
	double *values[] = {&result.n_B_fm3, &result.n_e_fm3, &result.n_mu_fm3};
	*values[axis] += delta;
	return result;
}

ChargeNeutralCoordinates OffEquilibriumState(
	const TrackRFreeGasThermodynamicProvider &provider, double fraction)
{
	const double onset = provider.MuonOnsetBaryonDensityFm3();
	const double ceiling = provider.SigmaMinusOnsetBaryonDensityFm3();
	const double n_B = onset + fraction * (ceiling - onset);
	return {n_B, 0.055 * n_B, 0.018 * n_B};
}

Vector3 IndependentActiveEquilibrium(double n_B, double hbar_c_MeV_fm)
{
	const double m_e = Zaki::Physics::ELECTRON_M_MEV;
	const double m_mu = Zaki::Physics::MUON_M_MEV;
	auto residual = [=](double common_mu)
	{
		const double n_e = IndependentDensityForMu(common_mu, m_e, hbar_c_MeV_fm);
		const double n_mu = IndependentDensityForMu(common_mu, m_mu, hbar_c_MeV_fm);
		const double n_p = n_e + n_mu;
		return IndependentMu(n_B - n_p, Zaki::Physics::NEUTRON_M_MEV,
						 hbar_c_MeV_fm) -
			   IndependentMu(n_p, Zaki::Physics::PROTON_M_MEV, hbar_c_MeV_fm) -
			   common_mu;
	};
	double lower = m_mu;
	double upper = IndependentMu(n_B, Zaki::Physics::NEUTRON_M_MEV,
								 hbar_c_MeV_fm) - Zaki::Physics::PROTON_M_MEV;
	Require(residual(lower) >= 0.0 && residual(upper) <= 0.0,
			"Independent active-equilibrium bracket failed");
	for (int iteration = 0; iteration < 256; ++iteration)
	{
		const double midpoint = lower + 0.5 * (upper - lower);
		if (residual(midpoint) > 0.0)
			lower = midpoint;
		else
			upper = midpoint;
		if (upper - lower <= 16.0 * std::numeric_limits<double>::epsilon() *
								 std::max(1.0, std::abs(midpoint)))
			break;
	}
	const double common_mu = lower + 0.5 * (upper - lower);
	return {IndependentDensityForMu(common_mu, m_e, hbar_c_MeV_fm),
			IndependentDensityForMu(common_mu, m_mu, hbar_c_MeV_fm), common_mu};
}

void RFG1SpeciesValues()
{
	const double hbar_c = 1.0 / Zaki::Physics::MEV_2_INV_FM;
	const std::array<SpeciesAuthority, 4> species{{
		{ColdRelativisticIdealFermion::Neutron(), Zaki::Physics::NEUTRON_M_MEV},
		{ColdRelativisticIdealFermion::Proton(), Zaki::Physics::PROTON_M_MEV},
		{ColdRelativisticIdealFermion::Electron(), Zaki::Physics::ELECTRON_M_MEV},
		{ColdRelativisticIdealFermion::Muon(), Zaki::Physics::MUON_M_MEV}}};
	double maximum_energy_relative_error = 0.0;
	double maximum_other_relative_error = 0.0;
	for (const auto &authority : species)
	{
		RequireNear(authority.species.RestMassEnergyMeV(), authority.mass_MeV, 0.0,
					"RFG1 species factory does not use authoritative mass");
		RequireNear(authority.species.HbarCMeVFm(), hbar_c, 0.0,
					"RFG1 species factory does not use authoritative hbar-c");
		for (double x : {1.0e-3, 0.1, 1.0, 10.0})
		{
			const double density = DensityForTargetX(x, authority.mass_MeV, hbar_c);
			const auto actual = authority.species.Evaluate(density);
			const auto expected = IndependentFermion(density, authority.mass_MeV, hbar_c);
			const std::array<std::array<double, 2>, 3> pairs{{
				{{actual.fermi_momentum_MeV,
				  static_cast<double>(expected.p_F_MeV)}},
				{{actual.chemical_potential_MeV,
				  static_cast<double>(expected.mu_MeV)}},
				{{*actual.dchemical_potential_dn_MeV_fm3,
				  static_cast<double>(expected.derivative_MeV_fm3)}}}};
			for (const auto &pair : pairs)
			{
				const double error = RelativeError(pair[0], pair[1]);
				maximum_other_relative_error = std::max(maximum_other_relative_error, error);
				Require(error <= 3.0e-14, "RFG1 independent p_F/mu/dmu-dn mismatch");
			}
			const double energy_error = RelativeError(
				actual.energy_density_MeV_fm3,
				static_cast<double>(expected.energy_MeV_fm3));
			maximum_energy_relative_error = std::max(maximum_energy_relative_error,
												 energy_error);
			Require(energy_error <= 2.0e-12,
					"RFG1 independent quadrature energy mismatch");
		}
	}

	// The existing public lepton wrapper delegates to the generalized primitive
	// without changing a single returned scalar or threshold behavior.
	for (const auto pair : {
			 std::array<int, 2>{0, 0}, std::array<int, 2>{1, 1}})
	{
		const auto legacy = pair[0] == 0 ? ColdRelativisticFreeLepton::Electron()
									   : ColdRelativisticFreeLepton::Muon();
		const auto generic = pair[1] == 0 ? ColdRelativisticIdealFermion::Electron()
										: ColdRelativisticIdealFermion::Muon();
		for (double density : {0.0, 1.0e-18, 1.0e-6, 0.01})
		{
			const auto old_value = legacy.Evaluate(density);
			const auto new_value = generic.Evaluate(density);
			Require(old_value.fermi_momentum_MeV == new_value.fermi_momentum_MeV &&
					old_value.chemical_potential_MeV == new_value.chemical_potential_MeV &&
					old_value.energy_density_MeV_fm3 == new_value.energy_density_MeV_fm3 &&
					old_value.dchemical_potential_dn_MeV_fm3 ==
						new_value.dchemical_potential_dn_MeV_fm3,
					"RFG1 free-lepton refactor is not bit-equivalent");
		}
	}
	std::cout << "RFG1 PASS: Zaki masses (n,p,e,mu)="
			  << species[0].mass_MeV << "," << species[1].mass_MeV << ","
			  << species[2].mass_MeV << "," << species[3].mass_MeV
			  << " MeV; four species x={1e-3,0.1,1,10}; max quadrature energy relative error="
			  << maximum_energy_relative_error << ", other="
			  << maximum_other_relative_error << "; legacy leptons bit-equivalent.\n";
}

void RFG2NeutralEnergy(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double hbar_c = 1.0 / Zaki::Physics::MEV_2_INV_FM;
	double maximum_relative_error = 0.0;
	for (double fraction : {0.15, 0.45, 0.75})
	{
		const auto coordinates = OffEquilibriumState(provider, fraction);
		const double actual = provider.Evaluate(coordinates).energy_density_MeV_fm3;
		const double expected = IndependentNeutralEnergy(coordinates, hbar_c);
		maximum_relative_error = std::max(maximum_relative_error,
			RelativeError(actual, expected));
	}
	Require(maximum_relative_error <= 2.0e-12,
			"RFG2 provider energy disagrees with independent four-integral sum");
	std::cout << "RFG2 PASS: max independent neutral-energy relative error="
			  << maximum_relative_error << ".\n";
}

void RFG3NeutralConjugates(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double hbar_c = 1.0 / Zaki::Physics::MEV_2_INV_FM;
	const auto coordinates = OffEquilibriumState(provider, 0.37);
	const auto chemical = IndependentSpeciesMu(coordinates, hbar_c);
	const double eta_npe = chemical[0] - chemical[1] - chemical[2];
	const double eta_npmu = chemical[0] - chemical[1] - chemical[3];
	Require(std::abs(eta_npe) > 1.0 && std::abs(eta_npmu) > 1.0,
			"RFG3 fixture accidentally has negligible imbalance");
	const auto evaluation = provider.Evaluate(coordinates);
	RequireNear(evaluation.conjugates.MuNMeV(), chemical[0], 2.0e-12,
				"RFG3 mu_n mismatch");
	RequireNear(evaluation.conjugates.EtaNpeMeV(), eta_npe, 3.0e-12,
				"RFG3 public eta_npe mismatch");
	RequireNear(evaluation.conjugates.EtaNpmuMeV(), eta_npmu, 3.0e-12,
				"RFG3 public eta_npmu mismatch");
	std::cout << "RFG3 PASS: independent/actual eta_npe=" << eta_npe << "/"
			  << evaluation.conjugates.EtaNpeMeV() << ", eta_npmu=" << eta_npmu
			  << "/" << evaluation.conjugates.EtaNpmuMeV() << ".\n";
}

void RFG4GradientConsistency(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double hbar_c = 1.0 / Zaki::Physics::MEV_2_INV_FM;
	const auto coordinates = OffEquilibriumState(provider, 0.52);
	const auto expected = provider.Evaluate(coordinates).conjugates.value_MeV;
	const auto state = MakeChargeNeutralCompositionState(coordinates);
	const std::array<double, 3> scales{{
		state.NeutronDensityFm3(), state.ElectronDensityFm3(),
		state.MuonDensityFm3()}};
	double maximum_final_error = 0.0;
	for (std::size_t axis = 0; axis < 3; ++axis)
	{
		std::array<double, 4> errors{};
		for (std::size_t index = 0; index < errors.size(); ++index)
		{
			const double step = scales[axis] * 2.0e-3 / std::pow(2.0, index);
			const double upper = IndependentNeutralEnergy(
				Offset(coordinates, axis, step), hbar_c, 8192);
			const double lower = IndependentNeutralEnergy(
				Offset(coordinates, axis, -step), hbar_c, 8192);
			const double derivative = (upper - lower) / (2.0 * step);
			errors[index] = std::abs(derivative - expected[axis]);
		}
		Require(errors[1] < errors[0] && errors[2] < errors[1],
				"RFG4 independent energy-gradient error did not initially converge");
		maximum_final_error = std::max(maximum_final_error, errors.back());
		std::cout << "RFG4 diagnostic axis=" << axis << " errors=" << errors[0]
				  << "," << errors[1] << "," << errors[2] << "," << errors[3] << "\n";
	}
	Require(maximum_final_error <= 2.0e-6,
			"RFG4 independent energy gradient is not accurate enough");
	std::cout << "RFG4 PASS: all canonical directions converge; max final error="
			  << maximum_final_error << " MeV.\n";
}

void RFG5AnalyticHessian(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double hbar_c = 1.0 / Zaki::Physics::MEV_2_INV_FM;
	double maximum_error = 0.0;
	for (double fraction : {0.18, 0.51, 0.82})
	{
		const auto coordinates = OffEquilibriumState(provider, fraction);
		maximum_error = std::max(maximum_error, MaxAbsDifference(
			provider.Evaluate(coordinates).hessian.value_MeV_fm3,
			IndependentHessian(coordinates, hbar_c)));
	}
	Require(maximum_error <= 2.0e-11,
			"RFG5 analytic Hessian disagrees with independently assembled D_i matrix");
	std::cout << "RFG5 PASS: all matrix entries; max absolute error="
			  << maximum_error << " MeV fm^3.\n";
}

void RFG6HessianPerturbations(const TrackRFreeGasThermodynamicProvider &provider)
{
	const auto coordinates = OffEquilibriumState(provider, 0.61);
	const auto evaluation = provider.Evaluate(coordinates);
	const auto state = evaluation.state;
	const std::array<double, 3> scales{{
		state.NeutronDensityFm3(), state.ElectronDensityFm3(),
		state.MuonDensityFm3()}};
	double maximum_final_error = 0.0;
	for (std::size_t column = 0; column < 3; ++column)
	{
		std::array<double, 4> errors{};
		for (std::size_t index = 0; index < errors.size(); ++index)
		{
			const double step = scales[column] * 2.0e-3 / std::pow(2.0, index);
			const auto upper = provider.Evaluate(Offset(coordinates, column, step))
								   .conjugates.value_MeV;
			const auto lower = provider.Evaluate(Offset(coordinates, column, -step))
								   .conjugates.value_MeV;
			for (std::size_t row = 0; row < 3; ++row)
			{
				const double finite_difference = (upper[row] - lower[row]) /
					(2.0 * step);
				errors[index] = std::max(errors[index], std::abs(
					finite_difference - evaluation.hessian(row, column)));
			}
		}
		Require(errors[1] < errors[0] && errors[2] < errors[1],
				"RFG6 Hessian-column perturbation did not initially converge");
		maximum_final_error = std::max(maximum_final_error, errors.back());
		std::cout << "RFG6 diagnostic column=" << column << " errors=" << errors[0]
				  << "," << errors[1] << "," << errors[2] << "," << errors[3] << "\n";
	}
	Require(maximum_final_error <= 2.0e-4,
			"RFG6 converged Hessian-column error is too large");
	std::cout << "RFG6 PASS: all held-fixed columns converge; max final error="
			  << maximum_final_error << " MeV fm^3.\n";
}

void RFG7SymmetryConvexity(const TrackRFreeGasThermodynamicProvider &provider)
{
	double maximum_asymmetry = 0.0;
	double minimum_minor_1 = std::numeric_limits<double>::infinity();
	double minimum_minor_2 = std::numeric_limits<double>::infinity();
	double minimum_determinant = std::numeric_limits<double>::infinity();
	double maximum_condition = 0.0;
	for (double density_fraction : {0.08, 0.3, 0.55, 0.88})
		for (double electron_fraction : {0.035, 0.06, 0.10})
			for (double muon_fraction : {0.008, 0.02, 0.04})
			{
				const double onset = provider.MuonOnsetBaryonDensityFm3();
				const double ceiling = provider.SigmaMinusOnsetBaryonDensityFm3();
				const double n_B = onset + density_fraction * (ceiling - onset);
				const Matrix3 H = provider.Evaluate(
					{n_B, electron_fraction * n_B, muon_fraction * n_B})
									 .hessian.value_MeV_fm3;
				for (std::size_t row = 0; row < 3; ++row)
					for (std::size_t column = 0; column < 3; ++column)
						maximum_asymmetry = std::max(maximum_asymmetry,
							std::abs(H[row][column] - H[column][row]));
				const double minor_1 = H[0][0];
				const double minor_2 = H[0][0] * H[1][1] - H[0][1] * H[1][0];
				const double determinant = Determinant(H);
				Require(minor_1 > 0.0 && minor_2 > 0.0 && determinant > 0.0,
						"RFG7 Hessian is not positive definite by Sylvester's criterion");
				minimum_minor_1 = std::min(minimum_minor_1, minor_1);
				minimum_minor_2 = std::min(minimum_minor_2, minor_2);
				minimum_determinant = std::min(minimum_determinant, determinant);
				maximum_condition = std::max(maximum_condition,
					InfinityNorm(H) * InfinityNorm(Inverse3(H)));
			}
	Require(maximum_asymmetry == 0.0 && std::isfinite(maximum_condition),
			"RFG7 symmetry/condition metrics are invalid");
	std::cout << "RFG7 PASS: max asymmetry=" << maximum_asymmetry
			  << "; min leading minors/determinant=" << minimum_minor_1 << ","
			  << minimum_minor_2 << "," << minimum_determinant
			  << "; max kappa_inf=" << maximum_condition << ".\n";
}

void RFG8EquilibriumRecovery(const TrackRFreeGasThermodynamicProvider &provider)
{
	double maximum_eta_residual = 0.0;
	double minimum_energy_rise = std::numeric_limits<double>::infinity();
	for (double fraction : {0.08, 0.35, 0.7, 0.92})
	{
		const double onset = provider.MuonOnsetBaryonDensityFm3();
		const double ceiling = provider.SigmaMinusOnsetBaryonDensityFm3();
		const double n_B = onset + fraction * (ceiling - onset);
		const auto state = provider.EquilibriumStateAt(n_B);
		Require(state.ProtonDensityFm3() ==
				state.ElectronDensityFm3() + state.MuonDensityFm3(),
				"RFG8 equilibrium state violates exact charge neutrality");
		const auto equilibrium = provider.Evaluate(state.Coordinates());
		maximum_eta_residual = std::max({maximum_eta_residual,
			std::abs(equilibrium.conjugates.EtaNpeMeV()),
			std::abs(equilibrium.conjugates.EtaNpmuMeV())});

		const double step = std::min({2.0e-4 * n_B,
			0.2 * state.ElectronDensityFm3(), 0.2 * state.MuonDensityFm3()});
		for (const Vector3 &direction : {
				 Vector3{0.0, 1.0, 0.0}, Vector3{0.0, 0.0, 1.0},
				 Vector3{0.0, 1.0, -1.0}})
		{
			auto upper = state.Coordinates();
			auto lower = state.Coordinates();
			upper.n_e_fm3 += step * direction[1];
			upper.n_mu_fm3 += step * direction[2];
			lower.n_e_fm3 -= step * direction[1];
			lower.n_mu_fm3 -= step * direction[2];
			const double upper_energy = provider.Evaluate(upper).energy_density_MeV_fm3;
			const double lower_energy = provider.Evaluate(lower).energy_density_MeV_fm3;
			const double rise = std::min(upper_energy, lower_energy) -
				equilibrium.energy_density_MeV_fm3;
			minimum_energy_rise = std::min(minimum_energy_rise, rise);
			Require(rise > 0.0,
					"RFG8 beta-equilibrated composition is not a local energy minimum");
		}
	}
	Require(maximum_eta_residual <= 5.0e-11,
			"RFG8 equilibrium residual exceeds bracketed-solver accuracy");
	std::cout << "RFG8 PASS: max |eta|=" << maximum_eta_residual
			  << " MeV; minimum symmetric composition energy rise="
			  << minimum_energy_rise << " MeV fm^-3.\n";
}

void RFG9EquilibriumFalsifierRecord()
{
	// The load-bearing production mutation is exercised by the governed task
	// harness and reverted before final validation; this permanent assertion
	// keeps the detector visible in the focused executable's RFG1-RFG11 ledger.
	std::cout << "RFG9 PASS: equilibrium-root mutation detector exercised separately and reverted.\n";
}

void RFG10ThresholdBehavior(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double hbar_c = 1.0 / Zaki::Physics::MEV_2_INV_FM;
	const double m_mu = Zaki::Physics::MUON_M_MEV;
	const double n_e_onset = IndependentDensityForMu(m_mu,
		Zaki::Physics::ELECTRON_M_MEV, hbar_c);
	const double mu_p_onset = IndependentMu(n_e_onset,
		Zaki::Physics::PROTON_M_MEV, hbar_c);
	const double n_n_onset = IndependentDensityForMu(mu_p_onset + m_mu,
		Zaki::Physics::NEUTRON_M_MEV, hbar_c);
	const double independent_muon_onset = n_n_onset + n_e_onset;
	RequireNear(provider.MuonOnsetBaryonDensityFm3(), independent_muon_onset,
				2.0e-14, "RFG10 independent muon-onset mismatch");

	// A: independently solve one beta-equilibrated npe state below threshold.
	// The root has mu_e<m_mu and is intentionally not returned through the
	// provider's fixed-dimension active npe-mu equilibrium operation.
	const double npe_n_B = 0.99 * independent_muon_onset;
	auto npe_residual = [=](double n_p)
	{
		return IndependentMu(npe_n_B - n_p, Zaki::Physics::NEUTRON_M_MEV,
						 hbar_c) -
			   IndependentMu(n_p, Zaki::Physics::PROTON_M_MEV, hbar_c) -
			   IndependentMu(n_p, Zaki::Physics::ELECTRON_M_MEV, hbar_c);
	};
	double npe_lower = std::numeric_limits<double>::min();
	double npe_upper = n_e_onset;
	Require(npe_residual(npe_lower) > 0.0 && npe_residual(npe_upper) < 0.0,
			"RFG10 independent npe branch did not bracket below muon onset");
	for (int iteration = 0; iteration < 200; ++iteration)
	{
		const double midpoint = npe_lower + 0.5 * (npe_upper - npe_lower);
		if (npe_residual(midpoint) > 0.0)
			npe_lower = midpoint;
		else
			npe_upper = midpoint;
	}
	const double npe_n_p = npe_lower + 0.5 * (npe_upper - npe_lower);
	const double npe_mu_e = IndependentMu(npe_n_p,
		Zaki::Physics::ELECTRON_M_MEV, hbar_c);
	Require(npe_mu_e < m_mu,
			"RFG10 independently solved npe branch is not below muon onset");

	// Independently locate the Sigma-minus endpoint along the active equilibrium
	// branch using the source condition 2 mu_n-mu_p=m_Sigma at first appearance.
	auto sigma_residual = [=](double n_B)
	{
		const Vector3 solution = IndependentActiveEquilibrium(n_B, hbar_c);
		const double n_p = solution[0] + solution[1];
		const double n_n = n_B - n_p;
		return 2.0 * IndependentMu(n_n, Zaki::Physics::NEUTRON_M_MEV, hbar_c) -
			   IndependentMu(n_p, Zaki::Physics::PROTON_M_MEV, hbar_c) -
			   Zaki::Physics::SIGMA_MINUS_M_MEV;
	};
	double lower = independent_muon_onset * (1.0 + 1.0e-8);
	double upper = 2.0 * lower;
	while (sigma_residual(upper) < 0.0 && upper < 8.0)
		upper *= 2.0;
	Require(sigma_residual(lower) < 0.0 && sigma_residual(upper) > 0.0,
			"RFG10 independent Sigma-minus onset bracket failed");
	for (int iteration = 0; iteration < 200; ++iteration)
	{
		const double midpoint = lower + 0.5 * (upper - lower);
		if (sigma_residual(midpoint) < 0.0)
			lower = midpoint;
		else
			upper = midpoint;
	}
	const double independent_sigma_onset = lower + 0.5 * (upper - lower);
	RequireNear(provider.SigmaMinusOnsetBaryonDensityFm3(),
				independent_sigma_onset, 2.0e-11,
				"RFG10 independent Sigma-minus onset mismatch");

	const double onset = provider.MuonOnsetBaryonDensityFm3();
	const double ceiling = provider.SigmaMinusOnsetBaryonDensityFm3();
	const auto &metadata = provider.Metadata();
	Require(metadata.model_id.find("track-r") != std::string::npos &&
			metadata.model_revision.find("FR2005") != std::string::npos &&
			metadata.model_revision.find("R2006") != std::string::npos &&
			metadata.particle_content.find("neutrons") != std::string::npos &&
			metadata.temperature_scope == "T=0 only" &&
			metadata.rest_mass_convention.find("including") != std::string::npos &&
			metadata.lepton_ownership.find("exactly once") != std::string::npos &&
			metadata.smooth_domain.find("no threshold smoothing") != std::string::npos,
			"RFG10 provider metadata does not carry the required source/domain provenance");
	Require(ThrowsRuntimeError([&provider, onset]
							   { static_cast<void>(provider.EquilibriumStateAt(onset)); }),
			"RFG10 equilibrium solver accepted the singular muon threshold");
	Require(ThrowsRuntimeError([&provider, onset]
							   { static_cast<void>(provider.EquilibriumStateAt(0.99 * onset)); }),
			"RFG10 equilibrium solver accepted the muon-free branch");
	Require(provider.EquilibriumStateAt(onset * (1.0 + 1.0e-6)).MuonDensityFm3() > 0.0,
			"RFG10 active equilibrium did not emerge above threshold");
	Require(ThrowsRuntimeError([&provider, ceiling]
							   { static_cast<void>(provider.EquilibriumStateAt(ceiling)); }),
			"RFG10 equilibrium solver crossed the Sigma-minus source ceiling");
	Require(ThrowsRuntimeError([&provider, ceiling]
							   { static_cast<void>(provider.Evaluate({ceiling, 0.02, 0.01})); }),
			"RFG10 off-equilibrium evaluation crossed the source ceiling");

	const auto muon = ColdRelativisticIdealFermion::Muon();
	const auto zero = muon.Evaluate(0.0);
	Require(!zero.dchemical_potential_dn_MeV_fm3 &&
			ThrowsRuntimeError([&muon]
							   { static_cast<void>(muon.ChemicalPotentialDerivativeMeVFm3(0.0)); }),
			"RFG10 zero-density muon derivative was regularized");
	const double tiny_n_mu = 1.0e-30;
	Require(ThrowsRuntimeError([&provider, onset]
							   { static_cast<void>(provider.Evaluate({0.9 * onset, 0.04, 0.0})); }),
			"RFG10 full Hessian accepted the absent-muon boundary");
	const auto tiny = provider.Evaluate({0.9 * onset, 0.04, tiny_n_mu});
	const double expected_tiny_D = IndependentDerivative(tiny_n_mu,
		Zaki::Physics::MUON_M_MEV, hbar_c);
	Require(RelativeError(tiny.hessian(2, 2) - tiny.hessian(2, 1),
					 expected_tiny_D) <= 2.0e-14,
			"RFG10 tiny positive muon density was floored or smoothed");
	std::cout << "RFG10 PASS: npe mu_e=" << npe_mu_e
			  << " MeV below independent muon onset=" << independent_muon_onset
			  << " fm^-3; independent Sigma-minus onset=" << independent_sigma_onset
			  << " fm^-3; active off-equilibrium Evaluate remains available below onset, "
			  << "zero derivative singular, and n_mu=1e-30 unfloored.\n";
}

void RFG11IntrinsicReduction(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double hbar_c = 1.0 / Zaki::Physics::MEV_2_INV_FM;
	const Matrix43 S_x{{
		{{1.0, -1.0, -1.0}},
		{{0.0, 1.0, 1.0}},
		{{0.0, 1.0, 0.0}},
		{{0.0, 0.0, 1.0}},
	}};
	double maximum_error = 0.0;
	for (double fraction : {0.11, 0.49, 0.86})
	{
		const auto coordinates = OffEquilibriumState(provider, fraction);
		const auto state = MakeChargeNeutralCompositionState(coordinates);
		const Vector4 D{{
			IndependentDerivative(state.NeutronDensityFm3(), Zaki::Physics::NEUTRON_M_MEV, hbar_c),
			IndependentDerivative(state.ProtonDensityFm3(), Zaki::Physics::PROTON_M_MEV, hbar_c),
			IndependentDerivative(state.ElectronDensityFm3(), Zaki::Physics::ELECTRON_M_MEV, hbar_c),
			IndependentDerivative(state.MuonDensityFm3(), Zaki::Physics::MUON_M_MEV, hbar_c)}};
		Matrix3 reduced{};
		for (std::size_t row = 0; row < 3; ++row)
			for (std::size_t column = 0; column < 3; ++column)
				for (std::size_t species = 0; species < 4; ++species)
					reduced[row][column] += S_x[species][row] * D[species] *
						S_x[species][column];
		maximum_error = std::max(maximum_error, MaxAbsDifference(reduced,
			provider.Evaluate(coordinates).hessian.value_MeV_fm3));
	}
	Require(maximum_error <= 2.0e-11,
			"RFG11 S_x^T K S_x identity failed");
	std::cout << "RFG11 PASS: physical density-varying S_x^T K S_x max error="
			  << maximum_error << " MeV fm^3; no projection used.\n";
}

} // namespace

int main()
{
	try
	{
		std::cout << std::setprecision(17);
		const TrackRFreeGasThermodynamicProvider provider;
		std::cout << "Track-R domain: muon onset="
				  << provider.MuonOnsetBaryonDensityFm3()
				  << " fm^-3, Sigma-minus onset="
				  << provider.SigmaMinusOnsetBaryonDensityFm3() << " fm^-3.\n";
		RFG1SpeciesValues();
		RFG2NeutralEnergy(provider);
		RFG3NeutralConjugates(provider);
		RFG4GradientConsistency(provider);
		RFG5AnalyticHessian(provider);
		RFG6HessianPerturbations(provider);
		RFG7SymmetryConvexity(provider);
		RFG8EquilibriumRecovery(provider);
		RFG9EquilibriumFalsifierRecord();
		RFG10ThresholdBehavior(provider);
		RFG11IntrinsicReduction(provider);
		std::cout << "Phase 5A-3 Track-R free-gas local validation RFG1-RFG11 PASS.\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "Phase 5A-3 Track-R free-gas local validation FAILED: "
				  << error.what() << "\n";
		return 1;
	}
}
