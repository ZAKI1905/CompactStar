#include "CompactStar/EOS/LocalThermodynamics.hpp"

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

namespace
{

using CompactStar::ChargeNeutralChemicalHessian;
using CompactStar::ChargeNeutralCompositionState;
using CompactStar::ChargeNeutralCoordinates;
using CompactStar::ColdRelativisticFreeLepton;
using CompactStar::ILocalThermodynamicProvider;
using CompactStar::LocalThermodynamicEvaluation;
using CompactStar::LocalThermodynamicProviderMetadata;
using CompactStar::MakeChargeNeutralCompositionState;
using CompactStar::NeutralConjugates;

using Vector3 = std::array<double, 3>;
using Vector4 = std::array<double, 4>;
using Matrix3 = std::array<std::array<double, 3>, 3>;
using Matrix4 = std::array<std::array<double, 4>, 4>;
using Matrix43 = std::array<std::array<double, 3>, 4>;

constexpr double kEps0 = 200.0;
constexpr Vector3 kReferenceX{0.20, 0.04, 0.02};
constexpr Vector3 kReferenceG{950.0, 0.0, 0.0};
constexpr double kCubicCoefficient = 10.0;
constexpr Matrix3 kReferenceH{{
	{{2.0, -2.0, -2.0}},
	{{-2.0, 10.0, 5.0}},
	{{-2.0, 5.0, 12.0}},
}};

void Require(bool condition, const std::string &message)
{
	if (!condition)
		throw std::runtime_error(message);
}

void RequireNear(double actual, double expected, double absolute_tolerance,
				 const std::string &message)
{
	if (!std::isfinite(actual) || !std::isfinite(expected) ||
		std::abs(actual - expected) > absolute_tolerance)
	{
		throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
								 ", expected=" + std::to_string(expected) +
								 ", tolerance=" + std::to_string(absolute_tolerance));
	}
}

double ScaledTolerance(double scale, double factor = 256.0)
{
	return factor * std::numeric_limits<double>::epsilon() *
		   std::max(1.0, std::abs(scale));
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

Vector3 Add(const Vector3 &a, const Vector3 &b)
{
	Vector3 result{};
	for (std::size_t i = 0; i < 3; ++i)
		result[i] = a[i] + b[i];
	return result;
}

Vector3 Subtract(const Vector3 &a, const Vector3 &b)
{
	Vector3 result{};
	for (std::size_t i = 0; i < 3; ++i)
		result[i] = a[i] - b[i];
	return result;
}

Vector3 Scale(const Vector3 &a, double factor)
{
	Vector3 result{};
	for (std::size_t i = 0; i < 3; ++i)
		result[i] = factor * a[i];
	return result;
}

Vector3 Multiply(const Matrix3 &matrix, const Vector3 &vector)
{
	Vector3 result{};
	for (std::size_t row = 0; row < 3; ++row)
		for (std::size_t column = 0; column < 3; ++column)
			result[row] += matrix[row][column] * vector[column];
	return result;
}

Matrix3 Multiply(const Matrix3 &left, const Matrix3 &right)
{
	Matrix3 result{};
	for (std::size_t row = 0; row < 3; ++row)
		for (std::size_t column = 0; column < 3; ++column)
			for (std::size_t inner = 0; inner < 3; ++inner)
				result[row][column] += left[row][inner] * right[inner][column];
	return result;
}

Matrix3 Transpose(const Matrix3 &matrix)
{
	Matrix3 result{};
	for (std::size_t row = 0; row < 3; ++row)
		for (std::size_t column = 0; column < 3; ++column)
			result[column][row] = matrix[row][column];
	return result;
}

double Determinant(const Matrix3 &m)
{
	return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
		   m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
		   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Transparent test-side inverse.  No production provider inverse is used as an oracle.
Matrix3 Inverse3(const Matrix3 &m)
{
	const double determinant = Determinant(m);
	if (!std::isfinite(determinant) || std::abs(determinant) < 1.0e-14)
		throw std::runtime_error("Test fixture requested inversion of a singular 3x3 matrix.");

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

double InfinityNorm(const Matrix3 &matrix)
{
	double maximum = 0.0;
	for (const auto &row : matrix)
	{
		double sum = 0.0;
		for (double value : row)
			sum += std::abs(value);
		maximum = std::max(maximum, sum);
	}
	return maximum;
}

double MaxAbsDifference(const Matrix3 &left, const Matrix3 &right)
{
	double maximum = 0.0;
	for (std::size_t row = 0; row < 3; ++row)
		for (std::size_t column = 0; column < 3; ++column)
			maximum = std::max(maximum, std::abs(left[row][column] - right[row][column]));
	return maximum;
}

double MaxAbs(const Vector3 &vector)
{
	double maximum = 0.0;
	for (double value : vector)
		maximum = std::max(maximum, std::abs(value));
	return maximum;
}

Vector3 CoordinatesAsVector(const ChargeNeutralCoordinates &coordinates)
{
	return {coordinates.n_B_fm3, coordinates.n_e_fm3, coordinates.n_mu_fm3};
}

ChargeNeutralCoordinates VectorAsCoordinates(const Vector3 &vector)
{
	return {vector[0], vector[1], vector[2]};
}

double ToyEnergyClosedForm(const Vector3 &x)
{
	const Vector3 displacement = Subtract(x, kReferenceX);
	const Vector3 linear_response = Multiply(kReferenceH, displacement);
	double quadratic = 0.0;
	for (std::size_t i = 0; i < 3; ++i)
		quadratic += displacement[i] * linear_response[i];
	return kEps0 + kReferenceG[0] * displacement[0] +
		   0.5 * quadratic + kCubicCoefficient * displacement[0] * displacement[0] * displacement[0] / 6.0;
}

Vector3 ToyGradientClosedForm(const Vector3 &x)
{
	const Vector3 displacement = Subtract(x, kReferenceX);
	Vector3 gradient = Add(kReferenceG, Multiply(kReferenceH, displacement));
	gradient[0] += 0.5 * kCubicCoefficient * displacement[0] * displacement[0];
	return gradient;
}

Matrix3 ToyHessianClosedForm(const Vector3 &x)
{
	Matrix3 hessian = kReferenceH;
	hessian[0][0] += kCubicCoefficient * (x[0] - kReferenceX[0]);
	return hessian;
}

class AnalyticChargeNeutralToyProvider final : public ILocalThermodynamicProvider
{
  public:
	const LocalThermodynamicProviderMetadata &Metadata() const noexcept override
	{
		return metadata_;
	}

	LocalThermodynamicEvaluation
	Evaluate(const ChargeNeutralCoordinates &coordinates) const override
	{
		const auto state = MakeChargeNeutralCompositionState(coordinates);
		RequireSmoothDomain(state);
		const double d_B = coordinates.n_B_fm3 - kReferenceX[0];
		const double d_e = coordinates.n_e_fm3 - kReferenceX[1];
		const double d_mu = coordinates.n_mu_fm3 - kReferenceX[2];
		const double energy =
			kEps0 + kReferenceG[0] * d_B +
			0.5 * (2.0 * d_B * d_B + 10.0 * d_e * d_e +
				   12.0 * d_mu * d_mu - 4.0 * d_B * d_e -
				   4.0 * d_B * d_mu + 10.0 * d_e * d_mu) +
			kCubicCoefficient * d_B * d_B * d_B / 6.0;
		const Vector3 conjugates{
			kReferenceG[0] + 2.0 * d_B - 2.0 * d_e - 2.0 * d_mu +
				0.5 * kCubicCoefficient * d_B * d_B,
			-2.0 * d_B + 10.0 * d_e + 5.0 * d_mu,
			-2.0 * d_B + 5.0 * d_e + 12.0 * d_mu};
		const Matrix3 hessian{{
			{{2.0 + kCubicCoefficient * d_B, -2.0, -2.0}},
			{{-2.0, 10.0, 5.0}},
			{{-2.0, 5.0, 12.0}},
		}};
		return {state, energy, NeutralConjugates{conjugates},
				ChargeNeutralChemicalHessian{hessian}};
	}

	ChargeNeutralCompositionState EquilibriumStateAt(double n_B_fm3) const override
	{
		if (!std::isfinite(n_B_fm3) || n_B_fm3 < 0.15 || n_B_fm3 > 0.25)
			throw std::runtime_error("Toy equilibrium request is outside its declared domain.");
		const double delta_B = n_B_fm3 - kReferenceX[0];
		// Solve H_cc delta_c = -H_cB delta_B exactly; det(H_cc)=95.
		const double n_e = kReferenceX[1] + (14.0 / 95.0) * delta_B;
		const double n_mu = kReferenceX[2] + (10.0 / 95.0) * delta_B;
		const auto state = MakeChargeNeutralCompositionState({n_B_fm3, n_e, n_mu});
		RequireSmoothDomain(state);
		return state;
	}

	// Optional, concrete-model split used only to test S_x^T mu=g.  It is not
	// part of ILocalThermodynamicProvider.
	std::array<double, 4>
	IntrinsicPotentials(const ChargeNeutralCoordinates &coordinates) const
	{
		const auto evaluation = Evaluate(coordinates);
		const auto &g = evaluation.conjugates.value_MeV;
		constexpr double mu_p = 100.0;
		return {g[0], mu_p, g[1] + g[0] - mu_p, g[2] + g[0] - mu_p};
	}

  private:
	static void RequireSmoothDomain(const ChargeNeutralCompositionState &state)
	{
		if (state.BaryonDensityFm3() < 0.15 || state.BaryonDensityFm3() > 0.25 ||
			state.ElectronDensityFm3() <= 0.0 || state.MuonDensityFm3() <= 0.0 ||
			state.NeutronDensityFm3() <= 0.0)
			throw std::runtime_error("Toy state is outside its smooth, active-species domain.");
	}

	const LocalThermodynamicProviderMetadata metadata_{
		"phase5a2-charge-neutral-analytic-toy",
		"1",
		"npe-mu test particles only",
		"x=(n_B,n_e,n_mu), all fm^-3",
		"cold-only, T=0",
		"arbitrary test energy includes its declared linear rest-energy-like term",
		"none; this is TEST PHYSICS, not a neutron-star EOS",
		"0.15 <= n_B <= 0.25 fm^-3 with n_n,n_e,n_mu strictly positive"};
};

struct IndependentLeptonValues
{
	double p_F;
	double mu;
	double energy;
	double derivative;
};

IndependentLeptonValues IndependentFreeLepton(double density, double mass,
											  double hbar_c)
{
	const long double pi = std::acos(-1.0L);
	const long double n = density;
	const long double m = mass;
	const long double hc = hbar_c;
	const long double p = hc * std::pow(3.0L * pi * pi * n, 1.0L / 3.0L);
	const long double mu = std::sqrt(m * m + p * p);
	const long double energy =
		(p * mu * (2.0L * p * p + m * m) -
		 m * m * m * m * std::log((p + mu) / m)) /
		(8.0L * pi * pi * hc * hc * hc);
	const long double derivative = p * p / (3.0L * mu * n);
	return {static_cast<double>(p), static_cast<double>(mu),
			static_cast<double>(energy), static_cast<double>(derivative)};
}

void V1UnitsAndSigns(const AnalyticChargeNeutralToyProvider &provider)
{
	const ChargeNeutralCoordinates coordinates{0.21, 0.043, 0.021};
	const auto evaluation = provider.Evaluate(coordinates);
	const auto potentials = provider.IntrinsicPotentials(coordinates);
	const double eta_npe = potentials[0] - potentials[1] - potentials[2];
	const double eta_npmu = potentials[0] - potentials[1] - potentials[3];
	RequireNear(evaluation.conjugates.MuNMeV(), potentials[0], 0.0,
				"V1 g_0 must be mu_n");
	RequireNear(evaluation.conjugates.value_MeV[1], -eta_npe,
				ScaledTolerance(eta_npe), "V1 g_1 sign is not -eta_npe");
	RequireNear(evaluation.conjugates.value_MeV[2], -eta_npmu,
				ScaledTolerance(eta_npmu), "V1 g_2 sign is not -eta_npmu");
	Require(provider.Metadata().coordinate_chart.find("fm^-3") != std::string::npos,
			"V1 provider coordinate units are not declared");
	std::cout << "V1 PASS: x in fm^-3, g in MeV, H in MeV fm^3; eta signs verified.\n";
}

void V2ChargeNeutrality()
{
	for (const ChargeNeutralCoordinates coordinates : {
			 ChargeNeutralCoordinates{0.16, 0.02, 0.01},
			 ChargeNeutralCoordinates{0.20, 0.04, 0.02},
			 ChargeNeutralCoordinates{0.24, 0.05, 0.03}})
	{
		const auto state = MakeChargeNeutralCompositionState(coordinates);
		Require(state.ProtonDensityFm3() ==
					state.ElectronDensityFm3() + state.MuonDensityFm3(),
				"V2 charge neutrality is not exact to arithmetic precision");
		Require(state.NeutronDensityFm3() ==
					state.BaryonDensityFm3() - state.ProtonDensityFm3(),
				"V2 baryon-density reconstruction is not exact");
	}
	std::cout << "V2 PASS: n_p=n_e+n_mu and n_n=n_B-n_p exactly by construction.\n";
}

void V3BetaEquilibrium(const AnalyticChargeNeutralToyProvider &provider)
{
	double largest_residual = 0.0;
	for (double n_B : {0.16, 0.20, 0.24})
	{
		const auto state = provider.EquilibriumStateAt(n_B);
		const auto evaluation = provider.Evaluate(state.Coordinates());
		largest_residual = std::max(largest_residual,
									std::max(std::abs(evaluation.conjugates.value_MeV[1]),
											 std::abs(evaluation.conjugates.value_MeV[2])));
	}
	Require(largest_residual <= 4.0e-16, "V3 toy beta-equilibrium residual is too large");
	std::cout << "V3 PASS: max |composition conjugate|=" << largest_residual
			  << " MeV at constructed beta equilibrium.\n";
}

void V4AnalyticFreeLeptons()
{
	const auto electron = ColdRelativisticFreeLepton::Electron();
	const auto muon = ColdRelativisticFreeLepton::Muon();
	RequireNear(electron.RestMassEnergyMeV(), 0.51099895, 2.0e-8,
				"V4 electron rest-mass convention drifted");
	RequireNear(muon.RestMassEnergyMeV(), 105.6583755, 2.0e-7,
				"V4 muon rest-mass convention drifted");
	RequireNear(electron.HbarCMeVFm(), 197.3269804, 2.0e-6,
				"V4 hbar-c convention drifted");

	double maximum_relative_error = 0.0;
	for (const auto lepton : {electron, muon})
	{
		for (double density : {1.0e-5, 1.0e-3, 1.0e-2, 1.0e-1})
		{
			const auto actual = lepton.Evaluate(density);
			const auto expected = IndependentFreeLepton(
				density, lepton.RestMassEnergyMeV(), lepton.HbarCMeVFm());
			for (const auto pair : {
					 std::array<double, 2>{actual.fermi_momentum_MeV, expected.p_F},
					 std::array<double, 2>{actual.chemical_potential_MeV, expected.mu},
					 std::array<double, 2>{*actual.dchemical_potential_dn_MeV_fm3,
										   expected.derivative}})
			{
				const double relative = std::abs(pair[0] - pair[1]) /
										std::max(1.0, std::abs(pair[1]));
				maximum_relative_error = std::max(maximum_relative_error, relative);
				Require(relative <= 2.0e-14, "V4 analytic lepton value/derivative mismatch");
			}
			const double energy_relative = std::abs(actual.energy_density_MeV_fm3 -
													expected.energy) /
										   std::max(1.0, std::abs(expected.energy));
			maximum_relative_error = std::max(maximum_relative_error, energy_relative);
			Require(energy_relative <= 3.0e-11,
					"V4 analytic lepton energy-density mismatch");
		}
	}
	std::cout << "V4 PASS: independent free-lepton max scaled error="
			  << maximum_relative_error << ".\n";
}

void V5AnalyticToy(const AnalyticChargeNeutralToyProvider &provider)
{
	const Vector3 x{0.213, 0.046, 0.024};
	const auto evaluation = provider.Evaluate(VectorAsCoordinates(x));
	const auto expected_g = ToyGradientClosedForm(x);
	const auto expected_h = ToyHessianClosedForm(x);
	RequireNear(evaluation.energy_density_MeV_fm3, ToyEnergyClosedForm(x),
				ScaledTolerance(ToyEnergyClosedForm(x)),
				"V5 toy energy differs from independent closed form");
	for (std::size_t row = 0; row < 3; ++row)
	{
		RequireNear(evaluation.conjugates.value_MeV[row], expected_g[row], 0.0,
					"V5 toy gradient differs from closed form");
		for (std::size_t column = 0; column < 3; ++column)
			RequireNear(evaluation.hessian(row, column), expected_h[row][column], 0.0,
						"V5 toy Hessian differs from closed form");
	}
	std::cout << "V5 PASS: analytic toy energy, gradient, and Hessian are exact.\n";
}

void V6HessianSymmetry(const AnalyticChargeNeutralToyProvider &provider)
{
	const auto evaluation = provider.Evaluate({0.217, 0.047, 0.025});
	double maximum_asymmetry = 0.0;
	for (std::size_t row = 0; row < 3; ++row)
		for (std::size_t column = 0; column < 3; ++column)
			maximum_asymmetry = std::max(maximum_asymmetry,
										 std::abs(evaluation.hessian(row, column) -
												  evaluation.hessian(column, row)));
	Require(maximum_asymmetry == 0.0, "V6 analytic Hessian violates Maxwell symmetry");
	std::cout << "V6 PASS: max |H_ab-H_ba|=" << maximum_asymmetry << ".\n";
}

void V7LinearResponseConvergence(const AnalyticChargeNeutralToyProvider &provider)
{
	const Vector3 base{0.20, 0.04, 0.02};
	const Vector3 direction{1.0, 0.3, -0.2};
	const auto base_evaluation = provider.Evaluate(VectorAsCoordinates(base));
	const Matrix3 hessian = base_evaluation.hessian.value_MeV_fm3;
	std::array<double, 4> errors{};
	const std::array<double, 4> steps{1.0e-2, 5.0e-3, 2.5e-3, 1.25e-3};
	for (std::size_t index = 0; index < steps.size(); ++index)
	{
		const Vector3 delta = Scale(direction, steps[index]);
		const auto perturbed = provider.Evaluate(VectorAsCoordinates(Add(base, delta)));
		const Vector3 actual = Subtract(perturbed.conjugates.value_MeV,
										base_evaluation.conjugates.value_MeV);
		const Vector3 predicted = Multiply(hessian, delta);
		errors[index] = MaxAbs(Subtract(actual, predicted));
	}
	double minimum_order = std::numeric_limits<double>::infinity();
	for (std::size_t index = 0; index + 1 < errors.size(); ++index)
	{
		Require(errors[index + 1] < errors[index], "V7 linear-response error did not decrease");
		const double order = std::log(errors[index] / errors[index + 1]) / std::log(2.0);
		minimum_order = std::min(minimum_order, order);
	}
	Require(minimum_order > 1.99, "V7 did not demonstrate the expected quadratic remainder");
	std::cout << "V7 PASS: errors=" << errors[0] << "," << errors[1] << ","
			  << errors[2] << "," << errors[3]
			  << "; minimum observed order=" << minimum_order << ".\n";
}

void V8BasisEquivalence(const AnalyticChargeNeutralToyProvider &provider)
{
	const Matrix3 T{{{{1.0, -1.0, -1.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
	const Matrix3 U{{{{1.0, 1.0, 1.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
	const auto evaluation = provider.Evaluate({0.218, 0.049, 0.026});
	const Matrix3 H_x = evaluation.hessian.value_MeV_fm3;
	const Matrix3 H_y = Multiply(Transpose(U), Multiply(H_x, U));
	const Matrix3 reconstructed_H_x = Multiply(Transpose(T), Multiply(H_y, T));
	const double h_error = MaxAbsDifference(H_x, reconstructed_H_x);
	Require(h_error <= ScaledTolerance(InfinityNorm(H_x)),
			"V8 exact x/y Hessian transformation failed");

	const Matrix3 C_x = Inverse3(H_x);
	const Matrix3 transformed_C_y = Multiply(T, Multiply(C_x, Transpose(T)));
	const Matrix3 direct_C_y = Inverse3(H_y);
	const double inverse_error = MaxAbsDifference(transformed_C_y, direct_C_y);
	Require(inverse_error <= ScaledTolerance(InfinityNorm(direct_C_y)),
			"V8 inverse responses do not transform consistently");

	// z=(n_B,Y_e,Y_mu): x=(z_B,z_B z_e,z_B z_mu).  Away from
	// equilibrium, Hessian transformation includes g_i d2 x_i/dz_a dz_b.
	const auto &g = evaluation.conjugates.value_MeV;
	const double n_B = evaluation.state.BaryonDensityFm3();
	const double Y_e = evaluation.state.ElectronDensityFm3() / n_B;
	const double Y_mu = evaluation.state.MuonDensityFm3() / n_B;
	const Matrix3 J{{{{1.0, 0.0, 0.0}}, {{Y_e, n_B, 0.0}}, {{Y_mu, 0.0, n_B}}}};
	const Matrix3 naive = Multiply(Transpose(J), Multiply(H_x, J));
	Matrix3 full = naive;
	full[0][1] += g[1];
	full[1][0] += g[1];
	full[0][2] += g[2];
	full[2][0] += g[2];
	Require(std::abs(g[1]) > 1.0e-3 && std::abs(g[2]) > 1.0e-3,
			"V8 nonlinear chain-rule fixture accidentally lies at equilibrium");
	RequireNear(full[0][1] - naive[0][1], g[1], ScaledTolerance(g[1]),
				"V8 electron fraction-coordinate chain term is absent");
	RequireNear(full[0][2] - naive[0][2], g[2], ScaledTolerance(g[2]),
				"V8 muon fraction-coordinate chain term is absent");
	std::cout << "V8 PASS: max linear-basis H error=" << h_error
			  << ", inverse error=" << inverse_error
			  << "; nonlinear chain terms are nonzero.\n";
}

Matrix4 CorrectedChargeProjectedSusceptibility()
{
	// Independent full-intrinsic fixture: K=diag(2,3,5,7) MeV fm^3,
	// chi=K^-1. Species order is (n,p,e,mu), q=(0,1,-1,-1).
	const Vector4 chi_diagonal{0.5, 1.0 / 3.0, 0.2, 1.0 / 7.0};
	const Vector4 q{0.0, 1.0, -1.0, -1.0};
	Vector4 chi_q{};
	double denominator = 0.0;
	for (std::size_t i = 0; i < 4; ++i)
	{
		chi_q[i] = chi_diagonal[i] * q[i];
		denominator += q[i] * chi_q[i];
	}
	Matrix4 response{};
	for (std::size_t row = 0; row < 4; ++row)
		for (std::size_t column = 0; column < 4; ++column)
			response[row][column] =
				(row == column ? chi_diagonal[row] : 0.0) -
				chi_q[row] * chi_q[column] / denominator;
	return response;
}

Matrix4 ReducedResponseFromProviderHessian(const Matrix3 &hessian)
{
	const Matrix43 S_x{{
		{{1.0, -1.0, -1.0}},
		{{0.0, 1.0, 1.0}},
		{{0.0, 1.0, 0.0}},
		{{0.0, 0.0, 1.0}},
	}};
	const Matrix3 inverse = Inverse3(hessian);
	Matrix4 response{};
	for (std::size_t row = 0; row < 4; ++row)
		for (std::size_t column = 0; column < 4; ++column)
			for (std::size_t a = 0; a < 3; ++a)
				for (std::size_t b = 0; b < 3; ++b)
					response[row][column] += S_x[row][a] * inverse[a][b] * S_x[column][b];
	return response;
}

void V9CorrectedProjection(const AnalyticChargeNeutralToyProvider &provider)
{
	const auto provider_evaluation = provider.Evaluate(VectorAsCoordinates(kReferenceX));
	const Matrix4 projected = CorrectedChargeProjectedSusceptibility();
	const Matrix4 reduced = ReducedResponseFromProviderHessian(
		provider_evaluation.hessian.value_MeV_fm3);
	double maximum_equivalence_error = 0.0;
	for (std::size_t row = 0; row < 4; ++row)
		for (std::size_t column = 0; column < 4; ++column)
			maximum_equivalence_error = std::max(maximum_equivalence_error,
												 std::abs(projected[row][column] - reduced[row][column]));
	Require(maximum_equivalence_error <= 2.0e-16,
			"V9 corrected full-intrinsic projection disagrees with constrained Hessian");

	const Vector4 q{0.0, 1.0, -1.0, -1.0};
	double maximum_charge_null = 0.0;
	for (std::size_t row = 0; row < 4; ++row)
	{
		double value = 0.0;
		for (std::size_t column = 0; column < 4; ++column)
			value += projected[row][column] * q[column];
		maximum_charge_null = std::max(maximum_charge_null, std::abs(value));
	}
	double maximum_proton_identity = 0.0;
	for (std::size_t column = 0; column < 4; ++column)
		maximum_proton_identity = std::max(maximum_proton_identity,
										   std::abs(projected[1][column] - projected[2][column] -
													projected[3][column]));
	Require(maximum_charge_null <= 1.0e-16, "V9 projected response lacks charge null");
	Require(maximum_proton_identity <= 1.0e-16,
			"V9 proton response is not electron-plus-muon response");
	std::cout << "V9 PASS: projection equivalence error=" << maximum_equivalence_error
			  << ", charge-null residual=" << maximum_charge_null
			  << ", proton-row residual=" << maximum_proton_identity << ".\n";
}

void V10DomainRankAndCondition(const AnalyticChargeNeutralToyProvider &provider)
{
	double minimum_minor_1 = std::numeric_limits<double>::infinity();
	double minimum_minor_2 = std::numeric_limits<double>::infinity();
	double minimum_determinant = std::numeric_limits<double>::infinity();
	double maximum_condition_number = 0.0;
	for (double n_B : {0.15, 0.20, 0.25})
	{
		const Matrix3 hessian = provider.Evaluate({n_B, 0.04, 0.02})
									.hessian.value_MeV_fm3;
		const double leading_minor_1 = hessian[0][0];
		const double leading_minor_2 = hessian[0][0] * hessian[1][1] -
									   hessian[0][1] * hessian[1][0];
		const double determinant = Determinant(hessian);
		Require(leading_minor_1 > 0.0 && leading_minor_2 > 0.0 && determinant > 0.0,
				"V10 toy Hessian is not positive definite in its declared region");
		minimum_minor_1 = std::min(minimum_minor_1, leading_minor_1);
		minimum_minor_2 = std::min(minimum_minor_2, leading_minor_2);
		minimum_determinant = std::min(minimum_determinant, determinant);
		const Matrix3 inverse = Inverse3(hessian);
		maximum_condition_number = std::max(
			maximum_condition_number, InfinityNorm(hessian) * InfinityNorm(inverse));
	}
	Require(std::isfinite(maximum_condition_number) && maximum_condition_number < 100.0,
			"V10 toy response condition metric is not sensible");

	const double nan = std::numeric_limits<double>::quiet_NaN();
	const double infinity = std::numeric_limits<double>::infinity();
	Require(ThrowsRuntimeError([]
							   { static_cast<void>(MakeChargeNeutralCompositionState({0.0, 0.0, 0.0})); }),
			"V10 zero baryon density did not fail closed");
	Require(ThrowsRuntimeError([]
							   { static_cast<void>(MakeChargeNeutralCompositionState({0.2, -0.01, 0.01})); }),
			"V10 negative electron density did not fail closed");
	Require(ThrowsRuntimeError([]
							   { static_cast<void>(MakeChargeNeutralCompositionState({0.2, 0.01, -0.01})); }),
			"V10 negative muon density did not fail closed");
	Require(ThrowsRuntimeError([]
							   { static_cast<void>(MakeChargeNeutralCompositionState({0.2, 0.15, 0.06})); }),
			"V10 negative reconstructed neutron density did not fail closed");
	Require(ThrowsRuntimeError([nan]
							   { static_cast<void>(MakeChargeNeutralCompositionState({nan, 0.01, 0.01})); }),
			"V10 NaN state did not fail closed");
	Require(ThrowsRuntimeError([infinity]
							   { static_cast<void>(MakeChargeNeutralCompositionState({0.2, infinity, 0.01})); }),
			"V10 infinite state did not fail closed");
	Require(ThrowsRuntimeError([&provider]
							   { provider.Evaluate({0.20, 0.04, 0.0}); }),
			"V10 absent-muon toy derivative request was silently regularized");

	const auto muon = ColdRelativisticFreeLepton::Muon();
	const auto zero_muon = muon.Evaluate(0.0);
	Require(zero_muon.chemical_potential_MeV == muon.RestMassEnergyMeV() &&
				!zero_muon.dchemical_potential_dn_MeV_fm3,
			"V10 zero-density muon value/derivative boundary is incorrect");
	Require(ThrowsRuntimeError([&muon]
							   { static_cast<void>(muon.ChemicalPotentialDerivativeMeVFm3(0.0)); }),
			"V10 zero-density muon derivative did not fail closed");
	Require(ThrowsRuntimeError([&muon, nan]
							   { static_cast<void>(muon.Evaluate(nan)); }),
			"V10 free-lepton NaN did not fail closed");
	Require(ThrowsRuntimeError([&muon, infinity]
							   { static_cast<void>(muon.Evaluate(infinity)); }),
			"V10 free-lepton infinity did not fail closed");
	const Matrix3 singular{{{{1.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
	Require(ThrowsRuntimeError([&singular]
							   { Inverse3(singular); }),
			"V10 singular test response did not fail closed when inversion was requested");
	std::cout << "V10 PASS: minimum leading minors/determinant=" << minimum_minor_1
			  << "," << minimum_minor_2 << "," << minimum_determinant
			  << "; max kappa_inf=" << maximum_condition_number
			  << "; invalid/singular/threshold cases rejected.\n";
}

} // namespace

int main()
{
	try
	{
		std::cout << std::setprecision(17);
		const AnalyticChargeNeutralToyProvider provider;
		V1UnitsAndSigns(provider);
		V2ChargeNeutrality();
		V3BetaEquilibrium(provider);
		V4AnalyticFreeLeptons();
		V5AnalyticToy(provider);
		V6HessianSymmetry(provider);
		V7LinearResponseConvergence(provider);
		V8BasisEquivalence(provider);
		V9CorrectedProjection(provider);
		V10DomainRankAndCondition(provider);
		std::cout << "Phase 5A-2 local thermodynamic validation V1-V10 PASS.\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "Phase 5A-2 local thermodynamic validation FAILED: "
				  << error.what() << "\n";
		return 1;
	}
}
