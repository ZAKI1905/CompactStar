#include "CompactStar/EOS/TrackRFreeGasThermodynamics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include <Zaki/Physics/Constants.hpp>

namespace
{
using namespace CompactStar;

const double mn = Zaki::Physics::NEUTRON_M_MEV;
const double mp = Zaki::Physics::PROTON_M_MEV;
const double me = Zaki::Physics::ELECTRON_M_MEV;
const double hc = 1.0 / Zaki::Physics::MEV_2_INV_FM;
const double pi = std::acos(-1.0);

void Require(bool ok, const std::string &message)
{
	if (!ok)
		throw std::runtime_error(message);
}

void Near(double actual, double expected, double tolerance, const std::string &message)
{
	if (!std::isfinite(actual) || !std::isfinite(expected) ||
		std::abs(actual - expected) > tolerance)
		throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
			", expected=" + std::to_string(expected));
}

template <class Error = std::runtime_error, class Function>
bool Throws(Function &&function)
{
	try
	{
		function();
	}
	catch (const Error &)
	{
		return true;
	}
	return false;
}

// Independent test-side ideal-gas expressions. No production fermion helper is
// used by the energy, chemical-potential, derivative, pressure, or onset oracles.
double Momentum(double density)
{
	return hc * std::cbrt(3.0 * pi * pi * density);
}

double Mu(double density, double mass)
{
	const double momentum = Momentum(density);
	return std::sqrt(mass * mass + momentum * momentum);
}

double Derivative(double density, double mass)
{
	const double momentum = Momentum(density);
	return pi * pi * hc * hc * hc /
		(std::sqrt(mass * mass + momentum * momentum) * momentum);
}

double Energy(double density, double mass)
{
	if (density == 0.0)
		return 0.0;
	// Direct phase-space quadrature in u=p/p_F, g=2. The rest mass is included.
	constexpr int panels = 16384;
	const double momentum = Momentum(density);
	double sum = 0.0;
	double correction = 0.0;
	for (int i = 0; i <= panels; ++i)
	{
		const double u = static_cast<double>(i) / panels;
		const double weight = (i == 0 || i == panels) ? 1.0 : (i % 2 ? 4.0 : 2.0);
		const double term = weight * u * u *
			std::sqrt(mass * mass + momentum * momentum * u * u) - correction;
		const double updated = sum + term;
		correction = (updated - sum) - term;
		sum = updated;
	}
	return density * sum / panels;
}

double PeEnergy(double density)
{
	return Energy(density, mp) + Energy(density, me);
}

double PeH(double density)
{
	return Mu(density, mp) + Mu(density, me);
}

double PeHessian(double density)
{
	return Derivative(density, mp) + Derivative(density, me);
}

long double IndependentOnsetExpression()
{
	// Independent closed form using the exact binary values of the Zaki double
	// constants. This does not call the provider or its threshold routine. The
	// separately supplied 80-digit review value is checked below; arm64 long
	// double is not assumed to provide more precision than double.
	const long double neutron_mass = Zaki::Physics::NEUTRON_M_MEV;
	const long double proton_mass = Zaki::Physics::PROTON_M_MEV;
	const long double electron_mass = Zaki::Physics::ELECTRON_M_MEV;
	const long double hbar_c = 1.0L / static_cast<long double>(Zaki::Physics::MEV_2_INV_FM);
	const long double mu_e = ((neutron_mass - proton_mass) *
		(neutron_mass + proton_mass) + electron_mass * electron_mass) /
		(2.0L * neutron_mass);
	const long double momentum = std::sqrt((mu_e - electron_mass) *
		(mu_e + electron_mass));
	const long double high_pi = std::acos(-1.0L);
	return momentum * momentum * momentum /
		(3.0L * high_pi * high_pi * hbar_c * hbar_c * hbar_c);
}

long double LdDensityFromMomentum(long double momentum)
{
	const long double high_pi = std::acos(-1.0L);
	const long double hbar_c = 1.0L / static_cast<long double>(Zaki::Physics::MEV_2_INV_FM);
	return momentum * momentum * momentum /
		(3.0L * high_pi * high_pi * hbar_c * hbar_c * hbar_c);
}

long double LdMomentum(long double density)
{
	const long double high_pi = std::acos(-1.0L);
	const long double hbar_c = 1.0L / static_cast<long double>(Zaki::Physics::MEV_2_INV_FM);
	return hbar_c * std::cbrt(3.0L * high_pi * high_pi * density);
}

long double LdEnergy(long double density, long double mass)
{
	constexpr int panels = 16384;
	const long double momentum = LdMomentum(density);
	long double sum = 0.0L;
	long double correction = 0.0L;
	for (int i = 0; i <= panels; ++i)
	{
		const long double u = static_cast<long double>(i) / panels;
		const long double weight = (i == 0 || i == panels) ? 1.0L : (i % 2 ? 4.0L : 2.0L);
		const long double term = weight * u * u *
			std::sqrt(mass * mass + momentum * momentum * u * u) - correction;
		const long double updated = sum + term;
		correction = (updated - sum) - term;
		sum = updated;
	}
	return density * sum / panels;
}

long double LdDensityFromMu(long double chemical_potential, long double mass)
{
	return LdDensityFromMomentum(std::sqrt((chemical_potential - mass) *
		(chemical_potential + mass)));
}

long double LdDerivative(long double density, long double mass)
{
	const long double high_pi = std::acos(-1.0L);
	const long double hbar_c = 1.0L / static_cast<long double>(Zaki::Physics::MEV_2_INV_FM);
	const long double momentum = hbar_c * std::cbrt(3.0L * high_pi * high_pi * density);
	return high_pi * high_pi * hbar_c * hbar_c * hbar_c /
		(std::hypot(mass, momentum) * momentum);
}

struct LongDoubleNpeState
{
	long double neutron_density;
	long double electron_density;
};

LongDoubleNpeState IndependentNpeLongDouble(long double baryon_density)
{
	const long double neutron_mass = Zaki::Physics::NEUTRON_M_MEV;
	const long double proton_mass = Zaki::Physics::PROTON_M_MEV;
	const long double electron_mass = Zaki::Physics::ELECTRON_M_MEV;
	long double lower = 0.0L;
	long double upper = 2.0L;
	for (int i = 0; i < 256; ++i)
	{
		const long double momentum = (lower + upper) / 2.0L;
		if (std::hypot(proton_mass, momentum) +
				std::hypot(electron_mass, momentum) < neutron_mass)
			lower = momentum;
		else
			upper = momentum;
	}
	lower = (lower + upper) / 2.0L;
	upper = 2.0L;
	for (int i = 0; i < 256; ++i)
	{
		const long double momentum = (lower + upper) / 2.0L;
		const long double electron_density = LdDensityFromMomentum(momentum);
		const long double mu_n = std::hypot(proton_mass, momentum) +
			std::hypot(electron_mass, momentum);
		const long double neutron_density = LdDensityFromMu(mu_n, neutron_mass);
		if (neutron_density + electron_density < baryon_density)
			lower = momentum;
		else
			upper = momentum;
	}
	const long double momentum = (lower + upper) / 2.0L;
	const long double electron_density = LdDensityFromMomentum(momentum);
	const long double mu_n = std::hypot(proton_mass, momentum) +
		std::hypot(electron_mass, momentum);
	return {LdDensityFromMu(mu_n, neutron_mass), electron_density};
}

template <class T, class = void>
struct HasHessian : std::false_type
{
};

template <class T>
struct HasHessian<T, std::void_t<decltype(std::declval<T>().hessian)>> : std::true_type
{
};

static_assert(std::tuple_size_v<PeChemicalHessian::Matrix> == 1);
static_assert(!std::is_convertible_v<PeChemicalHessian, NpeChemicalHessian>);
static_assert(!HasHessian<NeutronThresholdEvaluation>::value);
static_assert(!HasHessian<VacuumBoundaryEvaluation>::value);

void V1(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double onset = provider.NeutronOnsetBaryonDensityFm3();
	int count = 0;
	for (double factor : {1e-12, 1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.9})
	{
		const double density = factor * onset;
		const auto result = std::get<PeThermodynamicEvaluation>(provider.EquilibriumAt(density));
		const auto &state = result.state;
		Require(state.BaryonDensityFm3() == density &&
			state.ProtonDensityFm3() == density &&
			state.ElectronDensityFm3() == density &&
			state.NeutronDensityFm3() == 0.0 && state.MuonDensityFm3() == 0.0,
			"PE-V1 falsifier: p-e composition is not exact");
		Require(!result.active_particles.neutron && result.active_particles.proton &&
			result.active_particles.electron && !result.active_particles.muon &&
			result.response_dimension == 1 &&
			result.domain_status == LocalResponseDomainStatus::SmoothInterior,
			"PE-V1 falsifier: active-set or response dimension is wrong");
		++count;
	}
	std::cout << "PE-V1 PASS: " << count
		<< " logarithmic densities have exact n_p=n_e=n_B, n_n=n_mu=0 and a 1D active set\n";
}

void V2(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double onset = provider.NeutronOnsetBaryonDensityFm3();
	double max_energy_relative = 0.0;
	double max_h_error = 0.0;
	for (double factor : {1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.99})
	{
		const double density = factor * onset;
		const auto result = provider.EvaluatePe({density});
		const double energy = PeEnergy(density);
		max_energy_relative = std::max(max_energy_relative,
			std::abs(result.energy_density_MeV_fm3 / energy - 1.0));
		max_h_error = std::max(max_h_error,
			std::abs(result.conjugates.HPeMeV() - PeH(density)));
	}
	Require(max_energy_relative < 3e-12 && max_h_error < 3e-13,
		"PE-V2 falsifier: energy or conjugate differs from independent ideal-gas expressions");
	std::cout << "PE-V2 PASS: independent p/e phase-space energy relative="
		<< max_energy_relative << "; max |h_pe-(mu_p+mu_e)|=" << max_h_error << " MeV\n";
}

void V3(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double onset = provider.NeutronOnsetBaryonDensityFm3();
	double max_relative = 0.0;
	for (double factor : {1e-9, 1e-6, 1e-3, 0.1, 0.5, 0.99})
	{
		const double density = factor * onset;
		const double expected = PeHessian(density);
		const double actual = provider.EvaluatePe({density}).hessian(0, 0);
		max_relative = std::max(max_relative, std::abs(actual / expected - 1.0));
	}
	Require(max_relative < 3e-14,
		"PE-V3 falsifier: 1x1 Hessian is not independent D_p+D_e");
	std::cout << "PE-V3 PASS: 1x1 H_pe=D_p+D_e; max independent relative error="
		<< max_relative << "\n";
}

void V4(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double center = 0.3 * provider.NeutronOnsetBaryonDensityFm3();
	double previous = std::numeric_limits<double>::infinity();
	double minimum_order = std::numeric_limits<double>::infinity();
	double final_error = 0.0;
	for (double scale : {1e-2, 5e-3, 2.5e-3, 1.25e-3})
	{
		const double delta = scale * center;
		const double plus = provider.EvaluatePe({center + delta}).conjugates.HPeMeV();
		const double minus = provider.EvaluatePe({center - delta}).conjugates.HPeMeV();
		const double finite_derivative = (plus - minus) / (2.0 * delta);
		const double analytic = provider.EvaluatePe({center}).hessian(0, 0);
		const double error = std::abs(finite_derivative / analytic - 1.0);
		std::cout << "PE-V4 diagnostic scale=" << scale << " relative error=" << error << "\n";
		if (std::isfinite(previous))
			minimum_order = std::min(minimum_order, std::log(previous / error) / std::log(2.0));
		previous = error;
		final_error = error;
	}
	Require(minimum_order > 1.8 && final_error < 3e-7,
		"PE-V4 falsifier: Delta h does not converge to H_pe Delta n_B");
	std::cout << "PE-V4 PASS: centered finite perturbations converge with min order="
		<< minimum_order << "; final relative derivative error=" << final_error << "\n";
}

void V5(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double onset = provider.NeutronOnsetBaryonDensityFm3();
	double previous = std::numeric_limits<double>::infinity();
	double final_diagnostic = 0.0;
	for (double factor : {1e-6, 1e-4, 1e-2, 0.1, 0.5, 0.9, 0.99, 0.999999})
	{
		const double density = factor * onset;
		const auto result = provider.EvaluatePe({density});
		const double expected = mn - PeH(density);
		Near(result.eta_npe_threshold_diagnostic_MeV, expected, 3e-13,
			"PE-V5 falsifier: inactive-neutron diagnostic value");
		Require(result.eta_npe_threshold_diagnostic_MeV > 0.0 &&
			result.eta_npe_threshold_diagnostic_MeV < previous,
			"PE-V5 falsifier: neutron diagnostic is not positive and monotone");
		previous = result.eta_npe_threshold_diagnostic_MeV;
		final_diagnostic = previous;
	}
	Require(final_diagnostic < 1e-6,
		"PE-V5 falsifier: neutron diagnostic does not approach zero");
	std::cout << "PE-V5 PASS: m_n-mu_p-mu_e is positive, strictly decreasing, final="
		<< final_diagnostic << " MeV\n";
}

void V6(const TrackRFreeGasThermodynamicProvider &provider)
{
	const long double independent = IndependentOnsetExpression();
	const long double reviewed = 7.3567289037328326656352e-9L;
	Require(std::abs((independent - reviewed) / reviewed) < 2e-16L,
		"PE-V6 falsifier: independent onset expression misses reviewed high-precision reference");
	const double production = provider.NeutronOnsetBaryonDensityFm3();
	Require(std::abs(static_cast<long double>(production) / independent - 1.0L) < 3e-16L,
		"PE-V6 falsifier: production neutron onset differs from independent expression");
	std::cout << std::setprecision(22)
		<< "PE-V6 PASS: independent neutron-onset expression=" << independent
		<< " fm^-3; production=" << production << " fm^-3\n";
	std::cout << std::setprecision(17);
}

void V7(const TrackRFreeGasThermodynamicProvider &provider)
{
	const long double onset_high = IndependentOnsetExpression();
	const long double proton_mass = Zaki::Physics::PROTON_M_MEV;
	const long double electron_mass = Zaki::Physics::ELECTRON_M_MEV;
	const long double proton_mu = std::hypot(proton_mass, LdMomentum(onset_high));
	const long double electron_mu = std::hypot(electron_mass, LdMomentum(onset_high));
	const long double pressure_high = onset_high * (proton_mu + electron_mu) -
		LdEnergy(onset_high, proton_mass) - LdEnergy(onset_high, electron_mass);
	const double onset = static_cast<double>(onset_high);
	const double pressure = static_cast<double>(pressure_high);
	const double reviewed = 1.8964875026317866e-9;
	Near(pressure, reviewed, 1e-21,
		"PE-V7 falsifier: independently computed onset pressure");
	Require(pressure > 0.0, "PE-V7 falsifier: onset pressure is not positive");
	double previous_pressure = pressure;
	double final_ratio = 0.0;
	for (double factor : {0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6})
	{
		const double density = factor * onset;
		const double value = density * PeH(density) - PeEnergy(density);
		Require(value > 0.0 && value < previous_pressure,
			"PE-V7 falsifier: p-e pressure is not positive or does not approach zero");
		previous_pressure = value;
		final_ratio = value / density;
	}
	Require(final_ratio < 1e-4,
		"PE-V7 falsifier: P/n_B does not approach its zero-density limit");
	std::cout << "PE-V7 PASS: P_onset=" << pressure
		<< " MeV fm^-3 >0; P>0 below onset and P/n_B->0 (final " << final_ratio << " MeV)\n";
}

void V8(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double onset = provider.NeutronOnsetBaryonDensityFm3();
	const auto value = provider.EquilibriumAt(onset);
	const auto &threshold = std::get<NeutronThresholdEvaluation>(value);
	Require(threshold.state.NeutronDensityFm3() == 0.0 &&
		threshold.state.MuonDensityFm3() == 0.0 &&
		threshold.state.ProtonDensityFm3() == onset &&
		threshold.state.ElectronDensityFm3() == onset &&
		threshold.response_dimension == 0 &&
		threshold.domain_status == LocalResponseDomainStatus::SpeciesThreshold,
		"PE-V8 falsifier: exact neutron threshold is not value-only");
	Near(threshold.limiting_pe_conjugates.HPeMeV(), mn, 3e-13,
		"PE-V8 falsifier: limiting pe conjugate at neutron onset");
	Near(threshold.energy_density_MeV_fm3, PeEnergy(onset), 3e-17,
		"PE-V8 falsifier: neutron-threshold energy value");
	Near(threshold.eta_npe_threshold_diagnostic_MeV, 0.0, 3e-13,
		"PE-V8 falsifier: neutron threshold diagnostic");
	Require(std::holds_alternative<NeutronThresholdEvaluation>(
			provider.EvaluateActive(threshold.state.Coordinates())),
		"PE-V8 falsifier: exact threshold state does not roundtrip through active dispatch");
	Require(std::holds_alternative<PeThermodynamicEvaluation>(
			provider.EquilibriumAt(std::nextafter(onset, 0.0))) &&
		Throws<EquilibriumResolutionError>([&]
		{
			(void)provider.EquilibriumAt(std::nextafter(onset,
				std::numeric_limits<double>::infinity()));
		}), "PE-V8 falsifier: nextafter sides were tolerance-reclassified");
	Require(Throws([&] { (void)provider.EvaluatePe({onset}); }) &&
		Throws([&] { (void)provider.EvaluateNpe({onset, onset}); }),
		"PE-V8 falsifier: ordinary smooth Hessian exists at neutron threshold");
	std::cout << "PE-V8 PASS: exact value-only neutron threshold; lower nextafter is pe, "
		"upper nextafter remains Npe-classified and fails resolution explicitly\n";
}

void V9(const TrackRFreeGasThermodynamicProvider &provider)
{
	const long double onset = IndependentOnsetExpression();
	const long double proton_mass = Zaki::Physics::PROTON_M_MEV;
	const long double electron_mass = Zaki::Physics::ELECTRON_M_MEV;
	const long double A0 = LdDerivative(onset, proton_mass) +
		LdDerivative(onset, electron_mass);
	double max_tangent_relative = 0.0;
	double previous_inverse_error = std::numeric_limits<double>::infinity();
	double final_production_inverse_error = 0.0;
	for (double relative : {1e-3, 1e-4, 1e-5, 1e-6})
	{
		const double baryon_density = provider.NeutronOnsetBaryonDensityFm3() * (1.0 + relative);
		const auto result = std::get<NpeThermodynamicEvaluation>(provider.EquilibriumAt(baryon_density));
		const double h00 = result.hessian(0, 0);
		const double h01 = result.hessian(0, 1);
		const double h11 = result.hessian(1, 1);
		const double tangent = h00 + 2.0 * h01 + h11;
		const double expected_tangent = PeHessian(result.state.ElectronDensityFm3());
		max_tangent_relative = std::max(max_tangent_relative,
			std::abs(tangent / expected_tangent - 1.0));
		const double determinant = h00 * h11 - h01 * h01;
		const std::array<std::array<double, 2>, 2> inverse{{
			{{h11 / determinant, -h01 / determinant}},
			{{-h01 / determinant, h00 / determinant}}}};
		const double A = tangent;
		Near(inverse[0][0], 1.0 / A + 1.0 / h00, 5e-14 / A,
			"PE-V9 falsifier: exact inverse 00 identity");
		for (const auto &entry : {inverse[0][1], inverse[1][0], inverse[1][1]})
			Near(entry, 1.0 / A, 5e-14 / A,
				"PE-V9 falsifier: inverse tangent rank-one identity");
		const double inverse_error = std::abs(inverse[0][0] * static_cast<double>(A0) - 1.0);
		Require(inverse_error < previous_inverse_error,
			"PE-V9 falsifier: production inverse does not move toward pe tangent limit");
		previous_inverse_error = inverse_error;
		final_production_inverse_error = inverse_error;
	}
	Require(max_tangent_relative < 2e-14,
		"PE-V9 falsifier: t^T H_npe t is not H_pe");

	// Carry the independently constructed analytic matrix closer than the
	// production resolution boundary to falsify the limiting normalization.
	long double previous = std::numeric_limits<long double>::infinity();
	long double final = 0.0L;
	for (long double relative : {1e-4L, 1e-6L, 1e-8L, 1e-10L, 1e-12L})
	{
		const auto state = IndependentNpeLongDouble(onset * (1.0L + relative));
		const long double Dn = LdDerivative(state.neutron_density,
			static_cast<long double>(Zaki::Physics::NEUTRON_M_MEV));
		const long double A = LdDerivative(state.electron_density, proton_mass) +
			LdDerivative(state.electron_density, electron_mass);
		const long double inverse00 = 1.0L / A + 1.0L / Dn;
		const long double error = std::abs(inverse00 * A0 - 1.0L);
		Require(error < previous,
			"PE-V9 falsifier: independent inverse does not collapse onto tt^T/A0");
		previous = error;
		final = error;
	}
	Require(final < 0.02L,
		"PE-V9 falsifier: independently normalized inverse limit is not approached");
	std::cout << "PE-V9 PASS: t^T H_npe t=H_pe (max relative "
		<< max_tangent_relative << "); H^-1 -> tt^T/H_pe with no 1/2 factor; "
		<< "production/independent final normalized errors="
		<< final_production_inverse_error << "/" << static_cast<double>(final) << "\n";
}

void V10(const TrackRFreeGasThermodynamicProvider &provider)
{
	const auto value = provider.EquilibriumAt(0.0);
	const auto &vacuum = std::get<VacuumBoundaryEvaluation>(value);
	Require(vacuum.state.BaryonDensityFm3() == 0.0 &&
		vacuum.state.NeutronDensityFm3() == 0.0 &&
		vacuum.state.ProtonDensityFm3() == 0.0 &&
		vacuum.state.ElectronDensityFm3() == 0.0 &&
		vacuum.state.MuonDensityFm3() == 0.0 &&
		vacuum.energy_density_MeV_fm3 == 0.0 &&
		vacuum.response_dimension == 0 &&
		vacuum.domain_status == LocalResponseDomainStatus::VacuumBoundary,
		"PE-V10 falsifier: vacuum state/value/status is wrong");
	Near(vacuum.limiting_pe_conjugates.HPeMeV(), mp + me, 0.0,
		"PE-V10 falsifier: vacuum h_pe one-sided limit");
	Near(vacuum.eta_npe_threshold_diagnostic_MeV, mn - (mp + me), 0.0,
		"PE-V10 falsifier: vacuum neutron value diagnostic");
	Require(std::holds_alternative<VacuumBoundaryEvaluation>(
			provider.EvaluateActive({0.0, 0.0, 0.0})) &&
		Throws([&] { (void)provider.EvaluateActive({0.0, 1e-30, 0.0}); }),
		"PE-V10 falsifier: vacuum active dispatch accepts a non-vacuum composition");
	double previous_h_error = std::numeric_limits<double>::infinity();
	double previous_hessian = 0.0;
	for (double factor : {1e-3, 1e-5, 1e-7, 1e-9})
	{
		const double density = factor * provider.NeutronOnsetBaryonDensityFm3();
		const auto pe = provider.EvaluatePe({density});
		const double h_error = pe.conjugates.HPeMeV() - (mp + me);
		Require(h_error >= 0.0 && h_error < previous_h_error &&
			pe.hessian(0, 0) > previous_hessian,
			"PE-V10 falsifier: wrong vacuum one-sided h/H limits");
		previous_h_error = h_error;
		previous_hessian = pe.hessian(0, 0);
	}
	std::cout << "PE-V10 PASS: vacuum is value-only, epsilon=0, h_pe->m_p+m_e; "
		"H_pe grows without a fabricated finite boundary value\n";
}

void V11(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double onset = provider.NeutronOnsetBaryonDensityFm3();
	int available = 0;
	int failed_closed = 0;
	double max_h00_relative = 0.0;
	for (double relative : {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10})
	{
		const double density = onset * (1.0 + relative);
		Require(provider.EquilibriumDomainAt(density) == FreeGasEquilibriumDomain::Npe,
			"PE-V11 falsifier: physical classification changed by numerical policy");
		try
		{
			const auto result = std::get<NpeThermodynamicEvaluation>(provider.EquilibriumAt(density));
			const auto oracle = IndependentNpeLongDouble(static_cast<long double>(density));
			const long double expected = LdDerivative(oracle.neutron_density,
				static_cast<long double>(Zaki::Physics::NEUTRON_M_MEV));
			max_h00_relative = std::max(max_h00_relative,
				std::abs(result.hessian(0, 0) / static_cast<double>(expected) - 1.0));
			const double ulp = density - std::nextafter(density, 0.0);
			Require(result.state.NeutronDensityFm3() >=
				TrackRFreeGasThermodynamicProvider::MinimumResolvedNpeNeutronUlps() * ulp,
				"PE-V11 falsifier: response returned below documented ULP rule");
			++available;
		}
		catch (const EquilibriumResolutionError &)
		{
			const auto state = provider.EquilibriumStateAt(density);
			Require(state.BaryonDensityFm3() == density &&
				state.NeutronDensityFm3() > 0.0 && state.MuonDensityFm3() == 0.0 &&
				state.ProtonDensityFm3() == state.ElectronDensityFm3(),
				"PE-V11 falsifier: reliable composition values were replaced with another branch");
			++failed_closed;
		}
	}
	Require(available == 4 && failed_closed == 4 && max_h00_relative < 1e-6,
		"PE-V11 falsifier: N-1 policy boundary or achieved H00 accuracy changed");
	std::cout << "PE-V11 PASS: 4 resolved responses through 1e-6 above onset, "
		"4 explicit failures from 1e-7 through 1e-10; max allowed H00 relative error="
		<< max_h00_relative << "; no floor/substitution\n";
}

void V12(const TrackRFreeGasThermodynamicProvider &provider)
{
	const double neutron = provider.NeutronOnsetBaryonDensityFm3();
	const double muon = provider.MuonOnsetBaryonDensityFm3();
	const double sigma = provider.SigmaMinusOnsetBaryonDensityFm3();
	Require(Throws([&] { (void)provider.EquilibriumAt(-1e-12); }) &&
		std::holds_alternative<VacuumBoundaryEvaluation>(provider.EquilibriumAt(0.0)) &&
		std::holds_alternative<PeThermodynamicEvaluation>(provider.EquilibriumAt(0.5 * neutron)) &&
		std::holds_alternative<NeutronThresholdEvaluation>(provider.EquilibriumAt(neutron)) &&
		std::holds_alternative<NpeThermodynamicEvaluation>(provider.EquilibriumAt(0.1)) &&
		std::holds_alternative<MuonThresholdEvaluation>(provider.EquilibriumAt(muon)) &&
		std::holds_alternative<LocalThermodynamicEvaluation>(provider.EquilibriumAt(0.5)) &&
		Throws([&] { (void)provider.EquilibriumAt(sigma); }),
		"PE-V12 falsifier: full active-set dispatch ladder is incomplete or ambiguous");
	Require(provider.EquilibriumDomainAt(0.0) == FreeGasEquilibriumDomain::Vacuum &&
		provider.EquilibriumDomainAt(0.5 * neutron) == FreeGasEquilibriumDomain::ProtonElectron &&
		provider.EquilibriumDomainAt(neutron) == FreeGasEquilibriumDomain::NeutronOnset &&
		provider.EquilibriumDomainAt(0.1) == FreeGasEquilibriumDomain::Npe &&
		provider.EquilibriumDomainAt(muon) == FreeGasEquilibriumDomain::MuonOnset &&
		provider.EquilibriumDomainAt(0.5) == FreeGasEquilibriumDomain::NpeMuon,
		"PE-V12 falsifier: domain enum does not match typed dispatch");
	std::cout << "PE-V12 PASS: invalid / vacuum / pe / neutron threshold / npe / "
		"muon threshold / npe-mu / Sigma rejection dispatched explicitly\n";
}

void V13(const TrackRFreeGasThermodynamicProvider &provider)
{
	Near(provider.MuonOnsetBaryonDensityFm3(), 0.4569848054124199, 2e-15,
		"PE-V13 falsifier: muon onset regression");
	Near(provider.SigmaMinusOnsetBaryonDensityFm3(), 0.61735520796653, 2e-14,
		"PE-V13 falsifier: Sigma-minus onset regression");
	const auto npe = std::get<NpeThermodynamicEvaluation>(provider.EquilibriumAt(0.1));
	Require(npe.response_dimension == 2 && !npe.active_particles.muon,
		"PE-V13 falsifier: established npe response changed");
	const auto full = std::get<LocalThermodynamicEvaluation>(provider.EquilibriumAt(0.5));
	const auto direct = provider.Evaluate(full.state.Coordinates());
	Require(full.hessian.value_MeV_fm3 == direct.hessian.value_MeV_fm3 &&
		full.conjugates.value_MeV == direct.conjugates.value_MeV &&
		full.energy_density_MeV_fm3 == direct.energy_density_MeV_fm3,
		"PE-V13 falsifier: established npe-mu full path changed");
	std::cout << "PE-V13 PASS: muon/Sigma thresholds, npe 2D, and npe-mu 3D paths preserved; "
		"V1-V10, RFG1-RFG11 and R1-V1-R1-V10 remain separate regression executables\n";
}

} // namespace

int main()
{
	try
	{
		std::cout << std::setprecision(17);
		const TrackRFreeGasThermodynamicProvider provider;
		V1(provider);
		V2(provider);
		V3(provider);
		V4(provider);
		V5(provider);
		V6(provider);
		V7(provider);
		V8(provider);
		V9(provider);
		V10(provider);
		V11(provider);
		V12(provider);
		V13(provider);
		std::cout << "Phase 5A-5 PE-V1 through PE-V13 PASS\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << "\n";
		return 1;
	}
}
