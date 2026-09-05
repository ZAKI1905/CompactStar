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
using Matrix2 = std::array<std::array<double, 2>, 2>;
const double mn = Zaki::Physics::NEUTRON_M_MEV;
const double mp = Zaki::Physics::PROTON_M_MEV;
const double me = Zaki::Physics::ELECTRON_M_MEV;
const double mm = Zaki::Physics::MUON_M_MEV;
const double hc = 1.0 / Zaki::Physics::MEV_2_INV_FM;
const double pi = std::acos(-1.0);

void Require(bool ok, const std::string &message)
{
	if (!ok)
		throw std::runtime_error(message);
}
void Near(double a, double b, double tolerance, const std::string &message)
{
	if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tolerance)
		throw std::runtime_error(message + ": actual=" + std::to_string(a) +
								 ", expected=" + std::to_string(b));
}
template <class Error = std::runtime_error, class F> bool Throws(F &&f)
{
	try
	{
		f();
	}
	catch (const Error &)
	{
		return true;
	}
	return false;
}
// Test-only phase-space state counting, g=2. No production fermion is used
// in any equilibrium, derivative, energy, or threshold oracle below.
double DensityFromMomentum(double p)
{
	return 2.0 * (4.0 * pi * p * p * p / 3.0) / std::pow(2.0 * pi * hc, 3);
}
double Momentum(double n) { return hc * std::cbrt(3.0 * pi * pi * n); }
double Mu(double n, double m) { return std::sqrt(m * m + std::pow(Momentum(n), 2)); }
double DensityFromMu(double mu, double m)
{
	Require(mu >= m, "Independent inverse outside particle domain");
	return DensityFromMomentum(std::sqrt((mu - m) * (mu + m)));
}
double Derivative(double n, double m)
{
	// Reciprocal dn/dmu from the phase-space derivative.
	const double p = Momentum(n), mu = std::sqrt(m * m + p * p);
	return pi * pi * hc * hc * hc / (mu * p);
}
double Energy(double n, double m)
{
	// Direct Simpson quadrature in u=p/p_F, including rest mass.
	if (n == 0.0)
		return 0.0;
	constexpr int panels = 4096;
	const double p = Momentum(n);
	double sum = 0.0, correction = 0.0;
	for (int i = 0; i <= panels; ++i)
	{
		const double u = double(i) / panels;
		const double weight = (i == 0 || i == panels) ? 1.0 : (i % 2 ? 4.0 : 2.0);
		const double term = weight * u * u * std::sqrt(m * m + p * p * u * u) - correction;
		const double updated = sum + term;
		correction = (updated - sum) - term;
		sum = updated;
	}
	return n * sum / panels;
}
struct OracleState
{
	double nn, np, ne, nmu;
	double EnergyDensity() const
	{
		return Energy(nn, mn) + Energy(np, mp) + Energy(ne, me) + Energy(nmu, mm);
	}
};
// The lower boundary is independently solved in common proton/electron momentum.
double NeutronMomentum()
{
	double lo = 0.0, hi = 2.0;
	Require(std::hypot(mp, lo) + std::hypot(me, lo) < mn &&
				std::hypot(mp, hi) + std::hypot(me, hi) > mn,
			"R1-V9 neutron bracket");
	for (int i = 0; i < 200; ++i)
	{
		const double p = (lo + hi) / 2.0;
		if (std::hypot(mp, p) + std::hypot(me, p) < mn)
			lo = p;
		else
			hi = p;
	}
	return (lo + hi) / 2.0;
}
OracleState NpeAtMomentum(double p)
{
	const double ne = DensityFromMomentum(p);
	const double mun = std::hypot(mp, p) + std::hypot(me, p);
	return {DensityFromMu(mun, mn), ne, ne, 0.0};
}
OracleState MuonOnset() { return NpeAtMomentum(std::sqrt((mm - me) * (mm + me))); }
OracleState IndependentNpe(double b)
{
	// Different independent variable and equation from production's F(n_e)
	// at fixed b: invert the increasing equilibrium b(p)=n_n(mu_p+mu_e)+n_e(p).
	double lo = NeutronMomentum(), hi = std::sqrt((mm - me) * (mm + me));
	for (int i = 0; i < 200; ++i)
	{
		const double p = (lo + hi) / 2.0;
		const auto s = NpeAtMomentum(p);
		if (s.nn + s.ne < b)
			lo = p;
		else
			hi = p;
	}
	return NpeAtMomentum((lo + hi) / 2.0);
}
OracleState IndependentFull(double b)
{
	// Invert total equilibrium baryon density as a function of common mu.
	// This does not subtract n_p from a supplied b inside the root equation.
	auto at = [](double l)
	{
		const double ne = DensityFromMu(l, me), nmu = DensityFromMu(l, mm), np = ne + nmu;
		return OracleState{DensityFromMu(Mu(np, mp) + l, mn), np, ne, nmu};
	};
	double lo = mm, hi = 150.0;
	for (int i = 0; i < 200; ++i)
	{
		const double l = (lo + hi) / 2.0;
		const auto s = at(l);
		if (s.nn + s.np < b)
			lo = l;
		else
			hi = l;
	}
	return at((lo + hi) / 2.0);
}
Matrix2 IndependentH(double b, double ne)
{
	// S_z is the species-density Jacobian; contraction independently tests
	// matrix structure instead of copying the production 2x2 initializer.
	const std::array<std::array<double, 2>, 3> S{{{1, -1}, {0, 1}, {0, 1}}};
	const std::array<double, 3> D{{Derivative(b - ne, mn), Derivative(ne, mp), Derivative(ne, me)}};
	Matrix2 h{};
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			for (int k = 0; k < 3; ++k)
				h[i][j] += S[k][i] * D[k] * S[k][j];
	return h;
}

void V1(const TrackRFreeGasThermodynamicProvider &p)
{
	double max_eta = 0.0, max_density_error = 0.0, max_energy_relative = 0.0;
	for (double b : {1.01 * p.NeutronOnsetBaryonDensityFm3(), 2 * p.NeutronOnsetBaryonDensityFm3(),
					 1e-7, 1e-5, 1e-3, 0.01, 0.1, 0.3, 0.99 * p.MuonOnsetBaryonDensityFm3()})
	{
		const auto a = std::get<NpeThermodynamicEvaluation>(p.EquilibriumAt(b));
		const auto ref = IndependentNpe(b);
		const auto &s = a.state;
		Require(s.MuonDensityFm3() == 0.0 && s.NeutronDensityFm3() > 0.0 &&
					s.ProtonDensityFm3() == s.ElectronDensityFm3() && s.ElectronDensityFm3() > 0.0,
				"R1-V1 inactive-muon or active-density/neutrality failure");
		const double error = std::abs(s.ElectronDensityFm3() - ref.ne);
		Near(s.ElectronDensityFm3(), ref.ne, 2e-11 * ref.ne, "R1-V1 independent composition");
		max_density_error = std::max(max_density_error, error);
		max_eta = std::max(max_eta, std::abs(a.conjugates.EtaNpeMeV()));
		Near(a.conjugates.MuNMeV(), Mu(s.NeutronDensityFm3(), mn), 2e-12,
			 "R1-V1 intrinsic neutron mu");
		Near(a.conjugates.EtaNpeMeV(),
			 Mu(b - s.ElectronDensityFm3(), mn) - Mu(s.ElectronDensityFm3(), mp) -
				 Mu(s.ElectronDensityFm3(), me),
			 3e-12, "R1-V1 independent eta sign/value");
		// Preserve the requested baryon density in the energy comparison. The
		// oracle's inverse mu_n->n_n loses relative precision near neutron
		// appearance; its independently solved n_e is checked separately above.
		const double energy = OracleState{b - ref.ne, ref.ne, ref.ne, 0.0}.EnergyDensity();
		max_energy_relative =
			std::max(max_energy_relative, std::abs(a.energy_density_MeV_fm3 / energy - 1.0));
	}
	Require(max_eta < 5e-11 && max_energy_relative < 3e-12,
			"R1-V1 equilibrium or independent energy mismatch");
	std::cout << "R1-V1 PASS: 9 independent momentum solves; max |eta|=" << max_eta
			  << " MeV; max density error=" << max_density_error
			  << " fm^-3; energy relative=" << max_energy_relative << "\n";
}
void V2(const TrackRFreeGasThermodynamicProvider &p)
{
	const double onset = p.MuonOnsetBaryonDensityFm3();
	const auto ref = MuonOnset();
	Near(onset, ref.nn + ref.ne, 2e-14, "R1-V2 independent rest-mass muon onset");
	for (double factor : {0.99, 1.0, 1.01})
	{
		const auto s = p.EquilibriumStateAt(factor * onset);
		const double mue = Mu(s.ElectronDensityFm3(), me);
		if (factor < 1)
			Require(mue < mm, "R1-V2 muons absent above chemical threshold");
		else if (factor > 1)
			Require(mue > mm && s.MuonDensityFm3() > 0, "R1-V2 muons missing above threshold");
		else
			Near(mue, mm, 4e-14, "R1-V2 exact muon threshold value");
	}
	std::cout << "R1-V2 PASS: independent onset=" << ref.nn + ref.ne
			  << " fm^-3; strict inequalities and rest-mass equality\n";
}
void V3(const TrackRFreeGasThermodynamicProvider &p)
{
	double max_relative = 0.0;
	for (const NpeCoordinates z :
		 {NpeCoordinates{1e-6, 1e-8}, NpeCoordinates{0.05, 0.0002}, NpeCoordinates{0.3, 0.002}})
	{
		const auto a = p.EvaluateNpe(z);
		Require(std::abs(a.conjugates.EtaNpeMeV()) > 1e-4, "R1-V3 fixture must be off equilibrium");
		const auto ref = IndependentH(z.n_B_fm3, z.n_e_fm3);
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				max_relative = std::max(max_relative, std::abs(a.hessian(i, j) / ref[i][j] - 1.0));
	}
	Require(max_relative < 3e-14, "R1-V3 independent S_z^T K S_z Hessian mismatch");
	std::cout << "R1-V3 PASS: independent analytic contraction max relative error=" << max_relative
			  << "\n";
}
void V4(const TrackRFreeGasThermodynamicProvider &p)
{
	const NpeCoordinates z{0.18, 0.0012};
	const auto a = p.EvaluateNpe(z);
	double min_order = 10.0, max_final = 0.0, max_gradient = 0.0;
	for (const std::array<double, 2> direction :
		 {std::array<double, 2>{0.18, 0}, {0, 0.0012}, {0.18, -0.0012}})
	{
		std::array<double, 4> errors{};
		for (int k = 0; k < 4; ++k)
		{
			const double t = 1e-3 / std::pow(2.0, k), db = t * direction[0], de = t * direction[1];
			const auto b = p.EvaluateNpe({z.n_B_fm3 + db, z.n_e_fm3 + de});
			for (int row = 0; row < 2; ++row)
				errors[k] = std::max(
					errors[k], std::abs(b.conjugates.value_MeV[row] - a.conjugates.value_MeV[row] -
										a.hessian(row, 0) * db - a.hessian(row, 1) * de));
			if (k)
				min_order = std::min(min_order, std::log2(errors[k - 1] / errors[k]));
		}
		max_final = std::max(max_final, errors.back());
		std::cout << "R1-V4 diagnostic errors=" << errors[0] << "," << errors[1] << "," << errors[2]
				  << "," << errors[3] << " MeV\n";
	}
	// An independent energy integral additionally pins h=d epsilon/dz at this
	// off-equilibrium fixture, so a mutually consistent wrong h/H is rejected.
	auto energy = [](double b, double ne)
	{ return Energy(b - ne, mn) + Energy(ne, mp) + Energy(ne, me); };
	for (int axis = 0; axis < 2; ++axis)
	{
		const double step = (axis == 0 ? z.n_B_fm3 : z.n_e_fm3) * 1e-4;
		const double db = axis == 0 ? step : 0, de = axis == 1 ? step : 0;
		const double fd =
			(energy(z.n_B_fm3 + db, z.n_e_fm3 + de) - energy(z.n_B_fm3 - db, z.n_e_fm3 - de)) /
			(2 * step);
		max_gradient = std::max(max_gradient, std::abs(fd - a.conjugates.value_MeV[axis]));
	}
	Require(min_order > 1.95 && max_final < 2e-6 && max_gradient < 2e-6,
			"R1-V4 no quadratic remainder or wrong energy gradient");
	std::cout << "R1-V4 PASS: min order=" << min_order << "; max final=" << max_final
			  << " MeV; independent energy-gradient error=" << max_gradient << " MeV\n";
}
void V5(const TrackRFreeGasThermodynamicProvider &p)
{
	int count = 0;
	double min_d = 1e300, min_det = 1e300, max_asym = 0;
	const double low = 1.01 * p.NeutronOnsetBaryonDensityFm3(),
				 high = 0.999 * p.MuonOnsetBaryonDensityFm3();
	for (int k = 0; k <= 40; ++k)
	{
		const double b = low * std::pow(high / low, k / 40.0);
		const auto s = IndependentNpe(b);
		for (double factor : {0.8, 1.0, 1.2})
		{
			const double ne = s.ne * factor;
			if (ne >= b || Mu(ne, me) >= mm)
				continue; // grid restricted to declared open chart
			const auto h = p.EvaluateNpe({b, ne}).hessian;
			const double det = h(0, 0) * h(1, 1) - h(0, 1) * h(1, 0);
			Require(h(0, 0) > 0 && det > 0, "R1-V5 nonpositive Sylvester minor");
			Near(det, Derivative(b - ne, mn) * (Derivative(ne, mp) + Derivative(ne, me)),
				 3e-13 * det, "R1-V5 analytic determinant");
			min_d = std::min(min_d, h(0, 0));
			min_det = std::min(min_det, det);
			max_asym = std::max(max_asym, std::abs(h(0, 1) - h(1, 0)));
			++count;
		}
	}
	Require(count > 100 && max_asym == 0.0, "R1-V5 inadequate grid or asymmetric Hessian");
	std::cout << "R1-V5 PASS: " << count << " log-density/composition states; min H00=" << min_d
			  << "; min determinant=" << min_det << "; max asymmetry=" << max_asym << "\n";
}
void V6(const TrackRFreeGasThermodynamicProvider &p)
{
	const double onset = p.MuonOnsetBaryonDensityFm3();
	const auto limit = MuonOnset();
	const auto hlimit = IndependentH(limit.nn + limit.ne, limit.ne);
	const double elim = limit.EnergyDensity(), munlim = Mu(limit.nn, mn);
	double prev_error = 1e300, prev_D = 0, final_error = 0, final_block = 0, final_D = 0;
	for (double delta : {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9})
	{
		const auto a = std::get<NpeThermodynamicEvaluation>(p.EquilibriumAt(onset * (1 - delta)));
		const auto b = std::get<LocalThermodynamicEvaluation>(p.EquilibriumAt(onset * (1 + delta)));
		const auto ref_b = IndependentFull(onset * (1 + delta));
		Near(b.state.ElectronDensityFm3(), ref_b.ne, 3e-14,
			 "R1-V6 independent above-onset composition");
		Near(b.state.MuonDensityFm3(), ref_b.nmu, 3e-14, "R1-V6 independent above-onset muons");
		Require(a.state.MuonDensityFm3() == 0.0 && b.state.MuonDensityFm3() > 0,
				"R1-V6 active-set continuity");
		double error = 0, block = 0;
		for (const auto &pair :
			 {std::pair<double, double>{a.state.NeutronDensityFm3(), b.state.NeutronDensityFm3()},
			  {a.state.ProtonDensityFm3(), b.state.ProtonDensityFm3()},
			  {a.state.ElectronDensityFm3(), b.state.ElectronDensityFm3()}})
			error = std::max(error, std::abs(pair.first - pair.second) / onset);
		for (const auto &triple :
			 {std::array<double, 3>{a.state.NeutronDensityFm3(), b.state.NeutronDensityFm3(),
									limit.nn},
			  {a.state.ProtonDensityFm3(), b.state.ProtonDensityFm3(), limit.np},
			  {a.state.ElectronDensityFm3(), b.state.ElectronDensityFm3(), limit.ne}})
			for (int side = 0; side < 2; ++side)
				error = std::max(error, std::abs(triple[side] - triple[2]) / onset);
		const double mue_below = Mu(a.state.ElectronDensityFm3(), me),
					 mue_above = Mu(b.state.ElectronDensityFm3(), me);
		Require(mue_below < mm && mue_above > mm,
				"R1-V6 threshold inequality lost during approach");
		Near(mue_above, Mu(b.state.MuonDensityFm3(), mm), 3e-12, "R1-V6 common lepton mu");
		Require(std::abs(mue_below - mm) < 70 * delta && std::abs(mue_above - mm) < 70 * delta,
				"R1-V6 lepton potentials do not approach rest mass");
		for (double en : {a.energy_density_MeV_fm3, b.energy_density_MeV_fm3})
			error = std::max(error, std::abs(en / elim - 1));
		for (double mun : {a.conjugates.MuNMeV(), b.conjugates.MuNMeV()})
			error = std::max(error, std::abs(mun / munlim - 1));
		for (double eta :
			 {a.conjugates.EtaNpeMeV(), b.conjugates.EtaNpeMeV(), b.conjugates.EtaNpmuMeV()})
			Near(eta, 0, 5e-11, "R1-V6 beta-equilibrium continuity");
		Require(std::abs(a.eta_npmu_threshold_diagnostic_MeV) < 70 * delta,
				"R1-V6 inactive diagnostic not approaching zero");
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				block = std::max({block, std::abs(a.hessian(i, j) / hlimit[i][j] - 1),
								  std::abs(b.hessian(i, j) / hlimit[i][j] - 1)});
		const double D = b.hessian(2, 2) - b.hessian(2, 1);
		Near(D, Derivative(b.state.MuonDensityFm3(), mm), 3e-14 * D,
			 "R1-V6 inactive derivative clipped from above");
		Require(error < 3 * delta && block < 3 * delta && error < prev_error && D > prev_D,
				"R1-V6 discontinuity, shared block mismatch, or bounded D_mu");
		prev_error = error;
		prev_D = D;
		final_error = error;
		final_block = block;
		final_D = D;
	}
	std::cout << "R1-V6 PASS: two-sided delta=1e-2..1e-9; final scaled continuity error="
			  << final_error << "; shared block relative=" << final_block
			  << "; D_mu increased to=" << final_D << " MeV fm^3\n";
}
// Compile-time type checks falsify accidental unification of 2x2 and 3x3 responses.
template <class T, class = void> struct HasHessian : std::false_type
{
};
template <class T>
struct HasHessian<T, std::void_t<decltype(std::declval<T>().hessian)>> : std::true_type
{
};
static_assert(!std::is_convertible_v<NpeChemicalHessian, ChargeNeutralChemicalHessian>);
static_assert(!std::is_convertible_v<NpeThermodynamicEvaluation, LocalThermodynamicEvaluation>);
static_assert(!HasHessian<MuonThresholdEvaluation>::value);
static_assert(std::tuple_size_v<NpeChemicalHessian::Matrix> == 2);
static_assert(std::tuple_size_v<ChargeNeutralChemicalHessian::Matrix> == 3);
void V7(const TrackRFreeGasThermodynamicProvider &p)
{
	const ILocalThermodynamicProvider &generic = p;
	for (double b : {0.1, 0.99 * p.MuonOnsetBaryonDensityFm3(), p.MuonOnsetBaryonDensityFm3()})
	{
		const auto state = p.EquilibriumStateAt(b);
		Require(Throws([&] { (void)generic.Evaluate(state.Coordinates()); }),
				"R1-V7 ordinary full Hessian at n_mu=0");
		const auto a = generic.EquilibriumAt(b);
		Require(std::get_if<LocalThermodynamicEvaluation>(&a) == nullptr,
				"R1-V7 inactive full response in variant");
		std::visit(
			[](const auto &v)
			{
				Require(!v.active_particles.muon && v.response_dimension < 3,
						"R1-V7 per-result particle/dimension mismatch");
			},
			a);
	}
	Require(!ColdRelativisticIdealFermion::Muon().Evaluate(0).dchemical_potential_dn_MeV_fm3,
			"R1-V7 finite zero-density derivative");
	for (double nmu : {1e-18, 1e-24, 1e-30})
	{
		const auto a = generic.EvaluateActive({0.3, 0.002, nmu});
		const auto &full = std::get<LocalThermodynamicEvaluation>(a);
		Require(full.state.MuonDensityFm3() == nmu, "R1-V7 positive density floor");
		Near(full.hessian(2, 2) - full.hessian(2, 1), Derivative(nmu, mm),
			 3e-14 * Derivative(nmu, mm), "R1-V7 tiny active derivative");
	}
	std::cout << "R1-V7 PASS: compile-time separation, generic dispatch, full rejection at zero, "
				 "positive densities unfloored\n";
}
void V8(const TrackRFreeGasThermodynamicProvider &p)
{
	const double onset = p.MuonOnsetBaryonDensityFm3();
	const double below = std::nextafter(onset, 0.0), above = std::nextafter(onset, 1.0);
	Require(p.EquilibriumDomainAt(below) == FreeGasEquilibriumDomain::Npe &&
				p.EquilibriumDomainAt(onset) == FreeGasEquilibriumDomain::MuonOnset &&
				p.EquilibriumDomainAt(above) == FreeGasEquilibriumDomain::NpeMuon,
			"R1-V8 one-ULP domain ambiguity");
	const auto v = p.EquilibriumAt(onset);
	const auto &a = std::get<MuonThresholdEvaluation>(v);
	Require(a.domain_status == LocalResponseDomainStatus::SpeciesThreshold &&
				a.response_dimension == 0 && a.state.MuonDensityFm3() == 0.0,
			"R1-V8 threshold must have state but no smooth response");
	Near(a.limiting_npe_conjugates.EtaNpeMeV(), 0, 3e-12, "R1-V8 threshold eta_npe");
	Near(a.eta_npmu_threshold_diagnostic_MeV, 0, 3e-12, "R1-V8 threshold diagnostic");
	Require(
		std::holds_alternative<MuonThresholdEvaluation>(p.EvaluateActive(a.state.Coordinates())),
		"R1-V8 exact state roundtrip");
	Require(Throws([&] { (void)p.EvaluateNpe({onset, a.state.ElectronDensityFm3()}); }),
			"R1-V8 threshold ordinary npe Hessian");
	// The one-ULP side is physical and unambiguously classified, but a double
	// common chemical potential/root cannot certify that separation. No state
	// or Hessian is returned and no nearby density is substituted.
	for (double ne : {std::nextafter(a.state.ElectronDensityFm3(), 0.0),
					  std::nextafter(a.state.ElectronDensityFm3(), 1.0)})
		Require(Throws([&] { (void)p.EvaluateActive({onset, ne, 0.0}); }),
				"R1-V8 nearby composition silently replaced by onset state");
	for (double b : {below, above})
	{
		Require(Throws<EquilibriumResolutionError>([&] { (void)p.EquilibriumAt(b); }),
				"R1-V8 unresolvable near-threshold response fabricated");
		Require(Throws<EquilibriumResolutionError>([&] { (void)p.EquilibriumStateAt(b); }),
				"R1-V8 composition API bypassed precision failure");
	}
	for (double b : {onset * (1 - 1e-10), onset * (1 + 1e-10)})
		Require(!std::holds_alternative<MuonThresholdEvaluation>(p.EquilibriumAt(b)),
				"R1-V8 finite-width threshold switching");
	std::cout << "R1-V8 PASS: exact onset value-only; nextafter sides classified distinctly and "
				 "resolution errors explicit; 1e-10 relative sides smooth\n";
}
void V9(const TrackRFreeGasThermodynamicProvider &p)
{
	const double b0 = DensityFromMomentum(NeutronMomentum());
	Near(p.NeutronOnsetBaryonDensityFm3(), b0, 2e-12 * b0, "R1-V9 independent neutron onset");
	const double at = p.NeutronOnsetBaryonDensityFm3();
	Near(Mu(at, mp) + Mu(at, me), mn, 3e-13, "R1-V9 neutron zero-density condition");
	for (double b : {at * 0.01, at * 0.1, at * 0.99})
	{
		Require(Mu(b, mp) + Mu(b, me) < mn && Mu(b, me) < mm,
				"R1-V9 absent neutron/muon inequality");
		const double pressure = b * (Mu(b, mp) + Mu(b, me)) - Energy(b, mp) - Energy(b, me);
		Require(pressure > 0, "R1-V9 incorrectly placed zero-pressure surface");
		Require(p.EquilibriumDomainAt(b) == FreeGasEquilibriumDomain::ProtonElectron &&
					std::holds_alternative<PeThermodynamicEvaluation>(p.EquilibriumAt(b)),
				"R1-V9 p-e branch did not use its explicit 1D response");
	}
	Require(p.EquilibriumDomainAt(at) == FreeGasEquilibriumDomain::NeutronOnset &&
				std::holds_alternative<NeutronThresholdEvaluation>(p.EquilibriumAt(at)),
			"R1-V9 neutron onset must remain a distinct value-only gate");
	Require(p.EquilibriumDomainAt(std::nextafter(at, 0.0)) ==
					FreeGasEquilibriumDomain::ProtonElectron &&
				std::holds_alternative<PeThermodynamicEvaluation>(
					p.EquilibriumAt(std::nextafter(at, 0.0))) &&
				p.EquilibriumDomainAt(std::nextafter(at, 1.0)) == FreeGasEquilibriumDomain::Npe,
			"R1-V9 lower nextafter classification");
	Require(
		Throws<EquilibriumResolutionError>([&] { (void)p.EquilibriumAt(std::nextafter(at, 1.0)); }),
		"R1-V9 unresolved neutron appearance fabricated");
	const double pressure = at * (Mu(at, mp) + Mu(at, me)) - Energy(at, mp) - Energy(at, me);
	Require(pressure > 0, "R1-V9 neutron-onset pressure is not positive");
	std::cout << "R1-V9 PASS: neutron onset=" << at << " fm^-3; mu_e=" << Mu(at, me)
			  << " MeV; pressure=" << pressure
			  << " MeV fm^-3 >0; explicit p-e and neutron-threshold dispatch preserved\n";
}
void V10(const TrackRFreeGasThermodynamicProvider &p)
{
	const ILocalThermodynamicProvider &g = p;
	for (double b : {0.48, 0.52, 0.6})
	{
		const auto s = p.EquilibriumStateAt(b);
		const auto old = p.Evaluate(s.Coordinates());
		const auto active = std::get<LocalThermodynamicEvaluation>(g.EquilibriumAt(b));
		Require(old.hessian.value_MeV_fm3 == active.hessian.value_MeV_fm3 &&
					old.conjugates.value_MeV == active.conjugates.value_MeV &&
					old.energy_density_MeV_fm3 == active.energy_density_MeV_fm3 &&
					active.response_dimension == 3 && active.active_particles.muon &&
					active.domain_status == LocalResponseDomainStatus::SmoothInterior,
				"R1-V10 changed full-response semantics");
	}
	for (const ChargeNeutralCoordinates x : {ChargeNeutralCoordinates{0.1, -0.01, 0},
											 {0.1, 0.2, 0},
											 {0.1, 0.001, -1e-30},
											 {0.1, 0.001, std::numeric_limits<double>::quiet_NaN()},
											 {0.1, 0, 0}})
		Require(Throws([&] { (void)g.EvaluateActive(x); }),
				"R1-V10 invalid active coordinates accepted");
	Require(std::holds_alternative<VacuumBoundaryEvaluation>(g.EquilibriumAt(0.0)),
			"R1-V10 vacuum boundary dispatch regression");
	for (double b : {-1.0, std::numeric_limits<double>::infinity(),
					 std::numeric_limits<double>::quiet_NaN(), p.SigmaMinusOnsetBaryonDensityFm3()})
		Require(Throws([&] { (void)g.EquilibriumAt(b); }),
				"R1-V10 invalid equilibrium domain accepted");
	std::cout << "R1-V10 PASS: full-path equivalence and domain regressions; separate V1-V10 and "
				 "RFG1-RFG11 executables required\n";
}
} // namespace
int main()
{
	try
	{
		std::cout << std::setprecision(17);
		const TrackRFreeGasThermodynamicProvider p;
		V1(p);
		V2(p);
		V3(p);
		V4(p);
		V5(p);
		V6(p);
		V7(p);
		V8(p);
		V9(p);
		V10(p);
		std::cout << "Phase 5A-4 R1-V1 through R1-V10 PASS\n";
		return 0;
	}
	catch (const std::exception &e)
	{
		std::cerr << "FAIL: " << e.what() << "\n";
		return 1;
	}
}
