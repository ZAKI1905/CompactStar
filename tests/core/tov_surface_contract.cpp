// ADR-0009 validation. Diagnostic output only; never writes a governed baseline.
#include "../rotation/hartle_thorne_1968_hw_eos.hpp"
#include <CompactStar/Core/NStar.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/Geometry.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
using namespace CompactStar::Core;
namespace {
void Require(bool condition, const std::string &why) {
	if (!condition)
		throw std::runtime_error(why);
}
struct Probe : TOVSolver {
	using TOVSolver::PressureCutoff;
	size_t exports = 0;
	Probe() {
		n_exp_cond_f = [](const NStar &) { return true; };
	}
	void ExportNStarProfile(const size_t &, const Zaki::String::Directory &) override { ++exports; }
	void Limits(double lo, double hi) {
		r_min = lo;
		r_max = hi;
	}
	const auto &Labels() const { return eos_tab.extra_labels; }
	double Floor() const { return eos_tab.pre.front(); }
	std::pair<double, double> Range() const {
		return {central_eps_floor_factor * eos_tab.eps.front(), 0.999 * eos_tab.eps.back()};
	}
	double Prepare(double ec) {
		const auto range = Range();
		init_press = p_of_e(std::clamp(ec, range.first, range.second));
		return init_press;
	}
};
constexpr std::array<double, 4> densities = {454550405078491.75, 616488270506054.5,
											 731253342677476.12, 1298349261929558.8};
constexpr std::array<double, 4> targets = {1., 1.4, 1.6, 2.};
constexpr std::array<double, 4> oldM = {1.0000438094972419, 1.4000217830781583, 1.5999758341716293,
										2.0000952962048859};
constexpr std::array<double, 4> oldR = {13.426323081955102, 13.545323064955092, 13.468323075955098,
										12.712323183955165};
constexpr std::array<double, 4> expectedR = {1.2, 1., 5., 3.7};
constexpr std::array<double, 5> former = {1321824566394199.8, 1351861473269815.2, 1413998457266518.,
										  1429973936092014., 1529689333192825.5};
constexpr std::array<size_t, 5> resolutions = {2500, 5000, 10000, 20000, 40000};
void Validate(Probe &s, const std::vector<TOVPoint> &p) {
	Require(s.LastSolveStatus() == TOVSolveStatus::SURFACE_REACHED && !p.empty(),
			"not SURFACE_REACHED");
	const double cut = s.PressureCutoff();
	Require(std::fabs(p.back().p / cut - 1) <= 1e-8, "surface pressure");
	size_t exact = 0;
	for (size_t i = 0; i < p.size(); ++i) {
		const auto &q = p[i];
		exact += q.p == cut;
		Require(q.p >= cut && (!i || q.r > p[i - 1].r), "sub-cutoff / non-increasing node");
		Require(std::isfinite(q.m) && std::isfinite(q.dedp) && q.nu == 0., "invalid point field");
	}
	Require(exact == 1, "surface must be stored exactly once");
	const auto &q = p.back();
	Require(q.e == s.GetEDens(cut) && q.rho == s.GetRho(cut) && q.rho_i == s.GetRho_i(cut) &&
				q.dedp == s.GetEDensDeriv(cut),
			"surface EOS field mismatch");
	const double nu = s.GetNuDer(q.r * 1e5, {q.m * GSL_CONST_CGSM_SOLAR_MASS, cut});
	Require(std::fabs(q.nu_der / nu - 1) < 1e-13, "surface nu derivative");
}
std::vector<TOVPoint> Solve(Probe &s, double ec, size_t res) {
	s.SetRadialRes(res);
	std::vector<TOVPoint> p;
	Require(s.SingleStarSolveToTOVPoints(ec, p) > 0,
			"primitive failure ec=" + std::to_string(ec) + " res=" + std::to_string(res) +
				" status=" + std::to_string(int(s.LastSolveStatus())));
	Validate(s, p);
	return p;
}
bool Equal(const std::vector<TOVPoint> &a, const std::vector<TOVPoint> &b) {
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); ++i) {
		const auto &x = a[i], &y = b[i];
		if (x.r != y.r || x.m != y.m || x.nu_der != y.nu_der || x.nu != y.nu || x.p != y.p ||
			x.e != y.e || x.rho != y.rho || x.rho_i != y.rho_i || x.dedp != y.dedp)
			return false;
	}
	return true;
}
void V7a(Probe &s, std::ostream &o) {
	o << "gate\tlabel\tec\told_M\tnew_M\tdM_rel\told_R\tnew_R\tdR_m\n";
	for (size_t i = 0; i < 4; ++i) {
		auto p = Solve(s, densities[i], 10000);
		const auto &q = p.back();
		const double dm = q.m / oldM[i] - 1, dr = (q.r - oldR[i]) * 1000;
		o << "V7a\t" << targets[i] << '\t' << densities[i] << '\t' << oldM[i] << '\t' << q.m << '\t'
		  << dm << '\t' << oldR[i] << '\t' << q.r << '\t' << dr << std::endl;
		Require(std::fabs(dm) <= 1e-9, "V7a mass impact");
		Require(std::fabs(dr - expectedR[i]) <= .5, "V7a radius impact");
	}
}
void V7b(Probe &s, std::ostream &o) {
	o << "gate\ttarget\tM\tresidual\tec\tstatus\tec_lo\tec_hi\tM_lo\tM_hi\texcluded_"
		 "coarse\tfallback\tdeterministic\n";
	for (double target : targets) {
		s.SetRadialRes(10000);
		std::vector<TOVPoint> a, b;
		std::vector<std::string> la, lb;
		Require(s.SolveToProfile(target, a, &la) > 0, "target solve failed");
		Validate(s, a);
		Require(std::fabs(a.back().m - target) < 1e-4, "V7b residual");
		Require(s.SolveToProfile(target, b, &lb) > 0, "repeat target solve failed");
		Validate(s, b);
		Require(Equal(a, b) && la == lb, "target nondeterminism");
		// Test-side audit of the unchanged N=24 scan. Only complete samples enter a bracket.
		auto range = s.Range();
		std::vector<double> ec, m;
		size_t failed = 0;
		for (int i = 0; i <= 24; ++i) {
			double e = std::pow(10., std::log10(range.first) +
										 (double(i) / 24) *
											 (std::log10(range.second) - std::log10(range.first)));
			std::vector<TOVPoint> p;
			int n = s.SingleStarSolveToTOVPoints(e, p);
			if (n <= 0) {
				Require(p.empty() && s.LastSolveStatus() != TOVSolveStatus::SURFACE_REACHED,
						"failed coarse publication");
				++failed;
				continue;
			}
			Validate(s, p);
			ec.push_back(e);
			m.push_back(p.back().m);
		}
		size_t k = 0;
		while (k + 1 < m.size() && !(m[k + 1] > m[k] && m[k] <= target && target <= m[k + 1]))
			++k;
		Require(k + 1 < m.size(), "no stable complete bracket");
		Require(a.front().e >= ec[k] && a.front().e <= ec[k + 1],
				"target outside selected bracket");
		o << "V7b\t" << target << '\t' << a.back().m << '\t' << std::fabs(a.back().m - target)
		  << '\t' << a.front().e << "\tSURFACE_REACHED\t" << ec[k] << '\t' << ec[k + 1] << '\t'
		  << m[k] << '\t' << m[k + 1] << '\t' << failed << "\tno\tyes\n";
	}
}
struct Surface {
	double R, M, p, pc;
};
using Driver = std::unique_ptr<gsl_odeiv2_driver, decltype(&gsl_odeiv2_driver_free)>;
Driver MakeDriver(gsl_odeiv2_system &sys) {
	return Driver(gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, .1, 1e-10, 1e-10),
				  gsl_odeiv2_driver_free);
}
// Independent test-only locator: radial re-integration and bisection, never the
// production pressure-coordinate helper. EOS and radial physics remain production-owned.
Surface Radial(Probe &s, double ec, size_t res, double uniform_cm = 0., bool reset = false) {
	double pc = s.Prepare(ec), cut = s.PressureCutoff(), r = 1.;
	double y[2] = {(4. / 3.) * M_PI * r * r * r * s.GetEDens(pc), pc};
	gsl_odeiv2_system sys = {TOVSolver::ODE, nullptr, 2, &s};
	auto driver = MakeDriver(sys);
	const double step = (70e5 - 1.) / res;
	for (double target = 1.; target <= 70e5;) {
		double r0 = r, m0 = y[0], p0 = y[1];
		if (reset)
			driver = MakeDriver(sys); // discard adaptive step size as well as RK state
		Require(gsl_odeiv2_driver_apply(driver.get(), &r, target, y) == GSL_SUCCESS,
				"radial reference GSL failure");
		if (y[1] <= cut) {
			double lo = r0, hi = r;
			double bestM = m0, bestP = p0;
			for (int j = 0; j < 64 && hi - lo > 1e-6; ++j) {
				double mid = lo + (hi - lo) / 2, rr = r0, yy[2] = {m0, p0};
				auto d = MakeDriver(sys);
				Require(gsl_odeiv2_driver_apply(d.get(), &rr, mid, yy) == GSL_SUCCESS,
						"radial root re-integration failure");
				if (yy[1] >= cut) {
					lo = mid;
					bestM = yy[0];
					bestP = yy[1];
				} else
					hi = mid;
			}
			Require(hi - lo <= 1e-6, "radial root bracket did not contract");
			return {lo / 1e5, bestM / GSL_CONST_CGSM_SOLAR_MASS, bestP, pc};
		}
		double scale = target < 100		 ? .005
					   : target < 1000	 ? .025
					   : target < 10000	 ? .05
					   : target < 100000 ? .25
										 : 1.;
		target += uniform_cm > 0 ? uniform_cm : step * scale;
	}
	throw std::runtime_error("radial reference exhausted r_max");
}
double Spread(const std::vector<double> &v, double reference) {
	return (*std::max_element(v.begin(), v.end()) - *std::min_element(v.begin(), v.end())) /
		   std::fabs(reference);
}
void Sweep(Probe &s, std::ostream &o) {
	std::vector<double> ecs;
	for (int i = 0; i < 150; ++i)
		ecs.push_back(3e14 * std::pow(1.6e15 / 3e14, double(i) / 149));
	ecs.insert(ecs.end(), former.begin(), former.end());
	size_t production = 0, radial = 0;
	double maxM = 0, maxR = 0, maxC = 0, maxZ = 0, maxCrossM = 0, maxCrossR = 0, maxP = 0;
	o << "kind\tec\tpartition\tR\tM\tpc\tp_final\tcut\tstatus\n";
	for (double ec : ecs) {
		std::vector<double> rs, ms, cs, zs;
		auto add = [&](const Surface &v) {
			rs.push_back(v.R);
			ms.push_back(v.M);
			double c = v.M * Zaki::Physics::SUN_M_KM / v.R;
			cs.push_back(c);
			zs.push_back(std::sqrt(1 - 2 * c));
		};
		for (size_t res : resolutions) {
			auto p = Solve(s, ec, res);
			auto q = p.back();
			double cut = s.PressureCutoff();
			++production;
			Surface v = {q.r, q.m, q.p, s.GetInitPress()};
			add(v);
			o << "production\t" << ec << '\t' << res << '\t' << q.r << '\t' << q.m << '\t' << v.pc
			  << '\t' << q.p << '\t' << cut << "\tSURFACE_REACHED\n";
			maxP = std::max(maxP, std::fabs(q.p / cut - 1));
			auto c = Radial(s, ec, res);
			++radial;
			maxCrossR = std::max(maxCrossR, std::fabs(c.R / q.r - 1));
			maxCrossM = std::max(maxCrossM, std::fabs(c.M / q.m - 1));
		}
		for (double cm : {175., 700., 2800.}) {
			auto v = Radial(s, ec, 10000, cm);
			++radial;
			add(v);
			o << "uniform_radial\t" << ec << '\t' << cm << '\t' << v.R << '\t' << v.M << '\t'
			  << v.pc << '\t' << v.p << '\t' << s.PressureCutoff() << "\tBRACKET_REFINED\n";
		}
		double dm = Spread(ms, ms[2]), dr = Spread(rs, rs[2]), dc = Spread(cs, cs[2]),
			   dz = Spread(zs, zs[2]);
		maxM = std::max(maxM, dm);
		maxR = std::max(maxR, dr);
		maxC = std::max(maxC, dc);
		maxZ = std::max(maxZ, dz);
		Require(dm <= 1e-9 && dr <= 1e-8 && dc <= 1e-8 && dz <= 1e-8,
				"partition invariance ec=" + std::to_string(ec));
		Require(maxCrossM <= 1e-9 && maxCrossR <= 1e-8, "independent locator disagreement");
	}
	o << "SUMMARY production=" << production << " independent=" << radial
	  << " GSL_failures=0 p_rel=" << maxP << " M_spread=" << maxM << " R_spread=" << maxR
	  << " compactness_spread=" << maxC << " lapse_spread=" << maxZ << " cross_M=" << maxCrossM
	  << " cross_R=" << maxCrossR << std::endl;
}
void AuditAndReset(Probe &s, std::ostream &o) {
	const std::vector<size_t> audit = {1000, 1250, 1500, 1750, 2000, 2250, 2350,  2400,	 2450, 2475,
									   2490, 2495, 2499, 2500, 2501, 2505, 2510,  2525,	 2550, 2600,
									   2750, 3000, 3500, 4000, 5000, 7500, 10000, 20000, 40000};
	o << "gate\tec\tres\tR\tM\n";
	for (size_t res : audit) {
		auto p = Solve(s, 7.312533e14, res);
		o << "V5\t7.312533e14\t" << res << '\t' << p.back().r << '\t' << p.back().m << '\n';
	}
	double dm = 0, dr = 0;
	for (double ec : {densities[0], 7.312533e14, former[0], densities[3]})
		for (size_t res : {2500, 3000, 10000}) {
			auto a = Radial(s, ec, res), b = Radial(s, ec, res, 0., true);
			dm = std::max(dm, std::fabs(a.M / b.M - 1));
			dr = std::max(dr, std::fabs(a.R / b.R - 1));
		}
	o << "RESET M=" << dm << " R=" << dr << std::endl;
	Require(dm <= 1e-9 && dr <= 1e-8, "driver reset sensitivity");
}
void Derivatives(Probe &s, std::ostream &o) {
	std::vector<double> ecs(densities.begin(), densities.end());
	ecs.push_back(7.312533e14);
	ecs.insert(ecs.end(), former.begin(), former.end());
	double maxPartition = 0, maxStep = 0;
	o << "ec\tstep\tpartition\tdM_dpc_km3\n";
	for (double ec : ecs) {
		std::array<double, 3> reference{};
		size_t j = 0;
		for (double delta : {5e-4, 1e-3, 2e-3}) {
			std::vector<double> ds;
			for (size_t res : {2000, 2500, 2550, 4000, 5000, 10000, 20000, 40000}) {
				auto a = Solve(s, ec * (1 - delta), res);
				double pa = s.GetInitPress();
				auto b = Solve(s, ec * (1 + delta), res);
				double pb = s.GetInitPress();
				double d = (b.back().m - a.back().m) * Zaki::Physics::SUN_M_KM /
						   ((pb - pa) *
							(Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2));
				ds.push_back(d);
				o << ec << '\t' << delta << '\t' << res << '\t' << d << '\n';
			}
			for (double cm : {175., 700., 2800.}) {
				auto a = Radial(s, ec * (1 - delta), 10000, cm),
					 b = Radial(s, ec * (1 + delta), 10000, cm);
				double d = (b.M - a.M) * Zaki::Physics::SUN_M_KM /
						   ((b.pc - a.pc) *
							(Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2));
				ds.push_back(d);
				o << ec << '\t' << delta << '\t' << cm << '\t' << d << '\n';
			}
			reference[j++] = ds[5];
			double spread = Spread(ds, ds[5]);
			maxPartition = std::max(maxPartition, spread);
			Require(spread <= 1e-5, "V6 partition derivative instability");
		}
		double d = std::fabs(reference[0] / reference[1] - 1);
		maxStep = std::max(maxStep, d);
		Require(d <= 1e-5, "V6 step-pair instability ec=" + std::to_string(ec));
	}
	o << "SUMMARY partition=" << maxPartition << " step_pair=" << maxStep << std::endl;
}

void FailureContracts(Probe &s, const std::string &eos, std::ostream &o, bool forced) {
	auto failure = [&](const std::string &name, int n, const std::vector<TOVPoint> &p,
					   TOVSolveStatus status) {
		const bool rejected = n == 0 && p.empty() && status != TOVSolveStatus::SURFACE_REACHED;
		o << name << " n=" << n << " size=" << p.size() << " status=" << int(status)
		  << " rejected=" << rejected << std::endl;
		return rejected;
	};
	std::vector<TOVPoint> p;
	if (forced) {
		bool all = true;
		s.SetRadialRes(10000);
		int n = s.SingleStarSolveToTOVPoints(densities[2], p);
		all &= failure("primitive", n, p, s.LastSolveStatus());
		o << "primitive_GSL=" << (s.LastSolveStatus() == TOVSolveStatus::GSL_FAILURE) << std::endl;
		all &= s.LastSolveStatus() == TOVSolveStatus::GSL_FAILURE;
		std::vector<std::string> labels = {"stale"};
		n = s.SolveToProfile(1.6, p, &labels);
		all &= failure("target", n, p, s.LastSolveStatus()) && labels.empty();
		s.exports = 0;
		Zaki::Math::Axis ax{{densities[2], densities[2]}, 1, "Linear"};
		s.Solve(ax, "/", "forced-sequence");
		const bool seq = s.exports == 0 && s.LastSolveStatus() != TOVSolveStatus::SURFACE_REACHED;
		all &= seq;
		o << "Solve rejected=" << seq << " exports=" << s.exports << std::endl;
		s.exports = 0;
		s.GenTestSequence(densities[2], "/", "forced-resolution");
		const bool gen = s.exports == 0 && s.LastSolveStatus() != TOVSolveStatus::SURFACE_REACHED;
		all &= gen;
		o << "GenTestSequence rejected=" << gen << " exports=" << s.exports << std::endl;
		NStar ns;
		ns.SetWrkDir("/tmp/compactstar-tov-surf-ir/");
		int nn = ns.SolveTOV_Profile(eos, 1.6);
		const bool wrapper = nn == 0 && !ns.IsSurfaceFinalized();
		all &= wrapper;
		o << "NStar rejected=" << wrapper << " n=" << nn << std::endl;
		bool derivativeRejected = false;
		try {
			(void)Solve(s, densities[2] * (1 + 5e-4), 10000);
		} catch (const std::runtime_error &) {
			derivativeRejected = true;
		}
		all &= derivativeRejected;
		o << "derivative_builder rejected=" << derivativeRejected << std::endl;
		Require(all, "forced-failure publication detector");
		return;
	}
	p = Solve(s, densities[2], 10000);
	s.SetRadialRes(0);
	int n = s.SingleStarSolveToTOVPoints(densities[2], p);
	Require(failure("zero partition", n, p, s.LastSolveStatus()) &&
				s.LastSolveStatus() == TOVSolveStatus::PARTITION_INVALID,
			"zero partition");
	s.SetRadialRes(10000);
	s.SetMaxRadius(std::numeric_limits<double>::quiet_NaN());
	n = s.SingleStarSolveToTOVPoints(densities[2], p);
	Require(failure("nonfinite limit", n, p, s.LastSolveStatus()) &&
				s.LastSolveStatus() == TOVSolveStatus::PARTITION_INVALID,
			"nonfinite limit");
	s.SetMaxRadius(1.);
	n = s.SingleStarSolveToTOVPoints(densities[2], p);
	Require(failure("radius exhaustion", n, p, s.LastSolveStatus()) &&
				s.LastSolveStatus() == TOVSolveStatus::R_MAX_EXHAUSTED,
			"r_max exhaustion");
	s.SetMaxRadius(70.);
	n = s.SingleStarSolveToTOVPoints(std::numeric_limits<double>::quiet_NaN(), p);
	Require(failure("invalid central", n, p, s.LastSolveStatus()) &&
				s.LastSolveStatus() == TOVSolveStatus::INVALID_INITIAL_STATE,
			"nonfinite initial");
	Probe empty;
	n = empty.SingleStarSolveToTOVPoints(densities[2], p);
	Require(failure("empty EOS", n, p, empty.LastSolveStatus()) &&
				empty.LastSolveStatus() == TOVSolveStatus::EOS_DOMAIN_FAILURE,
			"EOS absence");
	s.Limits(1e4, std::nextafter(1e4, std::numeric_limits<double>::infinity()));
	n = s.SingleStarSolveToTOVPoints(densities[2], p);
	Require(failure("nonincreasing target", n, p, s.LastSolveStatus()) &&
				s.LastSolveStatus() == TOVSolveStatus::PARTITION_INVALID,
			"nonincreasing partition");
	s.Limits(1., 70e5);
	p = Solve(s, densities[2], 10000);
	std::vector<std::string> labels = {"stale"};
	n = s.SolveToProfile(100., p, &labels);
	Require(failure("unbracketed mass", n, p, s.LastSolveStatus()) && labels.empty(),
			"no nearest fallback");
	p = Solve(s, densities[2], 10000);
	o << "status reset success=" << (s.LastSolveStatus() == TOVSolveStatus::SURFACE_REACHED)
	  << std::endl;
}
void SurfaceConsumers(Probe &s, std::ostream &o) {
	o << "label\tec\tM\tR\tcompactness_M_over_R\tlapse\tredshift\tarea_km2\tI_km3\tB\n";
	for (size_t i = 0; i < 4; ++i) {
		auto p = Solve(s, densities[i], 10000);
		NStar ns(p, s.Labels());
		const auto &prof = ns.Profile();
		const double r = (*prof.GetRadius())[-1], m = (*prof.GetMass())[-1],
					 nu = (*prof.GetMetricNu())[-1];
		const double f = CompactStar::Geometry::MetricDenominator(r, m), z = prof.ExpNuSurface();
		Require(r == p.back().r && m == p.back().m * Zaki::Physics::SUN_M_KM,
				"surface profile coordinate");
		Require(std::fabs(nu - .5 * std::log(f)) <= 1e-13 &&
					std::fabs(z / std::sqrt(f) - 1) <= 1e-13,
				"surface lapse boundary");
		Require(std::fabs(std::exp((*prof.GetMetricLambda())[-1]) * std::sqrt(f) - 1) <= 1e-13,
				"surface metric lambda");
		o << targets[i] << '\t' << densities[i] << '\t' << p.back().m << '\t' << r << '\t' << m / r
		  << '\t' << z << '\t' << 1 / z - 1 << '\t' << 4 * M_PI * r * r * z * z << '\t'
		  << ns.GetSequence().I << '\t' << ns.GetSequence().b << '\n';
	}
}

void ReferenceRoot(Probe &s, std::ostream &o) {
	double lo = 1297985389788764.5, hi = densities[3];
	auto prod = Solve(s, lo, 10000);
	double bestec = lo, bestm = prod.back().m;
	for (int i = 0; i < 40; ++i) {
		double mid = .5 * (lo + hi);
		auto p = Solve(s, mid, 10000);
		double m = p.back().m;
		if (std::fabs(m - 2) < std::fabs(bestm - 2)) {
			bestm = m;
			bestec = mid;
		}
		if (m < 2)
			lo = mid;
		else
			hi = mid;
		if (std::fabs(m - 2) < 1e-10)
			break;
	}
	o << "reference ec=" << bestec << " M=" << bestm << " residual=" << std::fabs(bestm - 2)
	  << " bracket_lo=" << lo << " bracket_hi=" << hi
	  << " production_ec=1297985389788764.5 production_residual=" << std::fabs(prod.back().m - 2)
	  << std::endl;
}

void KnownFailures(Probe &s, std::ostream &o) {
	bool complete = true;
	auto one = [&](double ec, size_t res) {
		s.SetRadialRes(res);
		std::vector<TOVPoint> p;
		int n = s.SingleStarSolveToTOVPoints(ec, p);
		bool ok = n > 0 && s.LastSolveStatus() == TOVSolveStatus::SURFACE_REACHED;
		o << "known ec=" << ec << " res=" << res << " status=" << int(s.LastSolveStatus())
		  << " points=" << n << std::endl;
		if (ok)
			Validate(s, p);
		complete &= ok;
	};
	for (size_t res : {2500, 2510, 2525, 3000})
		one(7.312533e14, res);
	for (double ec : former)
		one(ec, 10000);
	Require(complete, "known surface failures recurred");
}

} // namespace
int main(int argc, char **argv) {
	try {
		Require(argc == 4, "usage: tov_surface_contract <eos-file> <v7a|v7b> <diagnostic-output>");
		gsl_set_error_handler_off();
		std::string eos = argv[1];
		if (eos == "HW") {
			eos = std::string(argv[3]) + ".eos";
			Require(ht68::WriteDensifiedEos(eos), "HW fixture write failed");
		}
		Probe s;
		s.SetWrkDir("/tmp/compactstar-tov-surf-ir/");
		s.ImportEOS(eos, true);
		std::ofstream out(argv[3]);
		Require(bool(out), "report open failed");
		out << std::setprecision(17);
		std::string mode = argv[2];
		if (mode == "v7a")
			V7a(s, out);
		else if (mode == "v7b")
			V7b(s, out);
		else if (mode == "sweep")
			Sweep(s, out);
		else if (mode == "audit")
			AuditAndReset(s, out);
		else if (mode == "derivatives")
			Derivatives(s, out);
		else if (mode == "contracts")
			FailureContracts(s, eos, out, false);
		else if (mode == "forced-failure")
			FailureContracts(s, eos, out, true);
		else if (mode == "consumers")
			SurfaceConsumers(s, out);
		else if (mode == "reference-root")
			ReferenceRoot(s, out);
		else if (mode == "known-failures")
			KnownFailures(s, out);
		else
			throw std::runtime_error("unknown mode");
		std::cout << "SURFACE PASS " << mode << std::endl;
		return 0;
	} catch (const std::exception &e) {
		std::cerr << "SURFACE FAIL: " << e.what() << std::endl;
		return 1;
	}
}
