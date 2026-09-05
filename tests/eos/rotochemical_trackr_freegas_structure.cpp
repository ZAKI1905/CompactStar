// Candidate Structure-1 validation. No inverse-mass solve or production ODE changes.
#include "structure1/table.hpp"
#include "structure1/tov_oracle.hpp"
#include <CompactStar/Core/NStar.hpp>
#include <CompactStar/Core/StarProfile.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <chrono>
#include <gsl/gsl_errno.h>
#include <iostream>
using namespace structure1;
using namespace CompactStar::Core;
struct Probe : TOVSolver
{
	using TOVSolver::PressureCutoff;
	double upper() const { return .999 * eos_tab.eps.back(); }
	double floor() const { return eos_tab.pre.front(); }
	const auto &labels() const { return eos_tab.extra_labels; }
	std::vector<std::vector<std::vector<double>>> captured;
	Probe()
	{
		n_exp_cond_f = [](const NStar &)
		{ return true; };
	}
	static auto columns(const NStar &star)
	{
		std::vector<std::vector<double>> a;
		const auto &d = star.Profile().Radial();
		for (size_t k = 0; k < d.Dim().size(); ++k)
		{
			a.emplace_back();
			for (size_t j = 0; j < d[int(k)].Size(); ++j)
				a.back().push_back(d[int(k)][j]);
		}
		return a;
	}
	void ExportNStarProfile(const size_t &, const Zaki::String::Directory &) override { captured.push_back(columns(n_star)); }
};
double lapse(double M, double R, double mass_length) { return std::sqrt(1 - 2 * mass_length * M / R); }
struct Result
{
	double requested, achieved, analytic_rho, pc, cut, M, R, z, ri, Rlo, Rhi, R0, ri0, oracleM, oracleR;
};
Result star(Probe &s, EnthalpyOracle &oracle, double requested, size_t res, std::vector<TOVPoint> *save = nullptr)
{
	require(requested < s.upper(), "central state would be clamped");
	s.SetRadialRes(res);
	std::vector<TOVPoint> p;
	require(s.SingleStarSolveToTOVPoints(requested, p) > 0 && s.LastSolveStatus() == TOVSolveStatus::SURFACE_REACHED, "canonical TOV did not reach surface");
	require(std::abs(p.front().e - requested) < 5e10, "central inversion exceeds density budget");
	require(std::abs(p.back().p / s.PressureCutoff() - 1) < 1e-8, "surface pressure residual");
	const auto &q = p.back();
	double pc = s.GetInitPress(), cut = s.PressureCutoff();
	double arho = rho(oracle.at_H(oracle.H_for_pressure(pc)).e);
	require(std::abs(arho - requested) < 5e10, "analytic central-state error exceeds total budget");
	auto ref = oracle.solve(arho, cut, 1e-12, 1e-10);
	// Direct-EOS enthalpy lower/upper bounds; independent of TOV spline/ODE.
	auto sq = oracle.at_H(ref.Hs);
	double m = q.m * oracle.solar_length;
	double Rhi = 2 * m / (1 - (1 - 2 * m / q.r) * std::exp(2 * ref.Hs));
	double mu = m + 4 * M_PI / 3 * sq.e * oracle.to_geo * (Rhi * Rhi * Rhi - q.r * q.r * q.r);
	require(q.r > 2 * mu, "tail bound metric domain");
	double gu = (mu + 4 * M_PI * Rhi * Rhi * Rhi * sq.p * oracle.to_geo) / (q.r * (q.r - 2 * mu));
	double Rlo = q.r + ref.Hs / gu, R0 = (Rlo + Rhi) / 2;
	double z = lapse(q.m, q.r, Zaki::Physics::SUN_M_KM);
	Result out{requested, p.front().e, arho, pc, cut, q.m, q.r, z, q.r / z, Rlo, Rhi, R0, R0 / lapse(q.m, R0, Zaki::Physics::SUN_M_KM), ref.M, ref.R};
	if (save)
		*save = std::move(p);
	return out;
}
void emit(std::ostream &o, const std::string &group, size_t resolution, size_t partition, double floor, const Result &q)
{
	o << group << '\t' << resolution << '\t' << partition << '\t' << q.requested << '\t' << q.achieved << '\t' << q.analytic_rho << '\t' << q.pc << '\t' << floor << '\t' << q.cut << '\t' << q.cut / q.pc << '\t' << q.M << '\t' << q.R << '\t' << q.z << '\t' << q.ri << '\t' << q.Rlo << '\t' << q.Rhi << '\t' << q.R0 << '\t' << q.ri0 << '\t' << q.oracleM << '\t' << q.oracleR << "\tSURFACE_REACHED" << std::endl;
}
bool bins(const Result &q) { return q.M - 1e-7 >= .615 && q.M + 1e-7 < .625 && q.R0 - 1e-5 >= 12.765 && q.R0 + 1e-5 < 12.775 && q.ri0 - 1e-5 >= 13.795 && q.ri0 + 1e-5 < 13.805; }
int main(int argc, char **argv)
{
	try
	{
		gsl_set_error_handler_off();
		TrackRFreeGasThermodynamicProvider eos;
		EnthalpyOracle oracle;
		std::filesystem::path dir = std::filesystem::absolute(argc > 1 ? std::string(argv[1]) : "structure1-artifacts-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::create_directories(dir);
		std::ofstream o(dir / "stars.tsv");
		o << std::setprecision(17);
		o << "group\tgrid\tpartition\trequested_rho\tachieved_rho\tanalytic_rho\tpc\tfloor\tcut\tcut_over_pc\tM\tR_cut\tlapse_profile\tRinf_cut_profile\tR0_lower\tR0_upper\tR0_est\tRinf0_profile\toracle_M_cut\toracle_R_cut\tstatus\n";
		std::ofstream summary(dir / "summary.txt");
		summary << std::setprecision(17);
		std::ofstream oo(dir / "oracle.tsv");
		oo << std::setprecision(17) << "tol\tstart\tM\tR_0\n";
		OracleStar prev{};
		for (double tol : {1e-9, 1e-10, 1e-11, 1e-12})
		{
			auto q = oracle.solve(1.10e15, 0, tol, tol * 100);
			oo << tol << '\t' << tol * 100 << '\t' << q.M << '\t' << q.R << std::endl;
			if (tol == 1e-12)
				require(std::abs(q.R - prev.R) < 1e-7 && std::abs(q.M - prev.M) < 1e-8, "S12 independent oracle fails refinement");
			prev = q;
		}
		double nb = central_n(eos, 1.10e15);
		summary << "midpoint_nB=" << nb << "\nmidpoint_eps=" << eos.BarotropeAt(nb).energy_density_MeV_fm3 << "\nG_cgs=" << oracle.G << "\nsolar_mass_g=" << GSL_CONST_CGSM_SOLAR_MASS << "\nraw_solar_length_km=" << oracle.solar_length << "\nprofile_solar_length_km=" << Zaki::Physics::SUN_M_KM << std::endl;
		Result last{};
		std::filesystem::path finest;
		for (size_t grid : {1024, 2048, 4096, 8192})
		{
			auto tab = generate(eos, dir / ("eos-" + std::to_string(grid) + ".tsv"), grid);
			finest = tab.file;
			Probe s;
			s.ImportEOS(tab.file.string(), true);
			require(s.upper() > 1.105e15, "S11 publication interval clamp");
			Result base{};
			for (size_t res : {2500, 5000, 10000, 20000, 40000})
			{
				auto q = star(s, oracle, 1.10e15, res);
				emit(o, "radial", grid, res, s.floor(), q);
				if (res == 2500)
					base = q;
				else
					require(std::abs(q.M / base.M - 1) < 1e-9 && std::abs(q.R / base.R - 1) < 1e-8 && std::abs(q.z / base.z - 1) < 1e-8, "S6 partition convergence");
				if (grid == 8192)
				{
					require(std::abs(q.M - q.oracleM) < 1e-7 && std::abs(q.R - q.oracleR) < 1e-6, "S12 material canonical/oracle disagreement");
				}
				if (grid == 8192 && res == 40000)
				{
					require(std::abs(q.R - last.R) < 1e-6 && std::abs(q.M - last.M) < 1e-7, "EOS table convergence");
					last = q;
				}
			}
			if (grid < 8192)
				last = base;
			summary << "table=" << grid << " rows=" << tab.rows.size() << " upper_n=" << tab.rows.back().n_B_fm3 << " Sigma_margin_n=" << eos.SigmaMinusOnsetBaryonDensityFm3() - tab.rows.back().n_B_fm3 << " upper_rho=" << s.upper() / .999 << " clamp_rho=" << s.upper() << " publication_clamp_margin=" << s.upper() - 1.105e15 << '\n';
		}
		summary << "midpoint_M=" << last.M << "\nmidpoint_Rcut=" << last.R << "\nmidpoint_R0=" << last.R0 << "\nmidpoint_Rinfcut=" << last.ri << "\nmidpoint_Rinf0=" << last.ri0 << "\nmidpoint_all_bins=" << bins(last) << std::endl;
		// S5: real sequence finalization and bulk NStar, same canonical points.
		Probe s;
		s.SetWrkDir(dir.string());
		s.ImportEOS(finest.string(), true);
		std::vector<TOVPoint> points;
		auto mid = star(s, oracle, 1.10e15, 10000, &points);
		NStar ns(points, s.labels());
		auto cols = Probe::columns(ns);
		require(ns.Profile().RadiusSurface() == mid.R && std::abs(ns.Profile().ExpNuSurface() - mid.z) < 1e-14, "S9 profile metric identity");
		s.Solve(Zaki::Math::Axis{{1.10e15, 1.10e15}, 1, "Linear"}, "s5/", "star");
		require(s.captured.size() == 2, "S5 supported sequence did not finalize");
		for (const auto &a : s.captured)
			require(a == cols, "S5 raw radial/species/metric columns differ");
		summary << "S5_sequence_all_columns_bit_identical=1\nS9_profile_R_over_lapse_identity=1\n";
		// S7 floor sweep. The independent oracle matches the achieved pressure state.
		double final_R0 = mid.R0;
		for (int power : {11, 12, 13, 14, 15, 16, 18})
		{
			double target = mid.pc * std::pow(10., -power), a = 1e-20, b = eos.NeutronOnsetBaryonDensityFm3();
			for (int j = 0; j < 90; ++j)
			{
				double n = (a + b) / 2;
				if (n == a || n == b)
					break;
				double p = eos.BarotropeAt(n).pressure_MeV_fm3 * CompactStar::Units::MEV_FM3_TO_ERG_CM3;
				if (p < target)
					a = n;
				else
					b = n;
			}
			auto tab = generate(eos, dir / ("floor-" + std::to_string(power) + ".tsv"), 8192, (a + b) / 2);
			Probe cut;
			cut.ImportEOS(tab.file.string(), true);
			auto q = star(cut, oracle, 1.10e15, 10000);
			emit(o, "floor", power, 10000, cut.floor(), q);
			if (power >= 15)
			{
				require(q.Rhi - q.Rlo < 1e-5, "S7 tail uncertainty cannot reach 5 cm");
				require(std::abs(q.R0 - final_R0) < 1e-5, "S7 floor plateau/convergence");
			}
		}
		// S8: nested entire printed interval; upper endpoint is explicitly excluded.
		double lo = 1.095e15, hi = std::nextafter(1.105e15, lo), previous_min = 0, previous_max = 0;
		std::vector<Result> old_sequence;
		for (int count : {17, 33, 65})
		{
			double minR = 1e99, maxR = 0, matchlo = 1e99, matchhi = 0, pm = 0, pr = 1e99, curve_error = 0;
			std::vector<Result> sequence;
			for (int j = 0; j < count; ++j)
			{
				double r = lo + (hi - lo) * double(j) / (count - 1);
				auto q = star(s, oracle, r, 10000);
				emit(o, "rho-bin", count, 10000, s.floor(), q);
				require(std::abs(q.R - q.oracleR) < 1e-6 && std::abs(q.M - q.oracleM) < 1e-7, "S12 rho-bin oracle disagreement");
				sequence.push_back(q);
				if (!old_sequence.empty() && j % 2)
					curve_error = std::max(curve_error, std::abs(q.R0 - .5 * (old_sequence[j / 2].R0 + old_sequence[j / 2 + 1].R0)));
				require(q.M > pm && q.R0 < pr, "S8 mass/radius branch monotonicity");
				pm = q.M;
				pr = q.R0;
				minR = std::min(minR, q.R0);
				maxR = std::max(maxR, q.R0);
				if (bins(q))
				{
					matchlo = std::min(matchlo, r);
					matchhi = std::max(matchhi, r);
				}
			}
			if (count > 17)
				require(std::abs(minR - previous_min) < 1e-7 && std::abs(maxR - previous_max) < 1e-7, "S8 interval image convergence");
			previous_min = minR;
			previous_max = maxR;
			if (count == 65)
				require(curve_error < 1e-5, "S8 central sampling budget");
			old_sequence = sequence;
			summary << "curve_interpolation_error=" << curve_error << "\n";
			summary << "rho_grid=" << count << " R0_min=" << minR << " R0_max=" << maxR << " matching_sample_min=" << matchlo << " matching_sample_max=" << matchhi << '\n';
		}
		double pm = 0;
		for (int j = 0; j <= 16; ++j)
		{
			double n = .59 + (.615 - .59) * j / 16.;
			auto value = eos.BarotropeAt(n);
			auto q = star(s, oracle, rho(value.energy_density_MeV_fm3), 10000);
			emit(o, "high-branch", 17, 10000, s.floor(), q);
			require(q.M > pm, "S8 high-branch mass monotonicity");
			pm = q.M;
		}
		summary << "largest_sample_mass=" << pm << " at_nB=.615; source truncated, not unrestricted maximum\n";
		// Domain adapters reject before a TOV or sequence call. No clamped star is accepted.
		for (double n : {eos.SigmaMinusOnsetBaryonDensityFm3(), .7})
		{
			bool refused = false;
			try
			{
				(void)eos.BarotropeAt(n);
			}
			catch (const std::exception &)
			{
				refused = true;
			}
			require(refused, "S11 central/sequence source guard");
			refused = false;
			try
			{
				(void)generate(eos, dir / "forbidden.tsv", 1024, 1e-14, n);
			}
			catch (const std::exception &)
			{
				refused = true;
			}
			require(refused, "S11 generator source guard");
		}
		summary << "S6=PASS\nS7=PASS\nS8=PASS\nS9=PASS\nS11=PASS\nS12=PASS\nS10_midpoint=" << (bins(mid) ? "PASS" : "FAIL") << std::endl;
		return 0;
	}
	catch (const std::exception &e)
	{
		std::cerr << "STRUCTURE-1 FAILURE: " << e.what() << '\n';
		return 1;
	}
}
