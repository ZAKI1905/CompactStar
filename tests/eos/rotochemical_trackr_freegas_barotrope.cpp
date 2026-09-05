#include "structure1/local_oracle.hpp"
#include "structure1/table.hpp"
#include <CompactStar/Core/TOVSolver.hpp>
#include <gsl/gsl_errno.h>
#include <iostream>
using namespace structure1;
#include <chrono>
int main(int argc, char **argv)
{
	try
	{
		gsl_set_error_handler_off();
		TrackRFreeGasThermodynamicProvider eos;
		LocalOracle oracle;
		std::filesystem::path dir = argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path("structure1-local-artifacts-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::create_directories(dir);
		std::ofstream log(dir / "local.tsv");
		log << std::setprecision(17);
		auto vacuum = eos.BarotropeAt(0);
		require(vacuum.values_resolved && vacuum.energy_density_MeV_fm3 == 0 && vacuum.pressure_MeV_fm3 == 0, "S3 vacuum values");
		for (auto f : {CompactStar::ColdRelativisticIdealFermion::Neutron(), CompactStar::ColdRelativisticIdealFermion::Proton(), CompactStar::ColdRelativisticIdealFermion::Electron(), CompactStar::ColdRelativisticIdealFermion::Muon()})
		{
			double k = f.RestMassEnergyMeV() * 1e-6, n = std::pow(k / f.HbarCMeVFm(), 3) / (3 * M_PI * M_PI);
			auto v = f.Values(n);
			require(std::abs(v.pressure_MeV_fm3 / (n * k * k / (5 * f.RestMassEnergyMeV())) - 1) < 1e-10, "S1 low-density pressure series/degeneracy");
		}
		double maxE = 0, maxP = 0, maxIdentity = 0, maxSlope = 0, maxCs = 0;
		for (double n : {1e-24, 1e-18, 1e-14, 1e-11, 1e-9, 7e-9, 8e-9, 1e-7, 1e-5, .01, .2, .456, .457, .6, .616})
		{
			auto p = eos.BarotropeAt(n);
			auto q = oracle.at_n(n);
			maxE = std::max(maxE, std::abs(p.energy_density_MeV_fm3 / q.e - 1));
			maxP = std::max(maxP, std::abs(p.pressure_MeV_fm3 / q.p - 1));
			if (q.p / q.e > 1e-5)
			{
				maxIdentity = std::max(maxIdentity, std::abs((q.n * q.h - q.e) / q.p - 1));
			}
			const double d = 1e-5;
			auto a = oracle.at_n(n * (1 - d)), b = oracle.at_n(n * (1 + d));
			maxSlope = std::max(maxSlope, std::abs((b.e - a.e) / (b.p - a.p) / q.slope - 1));
			maxCs = std::max(maxCs, 1 / q.slope);
			require(p.energy_density_MeV_fm3 > 0 && p.pressure_MeV_fm3 > 0 && q.slope >= 3, "S4 positive causal values");
		}
		log << "quadrature_energy_rel\t" << maxE << "\nquadrature_pressure_rel\t" << maxP << "\npressure_identity_rel\t" << maxIdentity << "\nslope_difference_rel\t" << maxSlope << "\nmax_cs2\t" << maxCs << std::endl;
		require(maxE < 1e-10 && maxP < 1e-8 && maxIdentity < 1e-8 && maxSlope < 1e-5, "S1/S2 independent local discrepancy");
		for (double d : {1e-7, 1e-9, 1e-10})
		{
			double n = eos.NeutronOnsetBaryonDensityFm3() * (1 + d);
			auto p = eos.BarotropeAt(n);
			auto q = oracle.at_n(n);
			bool refused = false;
			try
			{
				(void)eos.EquilibriumAt(n);
			}
			catch (const CompactStar::EquilibriumResolutionError &)
			{
				refused = true;
			}
			require(refused && p.number_densities_fm3[0] > 0 && std::abs(p.pressure_MeV_fm3 / q.p - 1) < 1e-8, "P7/N1 separation");
			log << "P7\t" << d << '\t' << p.energy_density_MeV_fm3 << '\t' << p.pressure_MeV_fm3 << '\t' << refused << '\n';
		}
		for (double n : {eos.SigmaMinusOnsetBaryonDensityFm3(), .7})
		{
			bool failed = false;
			try
			{
				(void)eos.BarotropeAt(n);
			}
			catch (const std::exception &)
			{
				failed = true;
			}
			require(failed, "S11 domain escaped");
		}
		// Certify value interpolation through tiny roundoff bands adjacent to the exact
		// thresholds. The generator samples both thresholds and both relative 1e-8 sides.
		// Oracle has direct values, even where the response-oriented root refuses.
		for (double t : {eos.NeutronOnsetBaryonDensityFm3(), eos.MuonOnsetBaryonDensityFm3()})
		{
			auto l = eos.BarotropeAt(t * (1 - 1e-8)), m = eos.BarotropeAt(t), r = eos.BarotropeAt(t * (1 + 1e-8));
			double errE = 0, errP = 0, errN = 0;
			for (int j = -100; j <= 100; ++j)
			{
				double n = t * (1 + j * 1e-10);
				auto q = oracle.at_n(n);
				auto a = j < 0 ? l : m, b = j < 0 ? m : r;
				double w = (n - a.n_B_fm3) / (b.n_B_fm3 - a.n_B_fm3);
				errE = std::max(errE, std::abs((a.energy_density_MeV_fm3 * (1 - w) + b.energy_density_MeV_fm3 * w) / q.e - 1));
				errP = std::max(errP, std::abs((a.pressure_MeV_fm3 * (1 - w) + b.pressure_MeV_fm3 * w) / q.p - 1));
				for (int k = 0; k < 4; ++k)
					errN = std::max(errN, std::abs((a.number_densities_fm3[k] * (1 - w) + b.number_densities_fm3[k] * w - q.ns[k]) / n));
			}
			log << "threshold_bridge\t" << t << '\t' << errE << '\t' << errP << '\t' << errN << std::endl;
			require(errE < 1e-9 && errP < 1e-7 && errN < 1e-9, "S3 threshold bridge cannot be certified");
		}
		for (double t : {eos.NeutronOnsetBaryonDensityFm3(), eos.MuonOnsetBaryonDensityFm3()})
		{
			auto limit = oracle.at_n(t);
			double previous = 1e99;
			for (double d : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10})
			{
				auto left = oracle.at_n(t * (1 - d)), right = oracle.at_n(t * (1 + d));
				double error = std::max(std::abs(left.slope / limit.slope - 1), std::abs(right.slope / limit.slope - 1));
				require(error < previous, "S3 one-sided slope fails to approach common limit");
				previous = error;
				log << "one_sided_slope\t" << t << '\t' << d << '\t' << left.slope << '\t' << right.slope << '\t' << limit.slope << '\n';
			}
		}
		double prior_slope = 1, prior_energy = 1;
		for (double t : {eos.NeutronOnsetBaryonDensityFm3(), eos.MuonOnsetBaryonDensityFm3()})
		{
			for (double n : {std::nextafter(t, 0.), t, std::nextafter(t, 1.)})
			{
				auto p = eos.BarotropeAt(n);
				require(p.values_resolved && p.energy_density_MeV_fm3 > 0 && p.pressure_MeV_fm3 > 0, "S3 adjacent values unavailable");
				log << "adjacent_threshold\t" << n << '\t' << p.pressure_MeV_fm3 << '\t' << int(p.domain) << '\n';
			}
		}
		for (size_t grid : {1024, 2048, 4096, 8192})
		{
			auto tab = generate(eos, dir / ("eos-" + std::to_string(grid) + ".tsv"), grid);
			size_t refused_rows = 0;
			for (const auto &p : tab.rows)
			{
				if (p.domain == CompactStar::FreeGasEquilibriumDomain::Npe && p.n_B_fm3 < eos.NeutronOnsetBaryonDensityFm3() * (1 + 1e-7))
				{
					bool failed = false;
					try
					{
						(void)eos.EquilibriumAt(p.n_B_fm3);
					}
					catch (const CompactStar::EquilibriumResolutionError &)
					{
						failed = true;
					}
					require(failed && p.values_resolved, "P7 table rows must not weaken N1");
					++refused_rows;
				}
			}
			require(refused_rows > 0, "P7 table generation omitted refusal window");
			CompactStar::Core::TOVSolver s;
			s.ImportEOS(tab.file.string(), true);
			double emax = 0, smax = 0, nmax = 0, ymax = 0, identity = 0;
			for (size_t j = 1; j < tab.rows.size(); ++j)
			{
				double n = std::sqrt(tab.rows[j - 1].n_B_fm3 * tab.rows[j].n_B_fm3);
				auto q = oracle.at_n(n);
				require(q.e > 0 && q.p > 0 && q.slope >= 3, "S4 whole-domain causal quadrature");
				double pc = q.p * CompactStar::Units::MEV_FM3_TO_ERG_CM3;
				emax = std::max(emax, std::abs(s.GetEDens(pc) / rho(q.e) - 1));
				nmax = std::max(nmax, std::abs(s.GetRho(pc) / q.n - 1));
				smax = std::max(smax, std::abs(s.GetEDensDeriv(pc) / q.slope - 1));
				auto ys = s.GetRho_i(pc);
				identity = std::max({identity, std::abs(ys[0] + ys[1] - 1), std::abs(ys[1] - ys[2] - ys[3])});
				for (int k = 0; k < 4; ++k)
					ymax = std::max(ymax, std::abs(ys[k] - q.ns[k] / q.n));
			}
			require(smax < prior_slope && emax < prior_energy, "S2 nested interpolant error does not decrease");
			if (grid == 8192)
				require(smax < .006 && emax < 1e-6 && nmax < 1e-6 && ymax < 1e-7 && identity < 1e-7, "S2/S4 finest off-grid envelope");
			prior_slope = smax;
			prior_energy = emax;
			log << "interpolant\t" << grid << '\t' << emax << '\t' << nmax << '\t' << ymax << '\t' << smax << '\t' << identity << std::endl;
		}
		std::cout << "S1-S4 PASS; nested errors in local.tsv\n";
		return 0;
	}
	catch (const std::exception &e)
	{
		std::cerr << "BAROTROPE FAILURE: " << e.what() << '\n';
		return 1;
	}
}
