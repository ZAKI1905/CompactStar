// Candidate Structure-1 run artifact generator. Not a scientific baseline.
#pragma once
#include <CompactStar/EOS/TrackRFreeGasThermodynamics.hpp>
#include <CompactStar/Units.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <gsl/gsl_const_cgsm.h>
#include <iomanip>
#include <set>
#include <stdexcept>
#include <vector>
namespace structure1
{
using CompactStar::FreeGasBarotropePoint;
using CompactStar::TrackRFreeGasThermodynamicProvider;
inline void require(bool b, const char *s)
{
	if (!b)
		throw std::runtime_error(s);
}
inline double rho(double e)
{
	const double c = GSL_CONST_CGSM_SPEED_OF_LIGHT;
	return e * CompactStar::Units::MEV_FM3_TO_ERG_CM3 / (c * c);
}
inline double central_n(const TrackRFreeGasThermodynamicProvider &eos, double r)
{
	double a = .5, b = std::nextafter(eos.SigmaMinusOnsetBaryonDensityFm3(), 0.);
	require(r > rho(eos.BarotropeAt(a).energy_density_MeV_fm3) &&
				r < rho(eos.BarotropeAt(b).energy_density_MeV_fm3),
			"central rho outside high-density source domain");
	for (int j = 0; j < 80; ++j)
	{
		double m = (a + b) / 2;
		if (m == a || m == b)
			break;
		if (rho(eos.BarotropeAt(m).energy_density_MeV_fm3) < r)
			a = m;
		else
			b = m;
	}
	return (a + b) / 2;
}
struct Table
{
	std::vector<FreeGasBarotropePoint> rows;
	std::filesystem::path file;
	size_t resolution;
};
inline Table generate(const TrackRFreeGasThermodynamicProvider &eos,
					  const std::filesystem::path &file, size_t resolution, double lower = 1e-14, double upper = .616)
{
	require(resolution >= 64 && lower > 0 && lower < upper && upper < eos.SigmaMinusOnsetBaryonDensityFm3(), "invalid table domain/resolution");
	std::set<double> grid;
	for (size_t i = 0; i <= resolution; ++i)
		grid.insert(lower * std::pow(upper / lower, double(i) / resolution));
	grid.insert(lower);
	grid.insert(upper);
	for (double t : {eos.NeutronOnsetBaryonDensityFm3(), eos.MuonOnsetBaryonDensityFm3()})
	{
		if (t < lower || t > upper)
			continue;
		grid.insert(t);
		// Nested logarithmic distance grids; no interval-wide threshold smoothing.
		size_t per_decade = resolution / 128;
		int decades = t == eos.NeutronOnsetBaryonDensityFm3() ? 8 + int(std::log2(double(resolution) / 1024)) : 8;
		for (size_t i = per_decade; i <= size_t(decades) * per_decade; ++i)
		{
			double d = std::pow(10., -double(i) / per_decade);
			for (double n : {t * (1 - d), t * (1 + d)})
				if (n > lower && n < upper)
					grid.insert(n);
		}
	}
	Table out{{}, file, resolution};
	for (double n : grid)
	{
		if (n > upper)
			continue;
		auto p = eos.BarotropeAt(n);
		if (!out.rows.empty())
			require(p.pressure_MeV_fm3 > out.rows.back().pressure_MeV_fm3 && p.energy_density_MeV_fm3 > out.rows.back().energy_density_MeV_fm3, "unordered EOS values");
		out.rows.push_back(p);
	}
	std::filesystem::create_directories(file.parent_path());
	require(!std::filesystem::exists(file), "refusing to overwrite an EOS artifact");
	std::ofstream f(file);
	f << std::setprecision(17);
	f << "eps(g/cm^3)\tpre(dyn/cm^2)\trho(1/fm^3)\t10\t11\t0\t1\n";
	for (const auto &p : out.rows)
	{
		f << rho(p.energy_density_MeV_fm3) << '\t' << p.pressure_MeV_fm3 * CompactStar::Units::MEV_FM3_TO_ERG_CM3 << '\t' << p.n_B_fm3;
		for (double n : p.number_densities_fm3)
			f << '\t' << n / p.n_B_fm3;
		f << '\n';
	}
	f.close();
	require(bool(f), "table write failed");
	std::ofstream meta(file.string() + ".provenance.txt");
	meta << std::setprecision(17) << "generator=structure1-candidate-v1\nsource_commit=f4ae22d971e25bdd74530aec184f3fe0c3440b95\nprovider_commit=933494d86daf2cf8965079ece49fabd66d9390e5\nmodel=track-r-fernandez-reisenegger-2005-free-gas-local\nresolution=" << resolution << "\nrows=" << out.rows.size() << "\nfloor_n=" << lower << "\nupper_n=" << upper << "\nneutron_onset=" << eos.NeutronOnsetBaryonDensityFm3() << "\nmuon_onset=" << eos.MuonOnsetBaryonDensityFm3() << "\nSigma_ceiling=" << eos.SigmaMinusOnsetBaryonDensityFm3() << "\nMeV_fm3_to_erg_cm3=" << CompactStar::Units::MEV_FM3_TO_ERG_CM3 << "\nc_cgs=" << GSL_CONST_CGSM_SPEED_OF_LIGHT << "\n";
	return out;
}
} // namespace structure1
