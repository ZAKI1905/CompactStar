// -------------------------------------------------------------
//                      CompOSE_Thermo Class
// -------------------------------------------------------------
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/EOS/CompOSE_Thermo.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace CompactStar
{
namespace EOS
{

namespace
{
// -------------------------------------------------------------
inline std::string JoinPath_(const std::string &dir, const std::string &file)
{
	if (dir.empty())
		return file;
	if (dir.back() == '/' || dir.back() == '\\')
		return dir + file;
	return dir + "/" + file;
}

// -------------------------------------------------------------
inline double Lerp_(double a, double b, double w)
{
	return (1.0 - w) * a + w * b;
}

// -------------------------------------------------------------
inline bool IsFinite_(double x)
{
	return std::isfinite(x);
}

// -------------------------------------------------------------
/**
 * @brief Read a CompOSE axis file (eos.t, eos.nb, eos.yq).
 *
 * Common CompOSE format:
 *   <dim_integer>
 *   <N_integer>
 *   v1
 *   v2
 *   ...
 *   vN
 *
 * Some distributions may omit the first two integers and provide only values.
 * This routine attempts the CompOSE header parse first; if it fails, it falls
 * back to reading all doubles in the file as the grid.
 */
std::vector<double> ReadCompOSEAxis_(const std::string &path)
{
	std::ifstream in(path);
	if (!in)
		throw std::runtime_error("CompOSE_Thermo: failed to open file: " + path);

	// Slurp tokens as strings to allow robust "try parse" without consuming stream state.
	std::vector<std::string> tokens;
	tokens.reserve(4096);

	std::string line;
	while (std::getline(in, line))
	{
		const auto hash = line.find('#');
		if (hash != std::string::npos)
			line = line.substr(0, hash);

		std::istringstream iss(line);
		std::string tok;
		while (iss >> tok)
			tokens.push_back(tok);
	}

	if (tokens.empty())
		throw std::runtime_error("CompOSE_Thermo: axis file empty: " + path);

	auto try_parse_int = [](const std::string &s, long long &out) -> bool
	{
		char *end = nullptr;
		errno = 0;
		const long long v = std::strtoll(s.c_str(), &end, 10);
		if (errno != 0 || end == s.c_str() || *end != '\0')
			return false;
		out = v;
		return true;
	};

	auto try_parse_double = [](const std::string &s, double &out) -> bool
	{
		char *end = nullptr;
		errno = 0;
		const double v = std::strtod(s.c_str(), &end);
		if (errno != 0 || end == s.c_str() || *end != '\0')
			return false;
		out = v;
		return true;
	};

	// Attempt CompOSE header interpretation:
	// tokens[0] = dim (ignored), tokens[1] = N
	long long dim = 0;
	long long N = 0;
	bool header_ok = false;
	if (tokens.size() >= 2 && try_parse_int(tokens[0], dim) && try_parse_int(tokens[1], N) && N > 0)
	{
		// Need at least N doubles following.
		if (static_cast<std::size_t>(2 + N) <= tokens.size())
		{
			std::vector<double> grid;
			grid.reserve(static_cast<std::size_t>(N));

			bool ok = true;
			for (long long i = 0; i < N; ++i)
			{
				double x = 0.0;
				if (!try_parse_double(tokens[static_cast<std::size_t>(2 + i)], x))
				{
					ok = false;
					break;
				}
				grid.push_back(x);
			}

			if (ok)
			{
				header_ok = true;

				// Verify monotone non-decreasing
				for (std::size_t i = 1; i < grid.size(); ++i)
				{
					if (grid[i] + 1e-15 < grid[i - 1])
						throw std::runtime_error("CompOSE_Thermo: grid not monotone increasing in: " + path);
				}
				return grid;
			}
		}
	}

	// Fallback: interpret every token as a double and use all successfully parsed values.
	// This tolerates plain column files without headers.
	std::vector<double> vals;
	vals.reserve(tokens.size());
	for (const auto &tok : tokens)
	{
		double x = 0.0;
		if (try_parse_double(tok, x))
			vals.push_back(x);
	}

	if (vals.size() < 2)
		throw std::runtime_error("CompOSE_Thermo: could not parse a valid axis grid from: " + path);

	for (std::size_t i = 1; i < vals.size(); ++i)
	{
		if (vals[i] + 1e-15 < vals[i - 1])
			throw std::runtime_error("CompOSE_Thermo: grid not monotone increasing in: " + path);
	}

	(void)header_ok;
	return vals;
}

// -------------------------------------------------------------
/**
 * @brief Bilinear interpolation on a rectilinear grid.
 *
 * z is stored row-major: z[ix * Ny + iy], ix indexes x-grid, iy indexes y-grid.
 */
double Bilinear_(const std::vector<double> &xgrid,
				 const std::vector<double> &ygrid,
				 const std::vector<double> &z,
				 double x, double y)
{
	const std::size_t Nx = xgrid.size();
	const std::size_t Ny = ygrid.size();
	if (Nx < 2 || Ny < 2)
		throw std::runtime_error("CompOSE_Thermo: Bilinear requires Nx,Ny >= 2.");

	auto ix_it = std::upper_bound(xgrid.begin(), xgrid.end(), x);
	std::size_t ix = 0;
	if (ix_it == xgrid.begin())
		ix = 0;
	else if (ix_it == xgrid.end())
		ix = Nx - 2;
	else
		ix = static_cast<std::size_t>(ix_it - xgrid.begin() - 1);

	auto iy_it = std::upper_bound(ygrid.begin(), ygrid.end(), y);
	std::size_t iy = 0;
	if (iy_it == ygrid.begin())
		iy = 0;
	else if (iy_it == ygrid.end())
		iy = Ny - 2;
	else
		iy = static_cast<std::size_t>(iy_it - ygrid.begin() - 1);

	const double x0 = xgrid[ix];
	const double x1 = xgrid[ix + 1];
	const double y0 = ygrid[iy];
	const double y1 = ygrid[iy + 1];

	const double tx = (std::abs(x1 - x0) > 0.0) ? (x - x0) / (x1 - x0) : 0.0;
	const double ty = (std::abs(y1 - y0) > 0.0) ? (y - y0) / (y1 - y0) : 0.0;

	const auto idx = [&](std::size_t a, std::size_t b) -> std::size_t
	{
		return a * Ny + b;
	};

	const double z00 = z[idx(ix, iy)];
	const double z10 = z[idx(ix + 1, iy)];
	const double z01 = z[idx(ix, iy + 1)];
	const double z11 = z[idx(ix + 1, iy + 1)];

	const double z0 = Lerp_(z00, z10, tx);
	const double z1 = Lerp_(z01, z11, tx);
	return Lerp_(z0, z1, ty);
}

} // anonymous namespace

// -------------------------------------------------------------
CompOSE_Thermo::CompOSE_Thermo(const std::string &directory, const Options &opt)
{
	LoadFromDirectory(directory, opt);
}

// -------------------------------------------------------------
void CompOSE_Thermo::LoadFromDirectory(const std::string &directory, const Options &opt)
{
	m_loaded = false;
	m_opt = opt;

	m_T.clear();
	m_nb.clear();
	m_Yq.clear();
	m_Q2_plane_values.clear();

	ReadAxes_(directory);
	ReadThermoQ2_(directory);

	m_loaded = true;
}

// -------------------------------------------------------------
void CompOSE_Thermo::ReadAxes_(const std::string &directory)
{
	const std::string t_path = JoinPath_(directory, "eos.t");
	const std::string nb_path = JoinPath_(directory, "eos.nb");
	const std::string yq_path = JoinPath_(directory, "eos.yq");

	m_T = ReadCompOSEAxis_(t_path);
	m_nb = ReadCompOSEAxis_(nb_path);
	m_Yq = ReadCompOSEAxis_(yq_path);

	if (m_T.size() < 2)
		throw std::runtime_error("CompOSE_Thermo: eos.t needs >= 2 points.");
	if (m_nb.size() < 2)
		throw std::runtime_error("CompOSE_Thermo: eos.nb needs >= 2 points.");
	if (m_Yq.size() < 2)
		throw std::runtime_error("CompOSE_Thermo: eos.yq needs >= 2 points.");
}

// -------------------------------------------------------------
void CompOSE_Thermo::ReadThermoQ2_(const std::string &directory)
{
	const std::string thermo_path = JoinPath_(directory, "eos.thermo");
	std::ifstream in(thermo_path);
	if (!in)
		throw std::runtime_error("CompOSE_Thermo: failed to open file: " + thermo_path);

	// Header: mn mp Il
	{
		std::string header;
		if (!std::getline(in, header))
			throw std::runtime_error("CompOSE_Thermo: eos.thermo is empty.");

		std::istringstream iss(header);
		if (!(iss >> m_mn_MeV >> m_mp_MeV >> m_Il))
			throw std::runtime_error("CompOSE_Thermo: failed to parse eos.thermo header (mn mp Il).");
	}

	const std::size_t NT = m_T.size();
	const std::size_t NNb = m_nb.size();
	const std::size_t NYq = m_Yq.size();

	const double NaN = std::numeric_limits<double>::quiet_NaN();
	m_Q2_plane_values.assign(NT, std::vector<double>(NNb * NYq, NaN));

	std::string line;
	std::size_t data_lines = 0;

	while (std::getline(in, line))
	{
		// Skip empty/whitespace-only lines
		bool all_ws = true;
		for (char c : line)
		{
			if (!std::isspace(static_cast<unsigned char>(c)))
			{
				all_ws = false;
				break;
			}
		}
		if (all_ws)
			continue;

		std::istringstream iss(line);

		int iT = 0, inb_i = 0, iYq = 0;
		double Q1 = 0.0, Q2 = 0.0, Q3 = 0.0, Q4 = 0.0, Q5 = 0.0, Q6 = 0.0, Q7 = 0.0;
		int Nadd = 0;

		if (!(iss >> iT >> inb_i >> iYq >> Q1 >> Q2 >> Q3 >> Q4 >> Q5 >> Q6 >> Q7 >> Nadd))
		{
			throw std::runtime_error("CompOSE_Thermo: parse error in eos.thermo at data line " + std::to_string(data_lines + 2) + ".");
		}

		// Skip optional columns using Nadd
		for (int k = 0; k < Nadd; ++k)
		{
			double dummy = 0.0;
			if (!(iss >> dummy))
			{
				throw std::runtime_error("CompOSE_Thermo: failed to read optional column " + std::to_string(k + 1) + "/" + std::to_string(Nadd) + " at eos.thermo line " + std::to_string(data_lines + 2) + ".");
			}
		}

		// CompOSE indices are 1-based
		if (iT <= 0 || inb_i <= 0 || iYq <= 0)
			throw std::runtime_error("CompOSE_Thermo: invalid 1-based indices in eos.thermo.");

		const std::size_t t_idx = static_cast<std::size_t>(iT - 1);
		const std::size_t nb_idx = static_cast<std::size_t>(inb_i - 1);
		const std::size_t yq_idx = static_cast<std::size_t>(iYq - 1);

		if (t_idx >= NT || nb_idx >= NNb || yq_idx >= NYq)
		{
			throw std::runtime_error("CompOSE_Thermo: eos.thermo index out of bounds at data line " + std::to_string(data_lines + 2) + ".");
		}

		m_Q2_plane_values[t_idx][nb_idx * NYq + yq_idx] = Q2;
		++data_lines;
	}

	// Validate completeness (strict, but catches partial/corrupt tables immediately).
	std::size_t missing = 0;
	for (std::size_t t = 0; t < NT; ++t)
	{
		for (std::size_t idx = 0; idx < NNb * NYq; ++idx)
		{
			if (!IsFinite_(m_Q2_plane_values[t][idx]))
				++missing;
		}
	}
	if (missing != 0)
	{
		throw std::runtime_error("CompOSE_Thermo: eos.thermo Q2 grid incomplete; missing " + std::to_string(missing) + " points.");
	}
}

// -------------------------------------------------------------
std::size_t CompOSE_Thermo::BracketIndex_(const std::vector<double> &grid, double x)
{
	if (grid.size() < 2)
		throw std::runtime_error("CompOSE_Thermo: BracketIndex_ requires grid.size() >= 2.");

	auto it = std::upper_bound(grid.begin(), grid.end(), x);
	if (it == grid.begin())
		return 0;
	if (it == grid.end())
		return grid.size() - 2;
	return static_cast<std::size_t>(it - grid.begin() - 1);
}

// -------------------------------------------------------------
void CompOSE_Thermo::ClampToDomain_(double &T_MeV, double &nb_fm3, double &Yq) const
{
	const auto clamp = [](double &v, double lo, double hi)
	{
		if (v < lo)
			v = lo;
		if (v > hi)
			v = hi;
	};

	clamp(T_MeV, m_T.front(), m_T.back());
	clamp(nb_fm3, m_nb.front(), m_nb.back());
	clamp(Yq, m_Yq.front(), m_Yq.back());
}

// -------------------------------------------------------------
double CompOSE_Thermo::Q2_OnPlane(std::size_t iT, double nb_fm3, double Yq) const
{
	if (!m_loaded)
		throw std::runtime_error("CompOSE_Thermo::Q2_OnPlane: table not loaded.");
	if (iT >= m_T.size())
		throw std::runtime_error("CompOSE_Thermo::Q2_OnPlane: iT out of range.");

	return Bilinear_(m_nb, m_Yq, m_Q2_plane_values[iT], nb_fm3, Yq);
}

// -------------------------------------------------------------
double CompOSE_Thermo::Q2(double T_MeV, double nb_fm3, double Yq) const
{
	if (!m_loaded)
		throw std::runtime_error("CompOSE_Thermo::Q2: table not loaded.");

	if (m_opt.clamp_to_domain)
	{
		ClampToDomain_(T_MeV, nb_fm3, Yq);
	}
	else
	{
		if (T_MeV < m_T.front() || T_MeV > m_T.back() ||
			nb_fm3 < m_nb.front() || nb_fm3 > m_nb.back() ||
			Yq < m_Yq.front() || Yq > m_Yq.back())
			throw std::runtime_error("CompOSE_Thermo::Q2: query out of domain.");
	}

	const std::size_t i = BracketIndex_(m_T, T_MeV);
	const double T0 = m_T[i];
	const double T1 = m_T[i + 1];

	const double q0 = Q2_OnPlane(i, nb_fm3, Yq);
	const double q1 = Q2_OnPlane(i + 1, nb_fm3, Yq);

	const double w = (std::abs(T1 - T0) > 0.0) ? (T_MeV - T0) / (T1 - T0) : 0.0;
	return Lerp_(q0, q1, w);
}

// -------------------------------------------------------------
double CompOSE_Thermo::Q2_ForCooling(double T_MeV, double nb_fm3, double Yq) const
{
	if (!m_loaded)
		throw std::runtime_error("CompOSE_Thermo::Q2_ForCooling: table not loaded.");

	if (m_opt.clamp_to_domain)
	{
		ClampToDomain_(T_MeV, nb_fm3, Yq);
	}
	else
	{
		if (T_MeV < m_T.front() || T_MeV > m_T.back() ||
			nb_fm3 < m_nb.front() || nb_fm3 > m_nb.back() ||
			Yq < m_Yq.front() || Yq > m_Yq.back())
			throw std::runtime_error("CompOSE_Thermo::Q2_ForCooling: query out of domain.");
	}

	// If low-T modeling is disabled, fall back to table interpolation
	if (!m_opt.enable_lowT_fit)
		return Q2(T_MeV, nb_fm3, Yq);

	const double Tswitch = m_opt.lowT_switch_MeV;
	const double wblend = m_opt.lowT_blend_width_MeV;

	// Low-T fitted slope: Q2 ~ a0 * T
	const double a0 = LowT_Slope_dQ2dT0_(nb_fm3, Yq);
	const double Q2_lowT = a0 * T_MeV;

	// Table-based entropy
	const double Q2_tab = Q2(T_MeV, nb_fm3, Yq);

	// Blend between low-T model and table entropy
	const double w = BlendLowT_(T_MeV, Tswitch, Tswitch, wblend);
	return (1.0 - w) * Q2_lowT + w * Q2_tab;
}

// -------------------------------------------------------------
double CompOSE_Thermo::dQ2dT(double T_MeV, double nb_fm3, double Yq) const
{
	if (!m_loaded)
		throw std::runtime_error("CompOSE_Thermo::dQ2dT: table not loaded.");

	if (m_opt.clamp_to_domain)
	{
		ClampToDomain_(T_MeV, nb_fm3, Yq);
	}
	else
	{
		if (T_MeV < m_T.front() || T_MeV > m_T.back() ||
			nb_fm3 < m_nb.front() || nb_fm3 > m_nb.back() ||
			Yq < m_Yq.front() || Yq > m_Yq.back())
			throw std::runtime_error("CompOSE_Thermo::dQ2dT: query out of domain.");
	}

	const std::size_t NT = m_T.size();

	auto deriv_at_index = [&](std::size_t iT) -> double
	{
		// forward at lower boundary
		if (iT == 0)
		{
			const double T0 = m_T[0], T1 = m_T[1];
			const double q0 = Q2_OnPlane(0, nb_fm3, Yq);
			const double q1 = Q2_OnPlane(1, nb_fm3, Yq);
			return (q1 - q0) / (T1 - T0);
		}

		// backward at upper boundary
		if (iT >= NT - 1)
		{
			const double Tm1 = m_T[NT - 2], Tm0 = m_T[NT - 1];
			const double qm1 = Q2_OnPlane(NT - 2, nb_fm3, Yq);
			const double qm0 = Q2_OnPlane(NT - 1, nb_fm3, Yq);
			return (qm0 - qm1) / (Tm0 - Tm1);
		}

		// central in interior
		if (m_opt.use_central_difference && iT >= 1 && (iT + 1) < NT)
		{
			const double Ta = m_T[iT - 1], Tb = m_T[iT + 1];
			const double qa = Q2_OnPlane(iT - 1, nb_fm3, Yq);
			const double qb = Q2_OnPlane(iT + 1, nb_fm3, Yq);
			return (qb - qa) / (Tb - Ta);
		}

		// fallback: local forward
		const double T0 = m_T[iT], T1 = m_T[iT + 1];
		const double q0 = Q2_OnPlane(iT, nb_fm3, Yq);
		const double q1 = Q2_OnPlane(iT + 1, nb_fm3, Yq);
		return (q1 - q0) / (T1 - T0);
	};

	// Optional low-T anchoring: if caller requests, use the first-interval derivative for very small T.
	if (m_opt.Tmin_for_derivative_MeV > 0.0 && T_MeV <= m_opt.Tmin_for_derivative_MeV)
		return deriv_at_index(0);

	// Smooth derivative across intervals by interpolating between node derivatives.
	const std::size_t i = BracketIndex_(m_T, T_MeV);
	const double T0 = m_T[i];
	const double T1 = m_T[i + 1];
	const double w = (std::abs(T1 - T0) > 0.0) ? (T_MeV - T0) / (T1 - T0) : 0.0;

	const double d0 = deriv_at_index(i);
	const double d1 = deriv_at_index(i + 1);
	return Lerp_(d0, d1, w);
}

// -------------------------------------------------------------
// Natural units: d e / dT(MeV)  -> fm^-3
double CompOSE_Thermo::CvDensity_Natural(double T_MeV, double nb_fm3, double Yq) const
{
	return T_MeV * nb_fm3 * dQ2dT(T_MeV, nb_fm3, Yq);
}

// -------------------------------------------------------------
// cV in cgs: erg cm^-3 K^-1
double CompOSE_Thermo::CvDensity_cgs(double T_MeV, double nb_fm3, double Yq) const
{
	const double cv_nat_fm3 = CvDensity_Natural(T_MeV, nb_fm3, Yq); // fm^-3

	// kB in MeV/K (since T[MeV] = kB*T[K])
	constexpr double kB_MeV_per_K = 8.617333262e-11;

	// 1 (MeV fm^-3) = 1.602176634e33 (erg cm^-3)   [CompOSE manual]
	constexpr double MeVfm3_to_ergcm3 = 1.602176634e33;

	// (fm^-3) * (MeV/K) -> (MeV fm^-3 / K) -> erg cm^-3 / K
	return cv_nat_fm3 *
		   kB_MeV_per_K * MeVfm3_to_ergcm3;
}

// -------------------------------------------------------------
double CompOSE_Thermo::CvPerBaryon(double T_MeV, double nb_fm3, double Yq) const
{
	(void)nb_fm3; // kept for consistent signature + domain checks in dQ2dT
	return T_MeV * dQ2dT(T_MeV, nb_fm3, Yq);
}

// -------------------------------------------------------------
//------------------------------------------------------------------------------
// Smooth blending helpers
double CompOSE_Thermo::SmoothStep01_(double x)
{
	// Clamp to [0,1]
	if (x <= 0.0)
		return 0.0;
	if (x >= 1.0)
		return 1.0;
	// Classic smoothstep: 3x^2 - 2x^3 (C1 continuous)
	return x * x * (3.0 - 2.0 * x);
}

double CompOSE_Thermo::BlendLowT_(double T, double lowT, double highT, double w) const
{
	// Return blend weight in [0,1] that goes from 0 (use low-T model) to 1 (use high-T model)
	// over a window centered on [lowT, highT]. Here we use a symmetric blend around lowT_switch.
	(void)highT; // kept for readability if you later change policy
	if (w <= 0.0)
		return (T >= lowT ? 1.0 : 0.0);

	const double T0 = lowT - w;
	const double T1 = lowT + w;
	if (T <= T0)
		return 0.0;
	if (T >= T1)
		return 1.0;
	const double x = (T - T0) / (T1 - T0);
	return SmoothStep01_(x);
}

//------------------------------------------------------------------------------
// Low-T slope inference using constrained least squares fit Q2(T) ~ a*T
//
// We enforce Q2(0)=0 by construction and fit only the slope a using the first
// few temperature grid points. This is robust and avoids relying solely on the
// first interval [0,2] MeV.
//
// Given data (Ti, Qi), minimize sum_i (Qi - a*Ti)^2 =>
//   a = (sum Ti*Qi)/(sum Ti^2), excluding the T=0 point (since it contributes nothing).
//
double CompOSE_Thermo::LowT_Slope_dQ2dT0_(double nb_fm3, double Yq) const
{
	if (!m_loaded)
		throw std::runtime_error("CompOSE_Thermo::LowT_Slope_dQ2dT0_: table not loaded.");

	// Need at least two temperature points
	const std::size_t NT = m_T.size();
	if (NT < 2)
		throw std::runtime_error("CompOSE_Thermo::LowT_Slope_dQ2dT0_: NT < 2.");

	// Determine how many points to use
	int nfit = m_opt.lowT_fit_points;
	if (nfit < 2)
		nfit = 2;
	if (static_cast<std::size_t>(nfit) > NT)
		nfit = static_cast<int>(NT);

	// Accumulate least-squares slope a = (Σ Ti*Qi)/(Σ Ti^2), skipping Ti=0.
	double sum_TQ = 0.0;
	double sum_T2 = 0.0;

	for (int i = 0; i < nfit; ++i)
	{
		const double Ti = m_T[static_cast<std::size_t>(i)];
		if (Ti <= 0.0)
			continue; // skip T=0 point

		const double Qi = Q2_OnPlane(static_cast<std::size_t>(i), nb_fm3, Yq);

		sum_TQ += Ti * Qi;
		sum_T2 += Ti * Ti;
	}

	// Fallback: if fit is ill-posed (e.g., all Ti=0), revert to first-interval slope.
	if (sum_T2 <= 0.0)
	{
		const double T0 = m_T[0];
		const double T1 = m_T[1];
		const double Q0 = Q2_OnPlane(0, nb_fm3, Yq);
		const double Q1 = Q2_OnPlane(1, nb_fm3, Yq);
		return (Q1 - Q0) / (T1 - T0);
	}

	return sum_TQ / sum_T2; // units 1/MeV
}

//------------------------------------------------------------------------------
// Public low-T cooling derivative
double CompOSE_Thermo::dQ2dT_ForCooling(double T_MeV, double nb_fm3, double Yq) const
{
	if (!m_loaded)
		throw std::runtime_error("CompOSE_Thermo::dQ2dT_ForCooling: table not loaded.");

	if (m_opt.clamp_to_domain)
	{
		ClampToDomain_(T_MeV, nb_fm3, Yq);
	}
	else
	{
		if (T_MeV < m_T.front() || T_MeV > m_T.back() ||
			nb_fm3 < m_nb.front() || nb_fm3 > m_nb.back() ||
			Yq < m_Yq.front() || Yq > m_Yq.back())
			throw std::runtime_error("CompOSE_Thermo::dQ2dT_ForCooling: query out of domain.");
	}

	// If low-T fit disabled, just use the regular derivative
	if (!m_opt.enable_lowT_fit)
		return dQ2dT(T_MeV, nb_fm3, Yq);

	const double Tswitch = m_opt.lowT_switch_MeV;
	const double wblend = m_opt.lowT_blend_width_MeV;

	// Fit-based slope (T->0) for this (nB,Yq)
	const double a0 = LowT_Slope_dQ2dT0_(nb_fm3, Yq);

	// Table-based derivative at this T
	const double dtab = dQ2dT(T_MeV, nb_fm3, Yq);

	// If below the switch region, use the fitted slope; above, use table;
	// blend smoothly around Tswitch to avoid kinks.
	const double w = BlendLowT_(T_MeV, Tswitch, Tswitch, wblend);
	return (1.0 - w) * a0 + w * dtab;
}

//------------------------------------------------------------------------------
// Natural units: d e / dT(MeV)  -> fm^-3
double CompOSE_Thermo::CvDensity_Natural_ForCooling(double T_MeV, double nb_fm3, double Yq) const
{
	// c_V = T * nB * dQ2/dT
	return T_MeV * nb_fm3 * dQ2dT_ForCooling(T_MeV, nb_fm3, Yq);
}
// -------------------------------------------------------------
// cV in cgs: erg cm^-3 K^-1
double CompOSE_Thermo::CvDensity_cgs_ForCooling(double T_MeV, double nb_fm3, double Yq) const
{
	const double cv_nat_fm3 = CvDensity_Natural_ForCooling(T_MeV, nb_fm3, Yq); // fm^-3

	// kB in MeV/K (since T[MeV] = kB*T[K])
	constexpr double kB_MeV_per_K = 8.617333262e-11;

	// 1 (MeV fm^-3) = 1.602176634e33 (erg cm^-3)   [CompOSE manual]
	constexpr double MeVfm3_to_ergcm3 = 1.602176634e33;

	// (fm^-3) * (MeV/K) -> (MeV fm^-3 / K) -> erg cm^-3 / K
	return cv_nat_fm3 *
		   kB_MeV_per_K * MeVfm3_to_ergcm3;
}

//------------------------------------------------------------------------------

} // namespace EOS
} // namespace CompactStar