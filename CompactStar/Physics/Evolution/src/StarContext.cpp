// -*- lsst-c++ -*-
/**
 * @file StarContext.cpp
 * @brief See StarContext.hpp. Read-only adapter over StarProfile.
 */

#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"

#include <Zaki/Physics/Constants.hpp> // unit conversions
#include <Zaki/Vector/DataColumn.hpp>

#include <cmath>
#include <stdexcept>

//==============================================================
namespace
{
inline double Lerp(double a, double b, double w) { return (1.0 - w) * a + w * b; }

inline std::vector<double> LogSpace(double a, double b, std::size_t N)
{
	if (N < 2 || a <= 0.0 || b <= 0.0)
		throw std::runtime_error("LogSpace: invalid args.");
	const double la = std::log(a);
	const double lb = std::log(b);
	std::vector<double> out(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double w = double(i) / double(N - 1);
		out[i] = std::exp((1.0 - w) * la + w * lb);
	}
	return out;
}

inline std::size_t Bracket(const std::vector<double> &grid, double x, std::size_t hint)
{
	if (grid.size() < 2)
		throw std::runtime_error("Bracket: grid too small.");

	if (hint + 1 < grid.size() && grid[hint] <= x && x <= grid[hint + 1])
		return hint;

	auto it = std::upper_bound(grid.begin(), grid.end(), x);
	if (it == grid.begin())
		return 0;
	if (it == grid.end())
		return grid.size() - 2;
	return std::size_t(it - grid.begin() - 1);
}

// CompOSE strong-sector charge q/e for numeric species codes.
inline bool ComposeCharge(int code, double &q)
{
	switch (code)
	{
	case 10:
		q = 0.0;
		return true; // n
	case 11:
		q = +1.0;
		return true; // p
	case 20:
		q = -1.0;
		return true; // Δ-
	case 21:
		q = 0.0;
		return true; // Δ0
	case 22:
		q = +1.0;
		return true; // Δ+
	case 23:
		q = +2.0;
		return true; // Δ++
	case 100:
		q = 0.0;
		return true; // Λ0
	case 110:
		q = -1.0;
		return true; // Σ-
	case 111:
		q = 0.0;
		return true; // Σ0
	case 112:
		q = +1.0;
		return true; // Σ+
	case 120:
		q = -1.0;
		return true; // Ξ-
	case 121:
		q = 0.0;
		return true; // Ξ0
	case 500:
		q = +2.0 / 3.0;
		return true; // u
	case 501:
		q = -1.0 / 3.0;
		return true; // d
	case 502:
		q = -1.0 / 3.0;
		return true; // s
	default:
		return false;
	}
}
} // namespace
//==============================================================

namespace CompactStar::Physics::Evolution
{
// --------------------
// Helpers
// --------------------
std::uint64_t StarContext::ProfileVersion() const
{
	if (!m_prof)
		return 0;

	// Examples: m_prof->Version(), m_prof->GetVersion(), m_prof->version().
	return static_cast<std::uint64_t>(m_prof->Version());
}
//==============================================================
//                   StarContext Class
//==============================================================

//--------------------------------------------------------------
StarContext::StarContext(const CompactStar::Core::StarProfile &prof)
	: m_prof(&prof)
{
	BindColumnsOrThrow_();
	ValidateOrThrow_();

	// Initialize cache version snapshot to current (cache is empty until requested).
	m_cached_version = ProfileVersion();
}

//--------------------------------------------------------------
std::size_t StarContext::Size() const
{
	return (m_r ? m_r->Size() : 0);
}

//--------------------------------------------------------------
double StarContext::RadiusSurface() const
{
	if (!m_r || m_r->Size() == 0)
		return 0.0;

	return (*m_r)[-1]; // km
}

//--------------------------------------------------------------
double StarContext::MassSurface() const
{
	if (!m_m || m_m->Size() == 0)
		return 0.0;

	return (*m_m)[-1]; // km
}

//--------------------------------------------------------------
double StarContext::ExpNuSurface() const
{
	if (!m_nu || m_nu->Size() == 0)
		return 0.0;

	return std::exp((*m_nu)[-1]);
}

//==============================================================
//                   Private helpers
//==============================================================

//--------------------------------------------------------------
void StarContext::BindColumnsOrThrow_()
{
	// Mandatory geometry
	m_r = m_prof->GetRadius();
	m_m = m_prof->GetMass();

	if (!m_r)
		throw std::runtime_error("StarContext: missing radius column r(km)");
	if (!m_m)
		throw std::runtime_error("StarContext: missing mass column m(km)");

	// Optional but strongly expected
	m_nu = m_prof->GetMetricNu();
	m_lam = m_prof->GetMetricLambda();

	// Optional thermodynamics
	m_nb = m_prof->GetBaryonDensity();
	m_pre = m_prof->GetPressure();
	m_eps = m_prof->GetEnergyDensity();
}

//--------------------------------------------------------------
void StarContext::ValidateOrThrow_()
{
	const std::size_t n = m_r->Size();
	if (n == 0)
		throw std::runtime_error("StarContext: profile has zero rows");

	auto check = [&](const Zaki::Vector::DataColumn *c,
					 const char *name)
	{
		if (c && c->Size() != n)
		{
			throw std::runtime_error(
				std::string("StarContext: column '") + name +
				"' has size " + std::to_string(c->Size()) +
				", expected " + std::to_string(n));
		}
	};

	check(m_m, "m");
	check(m_nu, "nu");
	check(m_lam, "lambda");
	check(m_nb, "nB");
	check(m_pre, "p");
	check(m_eps, "eps");
}

//--------------------------------------------------------------
const Zaki::Vector::DataColumn *StarContext::MassDensity_gcm3() const
{
	if (!IsValid())
		return nullptr;

	// Need energy density to build rho
	if (!m_eps)
		return nullptr;

	RefreshDerivedCachesIfNeeded_();

	return m_rho_gcm3.get();
}

//--------------------------------------------------------------
void StarContext::RefreshDerivedCachesIfNeeded_() const
{
	// If the profile is not valid, do not attempt to build derived caches.
	if (!IsValid())
		return;

	const auto v = ProfileVersion();

	// If profile changed since last snapshot, invalidate derived caches.
	if (v != m_cached_version)
	{
		// Invalidate mass density cache
		m_rho_gcm3.reset();

		// Invalidate direct Urca caches
		m_durca_mask.reset();
		m_durca_last_allowed = -1;
		m_durca_boundary_r_km = 0.0;

		// Invalidate derived Yq cache
		m_Yq_cache.reset();

		// Invalidate heat capacity cache
		m_cv_cache = HeatCapacityCache{};

		// Update cached version
		m_cached_version = v;
	}

	// Build on demand
	if (!m_rho_gcm3 && m_eps)
		BuildMassDensityCache_();

	// Build direct-Urca mask only when needed, not every time.
	// Guard on required inputs
	if (!m_durca_mask)
	{
		if (m_nb && m_prof)
			BuildDirectUrcaMaskCache_();
	}
}

//--------------------------------------------------------------
// Build mass density cache in g/cm^3 from energy density eps (km^-2)
void StarContext::BuildMassDensityCache_() const
{
	// Defensive: ensure profile still valid
	if (!m_prof || !m_eps)
		return;

	const std::size_t n = static_cast<std::size_t>(m_eps->Size());
	if (n == 0)
	{
		m_rho_gcm3.reset();
		return;
	}

	// ---- Unit conversion factor ----
	// eps is stored as (km^-2). Convert to rho [g/cm^3].
	const double kKmMinus2_to_gcm3 = Zaki::Physics::MEV_FM3_2_G_CM3 /
									 Zaki::Physics::MEV_FM3_2_INV_KM2;

	// Build values
	std::vector<double> rho(n);
	for (std::size_t i = 0; i < n; ++i)
	{
		const double eps_km2 = (*m_eps)[i]; // adjust if DataColumn uses operator[]
		rho[i] = eps_km2 * kKmMinus2_to_gcm3;
	}

	// Construct column
	m_rho_gcm3 = std::make_unique<Zaki::Vector::DataColumn>("rho(g/cm^3)", rho);
}

// //--------------------------------------------------------------
// // Build direct Urca mask cache
// void StarContext::BuildDirectUrcaMaskCache_() const
// {
// 	// Defensive: ensure profile still valid
// 	if (!m_prof || !m_nb || !m_r)
// 		return;

// 	const std::size_t n = static_cast<std::size_t>(m_nb->Size());
// 	if (n == 0)
// 	{
// 		m_durca_mask.reset();
// 		m_durca_last_allowed = -1;
// 		m_durca_boundary_r_km = 0.0;
// 		return;
// 	}

// 	// Species fractions Yi = n_i / nB are stored as "species" columns.
// 	// Per CompOSE convention:
// 	//   10: neutron, 11: proton, 0: electron
// 	const auto *Yn_col = m_prof->GetSpeciesPtr("10");
// 	const auto *Yp_col = m_prof->GetSpeciesPtr("11");
// 	const auto *Ye_col = m_prof->GetSpeciesPtr("0");

// 	if (!Yn_col || !Yp_col || !Ye_col)
// 	{
// 		// If any required species fraction is missing, we cannot build the mask.
// 		m_durca_mask.reset();
// 		m_durca_last_allowed = -1;
// 		m_durca_boundary_r_km = 0.0;
// 		return;
// 	}

// 	// Sanity check on sizes
// 	if (static_cast<std::size_t>(Yn_col->Size()) != n ||
// 		static_cast<std::size_t>(Yp_col->Size()) != n ||
// 		static_cast<std::size_t>(Ye_col->Size()) != n ||
// 		static_cast<std::size_t>(m_r->Size()) != n)
// 	{
// 		m_durca_mask.reset();
// 		m_durca_last_allowed = -1;
// 		m_durca_boundary_r_km = 0.0;
// 		return;
// 	}

// 	// Store as 0/1 bytes.
// 	std::vector<uint8_t> mask(n, 0);

// 	// kF = (3*pi^2*n)^(1/3).
// 	const double three_pi2 = 3.0 * M_PI * M_PI;

// 	// Compute mask and boundary in ONE reverse scan to avoid a second pass.
// 	long last_allowed = -1;
// 	double boundary_r_km = 0.0;

// 	for (long il = static_cast<long>(n) - 1; il >= 0; --il)
// 	{
// 		const std::size_t i = static_cast<std::size_t>(il);

// 		const double nB = (*m_nb)[i];
// 		const double Yn = (*Yn_col)[i];
// 		const double Yp = (*Yp_col)[i];
// 		const double Ye = (*Ye_col)[i];

// 		// If nB>0 and Y>=0 then nn/np/ne are automatically >=0.
// 		if (!(nB > 0.0) || Yn < 0.0 || Yp < 0.0 || Ye < 0.0)
// 		{
// 			mask[i] = 0;
// 			continue;
// 		}

// 		const double nn = Yn * nB;
// 		const double np = Yp * nB;
// 		const double ne = Ye * nB;

// 		// If any is exactly 0, cbrt handles it fine and DU is not allowed anyway.
// 		const double kFn = std::cbrt(three_pi2 * nn);
// 		const double kFp = std::cbrt(three_pi2 * np);
// 		const double kFe = std::cbrt(three_pi2 * ne);

// 		const bool allowed = (kFn <= (kFp + kFe));
// 		mask[i] = allowed ? 1 : 0;

// 		// First allowed we encounter in reverse order is the outermost allowed point.
// 		if (allowed && last_allowed < 0)
// 		{
// 			last_allowed = il;
// 			boundary_r_km = (*m_r)[i];
// 			// We do NOT break because we still need to fill mask for smaller i.
// 		}
// 	}

// 	m_durca_mask = std::make_unique<std::vector<std::uint8_t>>(std::move(mask));
// 	m_durca_last_allowed = last_allowed;
// 	m_durca_boundary_r_km = boundary_r_km;
// }

//--------------------------------------------------------------
// Build direct Urca mask cache
//
// Direct Urca (DU) kinematic allowance is typically checked via the
// Fermi-momentum triangle condition:
//
//    kFn <= kFp + kFl
//
// where l is a lepton (electron, sometimes muon). In your current setup you
// use electrons only.
//
// IMPORTANT for your EOS layout:
// - index 0 is the core
// - the EOS includes a crust region where electrons exist but baryons
//   (neutrons/protons) may be effectively absent.
//
// If you blindly apply the triangle condition in an electron-only region,
// you get the pathological result:
//
//   nn = np = 0  => kFn = kFp = 0, while kFe > 0  =>  0 <= 0 + kFe  (always true)
//
// which incorrectly flags DU as allowed in the crust.
//
// Fix:
// 1) Require a *baryonic medium* with meaningful nn and np before testing DU.
// 2) Define the DU boundary as the end of the *contiguous core DU-allowed region*
//    as you scan outward from the core. This prevents any later numerical
//    “islands” from corrupting the boundary.
//
// Notes on thresholds:
// - nB_min: excludes crust/vacuum-like zones where baryons are not present.
// - n_min : excludes nn/np ~ 0 where the triangle test becomes meaningless.
//   These thresholds are NOT physics thresholds for DU; they are numerical/semantic
//   guards against applying the criterion outside its domain.
//--------------------------------------------------------------
void StarContext::BuildDirectUrcaMaskCache_() const
{
	// Defensive: ensure profile still valid and required core columns exist.
	if (!m_prof || !m_nb || !m_r)
		return;

	const std::size_t n = static_cast<std::size_t>(m_nb->Size());
	if (n == 0)
	{
		m_durca_mask.reset();
		m_durca_last_allowed = -1;
		m_durca_boundary_r_km = 0.0;
		return;
	}

	// Species fractions Yi = n_i / nB are stored as "species" columns.
	// Per CompOSE convention in your loader:
	//   "10": neutron, "11": proton, "0": electron
	const auto *Yn_col = m_prof->GetSpeciesPtr("10");
	const auto *Yp_col = m_prof->GetSpeciesPtr("11");
	const auto *Ye_col = m_prof->GetSpeciesPtr("0");

	if (!Yn_col || !Yp_col || !Ye_col)
	{
		// Missing species fractions => cannot build DU mask.
		m_durca_mask.reset();
		m_durca_last_allowed = -1;
		m_durca_boundary_r_km = 0.0;
		return;
	}

	// Sanity check: all relevant columns must have consistent length.
	if (static_cast<std::size_t>(Yn_col->Size()) != n ||
		static_cast<std::size_t>(Yp_col->Size()) != n ||
		static_cast<std::size_t>(Ye_col->Size()) != n ||
		static_cast<std::size_t>(m_r->Size()) != n)
	{
		m_durca_mask.reset();
		m_durca_last_allowed = -1;
		m_durca_boundary_r_km = 0.0;
		return;
	}

	// Output mask (0/1) for all radii.
	std::vector<std::uint8_t> mask(n, 0);

	// kF = (3*pi^2*n)^(1/3), with n in fm^-3 -> kF in fm^-1 (consistent).
	const double three_pi2 = 3.0 * M_PI * M_PI;

	// ---------------------------------------------------------------------
	// Guards against the "electron-only crust" false positive
	// ---------------------------------------------------------------------
	//
	// In regions with essentially no baryons (nB ~ 0) or no nucleons (nn~0 or np~0),
	// DU is not physically meaningful, and the triangle inequality degenerates.
	//
	// Choose a very small nB_min relative to any physical DU onset (~0.31 fm^-3)
	// so we only exclude vacuum/crust-like zones.
	constexpr double nB_min = 1.0e-6; // fm^-3 (numerical/semantic guard, not a DU threshold)
	constexpr double n_min = 1.0e-12; // fm^-3 (ensures nn and np are genuinely present)

	// ---------------------------------------------------------------------
	// Define the DU boundary robustly for your indexing convention:
	// - index 0 is the core
	// - scan outward and identify the last index in the *contiguous* DU-allowed region
	// ---------------------------------------------------------------------
	long last_allowed = -1;
	double boundary_r_km = 0.0;
	bool seen_allowed_core_region = false;

	for (std::size_t i = 0; i < n; ++i)
	{
		// Read background densities and fractions.
		const double nB = (*m_nb)[i];
		const double Yn = (*Yn_col)[i];
		const double Yp = (*Yp_col)[i];
		const double Ye = (*Ye_col)[i];

		// Basic validity checks (fractions should be non-negative; nB must be finite).
		// If this fails, mark not allowed and continue.
		if (!std::isfinite(nB) || !std::isfinite(Yn) || !std::isfinite(Yp) || !std::isfinite(Ye) ||
			Yn < 0.0 || Yp < 0.0 || Ye < 0.0)
		{
			mask[i] = 0;
			continue;
		}

		// Exclude crust/vacuum-like zones with negligible baryons.
		if (nB < nB_min)
		{
			mask[i] = 0;
			// Do not "break" here: the profile might (in principle) have non-monotonic nB,
			// but the contiguous-core boundary logic below will still protect last_allowed.
			continue;
		}

		// Convert fractions to number densities in fm^-3.
		const double nn = Yn * nB;
		const double np = Yp * nB;
		const double ne = Ye * nB;

		// Require nucleons present. This is the critical guard that prevents
		// 0 <= kFe from making DU appear allowed in an electron-only crust.
		if (nn < n_min || np < n_min)
		{
			mask[i] = 0;
			continue;
		}

		// (Optional) require leptons; typically ne>0 in npe matter, but keep it safe.
		if (ne < n_min)
		{
			mask[i] = 0;
			continue;
		}

		// Compute Fermi momenta (fm^-1). cbrt is safe for positive inputs.
		const double kFn = std::cbrt(three_pi2 * nn);
		const double kFp = std::cbrt(three_pi2 * np);
		const double kFe = std::cbrt(three_pi2 * ne);

		// Kinematic allowance.
		const bool allowed = (kFn <= (kFp + kFe));
		mask[i] = allowed ? 1u : 0u;

		// Boundary selection: end of the contiguous DU-allowed core region.
		if (allowed)
		{
			seen_allowed_core_region = true;
			last_allowed = static_cast<long>(i);
			boundary_r_km = (*m_r)[i];
		}
		else if (seen_allowed_core_region)
		{
			// Once we have left the DU-allowed region while scanning outward from the core,
			// we stop: DU should not reappear at larger radius in a physically consistent
			// monotone profile, and this avoids "islands" caused by numerical noise.
			break;
		}
	}

	m_durca_mask = std::make_unique<std::vector<std::uint8_t>>(std::move(mask));
	m_durca_last_allowed = last_allowed;
	m_durca_boundary_r_km = (last_allowed >= 0) ? boundary_r_km : 0.0;
}

//--------------------------------------------------------------
// Direct Urca kinematic allowance mask (cached).
// mask[i] = 1 if direct Urca is kinematically allowed at radius index i,
// else 0. Computed from Fermi-momentum triangle condition:
//   kFn <= kFp + kFe  with kF = (3*pi^2*n)^(1/3).
// Number densities must be provided in fm^-3 for n_n, n_p, n_e.
// Cache invalidation is based on StarProfile versioning.
const std::vector<std::uint8_t> *StarContext::DirectUrcaMask() const
{
	if (!IsValid())
		return nullptr;

	RefreshDerivedCachesIfNeeded_();

	return m_durca_mask.get();
}

//--------------------------------------------------------------
// Last index (largest r index) where direct Urca is allowed.
// Returns -1 if mask unavailable or no region allows DU.
long StarContext::DirectUrcaLastAllowedIndex() const
{
	if (!IsValid())
		return -1;

	RefreshDerivedCachesIfNeeded_();

	return m_durca_last_allowed;
}

//--------------------------------------------------------------
// Radius (km) at the last allowed index, or 0 if not available.
double StarContext::DirectUrcaBoundaryRadius_km() const
{
	if (!IsValid())
		return 0.0;

	RefreshDerivedCachesIfNeeded_();

	return m_durca_boundary_r_km;
}

//--------------------------------------------------------------
const Zaki::Vector::DataColumn *StarContext::ChargeFractionYq() const
{
	RefreshDerivedCachesIfNeeded_();
	if (!m_Yq_cache)
		BuildYqCache_();
	return m_Yq_cache.get();
}

//--------------------------------------------------------------
void StarContext::BuildYqCache_() const
{
	if (!m_prof || !m_nb)
	{
		m_Yq_cache.reset();
		return;
	}

	const std::size_t N = Size();
	auto Yq = std::make_unique<Zaki::Vector::DataColumn>(N);

	// Prefer: if the profile already has a Yq species column, just use it.
	// Many CompOSE exports do not include it directly, so we compute it from strong-sector species fractions.
	//
	// We assume species columns are fractions Y_i = n_i/nB, consistent with your DirectUrca cache comment.
	//
	// Known CompOSE strong-sector numeric labels (as strings):
	const int codes[] = {10, 11, 20, 21, 22, 23, 100, 110, 111, 112, 120, 121, 500, 501, 502};

	struct Term
	{
		const Zaki::Vector::DataColumn *col;
		double q;
	};
	std::vector<Term> terms;
	terms.reserve(sizeof(codes) / sizeof(codes[0]));

	for (int code : codes)
	{
		const std::string label = std::to_string(code);
		const auto *col = m_prof->GetSpeciesPtr(label);
		if (!col)
			continue;

		double q = 0.0;
		if (!ComposeCharge(code, q))
			continue;
		terms.push_back({col, q});
	}

	if (terms.empty())
	{
		m_Yq_cache.reset();
		return;
	}

	for (std::size_t i = 0; i < N; ++i)
	{
		double yq = 0.0;
		for (const auto &t : terms)
			yq += t.q * (*(t.col))[i];

		(*Yq)[i] = yq;
	}

	m_Yq_cache = std::move(Yq);
}

//--------------------------------------------------------------
double StarContext::HeatCapacityStar_Tinf(double Tinf_MeV,
										  const CompactStar::EOS::CompOSE_Thermo &thermo) const
{
	// Ensure derived caches are in sync with profile version
	RefreshDerivedCachesIfNeeded_();

	// Build/rebuild if needed
	if (!m_cv_cache.loaded ||
		m_cv_cache.prof_version != ProfileVersion() ||
		m_cv_cache.thermo_tag != static_cast<const void *>(&thermo))
	{
		BuildHeatCapacityCache_(thermo);
	}

	const auto &Tg = m_cv_cache.Tinf_MeV;
	const auto &Cg = m_cv_cache.C_star;
	if (Tg.size() < 2)
		throw std::runtime_error("HeatCapacityStar_Tinf: cache invalid (grid too small).");

	// Clamp
	if (Tinf_MeV <= Tg.front())
		return Cg.front();
	if (Tinf_MeV >= Tg.back())
		return Cg.back();

	// Interpolate in log(T) for stability across decades
	const std::size_t i = Bracket(Tg, Tinf_MeV, m_cv_cache.last_i);
	m_cv_cache.last_i = i;

	const double T0 = Tg[i];
	const double T1 = Tg[i + 1];
	const double w = (std::log(Tinf_MeV) - std::log(T0)) / (std::log(T1) - std::log(T0));

	return Lerp(Cg[i], Cg[i + 1], w);
}

void StarContext::BuildHeatCapacityCache_(const CompactStar::EOS::CompOSE_Thermo &thermo) const
{
	if (!IsValid() || !m_nb || !m_nu)
	{
		m_cv_cache = HeatCapacityCache{};
		return;
	}

	// Build GR geometry weights from existing columns
	const GeometryCache geo(*this);

	const auto &r = geo.R();   // km
	const auto &wv = geo.WV(); // 4*pi*r^2*exp(Lambda)
	const auto &eminusnu = geo.ExpMinusNu();

	const std::size_t N = geo.Size();
	if (N < 2)
	{
		m_cv_cache = HeatCapacityCache{};
		return;
	}

	// Need Yq(r)
	const auto *Yq_col = ChargeFractionYq();
	if (!Yq_col)
		throw std::runtime_error("BuildHeatCapacityCache_: missing strong-sector composition to compute Yq(r).");

	// Temperature grid for T_infty (MeV). You can tune these later.
	const double Tinf_min = 1e-5; // MeV
	const double Tinf_max = 1.0;  // MeV
	const std::size_t NT = 160;

	HeatCapacityCache cache;
	cache.loaded = true;
	cache.prof_version = ProfileVersion();
	cache.thermo_tag = static_cast<const void *>(&thermo);
	cache.last_i = 0;

	cache.Tinf_MeV = LogSpace(Tinf_min, Tinf_max, NT);
	cache.C_star.assign(NT, 0.0);

	// C(Tinf) = ∫ cV(Tlocal, nb, Yq) * WV(r) dr, with Tlocal = Tinf*exp(-nu).
	for (std::size_t k = 0; k < NT; ++k)
	{
		const double Tinf = cache.Tinf_MeV[k];
		double sum = 0.0;

		for (std::size_t i = 0; i < N - 1; ++i)
		{
			const double dr = r[i + 1] - r[i]; // km

			const double nb0 = (*m_nb)[i];
			const double nb1 = (*m_nb)[i + 1];

			const double yq0 = (*Yq_col)[i];
			const double yq1 = (*Yq_col)[i + 1];

			const double T0 = Tinf * eminusnu[i];
			const double T1 = Tinf * eminusnu[i + 1];

			const double cv0 = thermo.CvDensity_ForCooling(T0, nb0, yq0);
			const double cv1 = thermo.CvDensity_ForCooling(T1, nb1, yq1);

			const double f0 = cv0 * wv[i];
			const double f1 = cv1 * wv[i + 1];

			sum += 0.5 * (f0 + f1) * dr;
		}

		cache.C_star[k] = sum;
	}

	m_cv_cache = std::move(cache);
}
//--------------------------------------------------------------

//==============================================================
} // namespace CompactStar::Physics::Evolution