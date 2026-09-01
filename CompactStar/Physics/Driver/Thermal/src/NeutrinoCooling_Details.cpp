// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025
 * Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file NeutrinoCooling_Details.cpp
 * @brief Out-of-line implementation of NeutrinoCooling shared computations.
 *
 * Rationale:
 *  - Keep heavy computations and structure/microphysics dereferences out of headers.
 *  - Avoid incomplete-type problems when using ctx-provided objects.
 *  - Centralize neutrino luminosity computation so RHS and diagnostics agree.
 */

#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling_Details.hpp"
#include "CompactStar/Units.hpp"
#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp"

#include <cmath>
#include <stdexcept>

// Diagnostics packet API
#include "CompactStar/Physics/Evolution/Diagnostics/DiagnosticPacket.hpp"

// If you have these in your project, include them here; otherwise remove.
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
// #include "CompactStar/Physics/Microphysics/NeutrinoEmissivity.hpp"

#include <Zaki/Util/Instrumentor.hpp>	  // PROFILE_FUNCTION
#include <Zaki/Util/Logger.hpp>			  // Z_LOG_INFO/WARNING/ERROR
#include <Zaki/Vector/IntegrateTrapz.hpp> // for volume integrals

#include <Zaki/Physics/Constants.hpp> // unit conversions
namespace CompactStar::Physics::Driver::Thermal
{

// -----------------------------------------------------------------------------
namespace Detail
{

namespace // unnamed => internal linkage (preferred over `static`)
{
// -----------------------------------------------------------------------------
//  Build Neutrino-Cooling Cache
// -----------------------------------------------------------------------------
static void BuildNeutrinoCoolingCache(const CompactStar::Physics::Evolution::StarContext &sc,
									  const CompactStar::Physics::Evolution::DriverContext &ctx,
									  CompactStar::Physics::Driver::Thermal::NeutrinoCoolingCachePayload &out)
{
	// Reset output (overwrite semantics).
	out = {};
	out.valid = false;

	// Require geometry cache (metric + radial grid + pre-weighted volume factor).
	if (!ctx.geo)
		return;

	const auto &R_km = ctx.geo->R();
	const auto &expMinusNu = ctx.geo->ExpMinusNu();
	const auto &wVExp2Nu = ctx.geo->WVExp2Nu(); // 4π r^2 e^Λ e^{2ν}

	const std::size_t N = R_km.Size();
	if (N < 2)
		return;

	if (expMinusNu.Size() != N || wVExp2Nu.Size() != N)
		return;

	// Require mass density column (already cached by StarContext).
	const auto *rho_gcm3_col = sc.MassDensity_gcm3();
	if (!rho_gcm3_col || rho_gcm3_col->Size() != N)
		return;

	// DU boundary index (last allowed). -1 means nowhere allowed/unavailable.
	long durca_last = sc.DirectUrcaLastAllowedIndex();
	if (durca_last >= static_cast<long>(N))
		durca_last = static_cast<long>(N) - 1;

	out.n_zones = N;
	out.durca_last_allowed = durca_last;

	// ---------------------------------------------------------------------
	// Coefficient conventions
	// ---------------------------------------------------------------------
	// The km^3 -> cm^3 conversion applied after integrating the geometry weights
	// is CompactStar::Units::KM3_TO_CM3. NOTE: this is done once per profile
	// version, not in the per-step hot path.

	// (T/1e9)^n factors moved into coefficients => multiply by (1e-9)^n
	constexpr double INV_1E9_POW6 = 1.0e-54; // (1e-9)^6
	constexpr double INV_1E9_POW8 = 1.0e-72; // (1e-9)^8

	// Placeholder normalizations (must match your emissivity model).
	constexpr double Q0_DU = 1.0e27; // erg cm^-3 s^-1 * rho15 * T9^6
	constexpr double Q0_MU = 1.0e21; // erg cm^-3 s^-1 * rho15 * T9^8

	// Helper: e^{-nu}^6 and e^{-nu}^8 using expMinusNu (already e^{-nu}).
	auto pow6 = [](double a) noexcept
	{
		const double a2 = a * a;
		const double a4 = a2 * a2;
		return a4 * a2;
	};
	auto pow8 = [](double a) noexcept
	{
		const double a2 = a * a;
		const double a4 = a2 * a2;
		return a4 * a4;
	};

	// ---------------------------------------------------------------------
	// Trapezoid integrators for the *coefficients* (temperature factored out)
	// ---------------------------------------------------------------------
	auto integrate_coeff = [&](std::size_t i0, std::size_t i1, int n_power) -> double
	{
		// Integrate over indices [i0..i1] inclusive (requires i1>=i0+1).
		if (i1 <= i0 || i1 >= N)
			return 0.0;

		auto f_at = [&](std::size_t i) -> double
		{
			const double rho_gcm3 = (*rho_gcm3_col)[i];
			const double w = wVExp2Nu[i];
			const double emn = expMinusNu[i];

			if (!(rho_gcm3 > 0.0) || !std::isfinite(rho_gcm3))
				return 0.0;
			if (!(w > 0.0) || !std::isfinite(w))
				return 0.0;
			if (!(emn > 0.0) || !std::isfinite(emn))
				return 0.0;

			// rho15 = rho / 1e15
			const double rho15 = rho_gcm3 * 1.0e-15;

			// Need e^{-n nu} = (e^{-nu})^n
			double emn_n = 0.0;
			if (n_power == 6)
				emn_n = pow6(emn);
			else if (n_power == 8)
				emn_n = pow8(emn);
			else
				return 0.0;

			// Coefficient integrand (still needs dr integration):
			// rho15 * e^{-n nu} * wVExp2Nu
			return rho15 * emn_n * w;
		};

		double f0 = f_at(i0);
		double I_km3 = 0.0;

		for (std::size_t i = i0; i < i1; ++i)
		{
			const double dr = R_km[i + 1] - R_km[i];
			if (!(dr > 0.0) || !std::isfinite(dr))
			{
				// Grid should be monotone; if not, skip safely.
				f0 = f_at(i + 1);
				continue;
			}

			const double f1 = f_at(i + 1);
			I_km3 += 0.5 * (f0 + f1) * dr;
			f0 = f1;
		}

		return I_km3;
	};

	// ---------------------------------------------------------------------
	// Build MU coefficient over full star [0..N-1]
	// ---------------------------------------------------------------------
	const double I_MU_km3 = integrate_coeff(0, N - 1, 8);

	// K_MU: [erg/s/K^8] = Q0_MU * (1e-9)^8 * (km^3->cm^3) * I_MU_km3
	out.K_MU_erg_s_K8 = Q0_MU * INV_1E9_POW8 * CompactStar::Units::KM3_TO_CM3 * I_MU_km3;

	// ---------------------------------------------------------------------
	// Build DU coefficient over allowed region only [0..durca_last]
	// ---------------------------------------------------------------------
	double I_DU_km3 = 0.0;
	if (durca_last >= 1) // need at least two points to integrate
	{
		const std::size_t M = static_cast<std::size_t>(durca_last);
		I_DU_km3 = integrate_coeff(0, M, 6);
	}

	out.K_DU_erg_s_K6 = Q0_DU * INV_1E9_POW6 * CompactStar::Units::KM3_TO_CM3 * I_DU_km3;

	// PBF remains a hook (0 unless implemented)
	out.K_PBF = 0.0;

	// Basic sanity: allow zero (e.g., no DU region), but require finiteness.
	if (!std::isfinite(out.K_MU_erg_s_K8))
		out.K_MU_erg_s_K8 = 0.0;
	if (!std::isfinite(out.K_DU_erg_s_K6))
		out.K_DU_erg_s_K6 = 0.0;

	out.valid = true;

	// Record built version.
	out.built_version = sc.ProfileVersion();
}
// -----------------------------------------------------------------------------
// k_B in MeV/K
constexpr double MEV_PER_K = 8.617333262145e-11;
} // namespace
// -----------------------------------------------------------------------------
// namespace
// {

// // -----------------------------------------------------------------------------
// // Unit conversions
// // -----------------------------------------------------------------------------
// constexpr double KM3_TO_CM3 = 1.0e15; // (1e5 cm)^3

// // -----------------------------------------------------------------------------
// // Convert eps [km^-2] -> rho [g/cm^3] via MeV fm^-3 bridge
// // -----------------------------------------------------------------------------
// // inline double RhoFromEps_km2(double eps_km2) noexcept
// // {
// // 	// eps(MeV fm^-3) = eps(km^-2) / (MeV fm^-3 -> km^-2)
// // 	const double eps_MeV_fm3 = eps_km2 / Zaki::Physics::MEV_FM3_2_INV_KM2;
// // 	// rho(g/cm^3) = eps(MeV fm^-3) * (MeV fm^-3 -> g/cm^3)
// // 	return eps_MeV_fm3 * Zaki::Physics::MEV_FM3_2_G_CM3;
// // }
// // -----------------------------------------------------------------------------
// // /**
// //  * @brief Minimal placeholder emissivity model for wiring tests.
// //  *
// //  * This exists only so the base structure (catalog/packet/timeseries) can run
// //  * before you plug in real emissivities.
// //  *
// //  * Replace with: DUrca/MUrca/PBF emissivities integrated over proper volume.
// //  *
// //  * @param Tinf_K Redshifted interior temperature [K].
// //  * @return A tiny luminosity scale [erg/s] for non-zero RHS during smoke tests.
// //  */
// // double PlaceholderLuminosity(double Tinf_K)
// // {
// // 	// Deliberately small-ish for stability; choose any smooth T dependence.
// // 	// Typical neutrino luminosities can be huge; this is not physical.
// // 	// L ~ 1e30 * (T/1e8)^6 erg/s as a harmless placeholder.
// // 	const double T8 = Tinf_K / 1.0e8;
// // 	return 1.0e30 * std::pow(T8, 6.0);
// // }

// // -----------------------------------------------------------------------------
// // Minimal baseline emissivities (replace coefficients with your chosen fits)
// // Returns Q in erg cm^-3 s^-1.
// // -----------------------------------------------------------------------------
// double Q_ModifiedUrca_erg_cm3_s(double T_local_K, double rho_g_cm3) noexcept
// {
// 	// MUrca ~ T^8; include mild density dependence so profiles matter.
// 	const double T9 = T_local_K / 1.0e9;
// 	const double rho15 = rho_g_cm3 / 1.0e15;

// 	// Replace this coefficient with your preferred Yakovlev/Pethick fit.
// 	// Typical order for MUrca is ~ 1e21 * (rho/1e15) * T9^8 erg cm^-3 s^-1.
// 	const double Q0 = 1.0e21;

// 	// Fast T^8 via squaring
// 	const double T9_2 = T9 * T9;
// 	const double T9_4 = T9_2 * T9_2;
// 	const double T9_8 = T9_4 * T9_4;

// 	// return Q0 * rho15 * std::pow(T9, 8.0);
// 	return Q0 * rho15 * T9_8;
// }
// // -----------------------------------------------------------------------------
// double Q_DirectUrca_erg_cm3_s(double T_local_K, double rho_g_cm3) noexcept
// {
// 	// DUrca ~ T^6; much stronger normalization.
// 	const double T9 = T_local_K / 1.0e9;
// 	const double rho15 = rho_g_cm3 / 1.0e15;

// 	// Replace with your DUrca fit + threshold logic (proton fraction etc).
// 	const double Q0 = 1.0e27;

// 	// Fast T^6 = T^2 * T^4
// 	const double T9_2 = T9 * T9;
// 	const double T9_4 = T9_2 * T9_2;
// 	const double T9_6 = T9_4 * T9_2;

// 	// return Q0 * rho15 * std::pow(T9, 6.0);
// 	return Q0 * rho15 * T9_6;
// }
// // -----------------------------------------------------------------------------
// struct TPowers9
// {
// 	// Base powers (always cheap to compute once)
// 	double T9;	 // optional, but often useful for debugging/logging
// 	double T9_2; // T9^2
// 	double T9_4; // T9^4

// 	// Optional cached higher powers
// 	mutable double T9_6 = 0.0; // computed on demand
// 	mutable double T9_8 = 0.0; // computed on demand
// 	mutable bool has_T9_6 = false;
// 	mutable bool has_T9_8 = false;

// 	// Compute T9^6 = T9^4 * T9^2 (on demand)
// 	inline double Pow6() const noexcept
// 	{
// 		if (!has_T9_6)
// 		{
// 			T9_6 = T9_4 * T9_2;
// 			has_T9_6 = true;
// 		}
// 		return T9_6;
// 	}

// 	// Compute T9^8 = T9^4 * T9^4 (on demand)
// 	inline double Pow8() const noexcept
// 	{
// 		if (!has_T9_8)
// 		{
// 			T9_8 = T9_4 * T9_4;
// 			has_T9_8 = true;
// 		}
// 		return T9_8;
// 	}
// };

// // -----------------------------------------------------------------------------
// inline TPowers9 MakeTPowers9(double T_local_K) noexcept
// {
// 	const double T9 = T_local_K * 1.0e-9;
// 	const double T9_2 = T9 * T9;
// 	const double T9_4 = T9_2 * T9_2;

// 	TPowers9 p;
// 	p.T9 = T9;
// 	p.T9_2 = T9_2;
// 	p.T9_4 = T9_4;
// 	return p;
// }

// // -----------------------------------------------------------------------------
// inline double Q_ModifiedUrca_from_T9powers(double rho_g_cm3, const TPowers9 &p) noexcept
// {
// 	constexpr double Q0 = 1.0e21;
// 	const double rho15 = rho_g_cm3 * 1.0e-15;
// 	return Q0 * rho15 * p.Pow8();
// }

// // -----------------------------------------------------------------------------
// inline double Q_DirectUrca_from_T9powers(double rho_g_cm3, const TPowers9 &p) noexcept
// {
// 	constexpr double Q0 = 1.0e27;
// 	const double rho15 = rho_g_cm3 * 1.0e-15;
// 	return Q0 * rho15 * p.Pow6();
// }

// // -----------------------------------------------------------------------------
// } // namespace

// -----------------------------------------------------------------------------
//  ComputeDerived
// -----------------------------------------------------------------------------
// NeutrinoCooling_Details ComputeDerived(const NeutrinoCooling &drv,
// 									   const Evolution::StateVector &Y,
// 									   const Evolution::DriverContext &ctx)
// {
// 	PROFILE_FUNCTION();

// 	NeutrinoCooling_Details d;

// 	// ---------------------------------------------------------------------
// 	// 1) Extract evolved DOF: T_inf
// 	// ---------------------------------------------------------------------
// 	const auto &thermal = Y.GetThermal();
// 	if (thermal.NumComponents() == 0)
// 	{
// 		d.ok = false;
// 		d.message = "ThermalState has zero components.";
// 		return d;
// 	}

// 	d.Tinf_K = thermal.Tinf();
// 	if (!(d.Tinf_K > 0.0))
// 	{
// 		d.ok = false;
// 		d.message = "Tinf <= 0; neutrino cooling ill-defined.";
// 		return d;
// 	}

// 	// ---------------------------------------------------------------------
// 	// 2) Resolve/validate options used in RHS
// 	// ---------------------------------------------------------------------
// 	// Your NeutrinoCooling::Options currently includes global_scale and channel toggles.
// 	// You also need an effective heat capacity; choose one policy:
// 	//   A) store C_eff inside NeutrinoCooling::Options (recommended for symmetry w/ PhotonCooling),
// 	//   B) fetch from ctx (if you already have a thermal capacity model),
// 	//   C) reuse the same C_eff as photon cooling in a shared thermal config object.
// 	//
// 	// For now, we assume NeutrinoCooling has a C_eff accessor or that you will add it.
// 	//
// 	// If you don't have it yet, set a temporary constant here and move it later.
// 	//
// 	// IMPORTANT: make this match your driver header/cpp choices.
// 	//
// 	// 2) Heat capacity policy (keep your requested constant for now)
// 	d.C_eff_erg_K = 1.0e40;
// 	if (!(d.C_eff_erg_K > 0.0))
// 	{
// 		d.ok = false;
// 		d.message = "C_eff <= 0.";
// 		return d;
// 	}

// 	if (!(drv.GetOptions().global_scale > 0.0))
// 	{
// 		// Disabled but valid: no cooling contribution.
// 		d.ok = true;
// 		d.message = "cooling disabled: global_scale <= 0.";
// 		d.L_nu_inf_erg_s = 0.0;
// 		d.dTinf_dt_K_s = 0.0;
// 		d.dLnTinf_dt_1_s = 0.0;
// 		return d;
// 	}

// 	// If all channels disabled, treat as disabled-but-valid.
// 	if (!drv.GetOptions().include_direct_urca &&
// 		!drv.GetOptions().include_modified_urca &&
// 		!drv.GetOptions().include_pair_breaking)
// 	{
// 		d.ok = true;
// 		d.message = "cooling disabled: all neutrino channels disabled by options.";
// 		d.L_nu_inf_erg_s = 0.0;
// 		d.dTinf_dt_K_s = 0.0;
// 		d.dLnTinf_dt_1_s = 0.0;
// 		return d;
// 	}

// 	// ---------------------------------------------------------------------
// 	// 3) Compute neutrino luminosity at infinity (placeholder vs real)
// 	// ---------------------------------------------------------------------
// 	//
// 	// This is where the real work will go:
// 	//   - loop over radial zones
// 	//   - compute local T(r) (isothermal redshifted => T_local(r) = Tinf * e^{-nu(r)})
// 	//   - compute emissivities Q_nu (erg cm^-3 s^-1) per process via microphysics
// 	//   - multiply by proper volume element dV_proper (cm^3)
// 	//   - apply redshift factors to get L_inf
// 	//
// 	if (!ctx.geo)
// 	{
// 		d.ok = false;
// 		d.message = "ctx.geo == nullptr";
// 		return d;
// 	}

// 	const auto &R_km = ctx.geo->R();
// 	const auto &e_minus_nu = ctx.geo->ExpMinusNu();
// 	const auto &wV_e2nu = ctx.geo->WVExp2Nu(); // 4π r^2 e^Λ e^{2ν}

// 	const std::size_t N = R_km.Size();
// 	if (N < 2)
// 	{
// 		d.ok = true;
// 		d.n_zones = N;
// 		d.L_nu_inf_erg_s = 0.0;
// 		d.dTinf_dt_K_s = 0.0;
// 		d.dLnTinf_dt_1_s = 0.0;
// 		return d;
// 	}
// 	if (e_minus_nu.Size() != N || wV_e2nu.Size() != N)
// 	{
// 		d.ok = false;
// 		d.message = "GeometryCache size mismatch (R vs ExpMinusNu vs WVExp2Nu).";
// 		return d;
// 	}

// 	// 4) Require structure for eps(r)
// 	if (!ctx.star)
// 	{
// 		d.ok = false;
// 		d.message = "ctx.star == nullptr (StarContext required for eps(r)).";
// 		return d;
// 	}

// 	const auto *rho_g_cm3_col = ctx.star->MassDensity_gcm3(); // rho [g/cm^3]
// 	if (!rho_g_cm3_col || rho_g_cm3_col->Size() != N)
// 	{
// 		d.ok = false;
// 		d.message = "StarContext::MassDensity_gcm3 missing or size mismatch with geometry grid.";
// 		return d;
// 	}

// 	d.ok = true;
// 	d.has_structure = true;
// 	d.n_zones = N;

// 	// 5) Build integrands in (erg cm^-3 s^-1) * (km^2)?? no:
// 	// integrand[i] = Q(erg/cm^3/s) * WVExp2Nu(km^2?) * ... then dr(km) => km^3; convert km^3->cm^3 later
// 	// ---------------------------------------------------------------------
// 	// 5) Single-pass trapezoid integration over enabled channels
// 	// ---------------------------------------------------------------------

// 	// Hoist options locally (cuts repeated virtual/inline calls and branches).
// 	const auto &opt = drv.GetOptions();
// 	const bool do_DU = opt.include_direct_urca;
// 	const bool do_MU = opt.include_modified_urca;
// 	const bool do_PBF = opt.include_pair_breaking;
// 	const double gscale = opt.global_scale;

// 	// Fast exit if globally disabled (you already handled <=0 earlier, but keep it safe).
// 	if (!(gscale > 0.0))
// 	{
// 		d.L_nu_DU_inf_erg_s = 0.0;
// 		d.L_nu_MU_inf_erg_s = 0.0;
// 		d.L_nu_PBF_inf_erg_s = 0.0;
// 		d.L_nu_inf_erg_s = 0.0;
// 		d.dTinf_dt_K_s = 0.0;
// 		d.dLnTinf_dt_1_s = 0.0;
// 		return d;
// 	}

// 	// -------------------------------------------------
// 	// DUrca boundary (cached in StarContext). If <0, DUrca is nowhere allowed.
// 	long durca_last = -1;
// 	if (do_DU)
// 	{
// 		durca_last = ctx.star->DirectUrcaLastAllowedIndex();
// 		if (durca_last >= static_cast<long>(N))
// 			durca_last = static_cast<long>(N) - 1;
// 	}
// 	// -------------------------------------------------

// 	// Helper lambda: compute emissivity-weighted integrand at zone i.
// 	// Returns f_total(i) = (Q_DU + Q_MU + Q_PBF) * wV_e2nu[i]
// 	// and also provides optional per-channel f's for bookkeeping.
// 	// auto eval_f = [&](std::size_t i, double &f_DU, double &f_MU, double &f_PBF) -> double
// 	// {
// 	// 	f_DU = 0.0;
// 	// 	f_MU = 0.0;
// 	// 	f_PBF = 0.0;

// 	// 	// Local temperature (isothermal interior w/ redshift)
// 	// 	const double T_local = d.Tinf_K * e_minus_nu[i];
// 	// 	if (!(T_local > 0.0) || !std::isfinite(T_local))
// 	// 		return 0.0;

// 	// 	// const double rho_g_cm3 = RhoFromEps_km2((*eps_km2_col)[i]);
// 	// 	const double rho_g_cm3 = (*rho_g_cm3_col)[i];
// 	// 	if (!(rho_g_cm3 > 0.0) || !std::isfinite(rho_g_cm3))
// 	// 		return 0.0;

// 	// 	const auto p = MakeTPowers9(T_local);

// 	// 	// Compute only enabled channels.
// 	// 	// Multiply by geometric weight once at the end.
// 	// 	double Qsum = 0.0;

// 	// 	if (do_DU)
// 	// 	{
// 	// 		const double Q = Q_DirectUrca_from_T9powers(rho_g_cm3, p);
// 	// 		if (std::isfinite(Q) && Q > 0.0)
// 	// 		{
// 	// 			f_DU = Q;
// 	// 			Qsum += Q;
// 	// 		}
// 	// 	}

// 	// 	if (do_MU)
// 	// 	{
// 	// 		const double Q = Q_ModifiedUrca_from_T9powers(rho_g_cm3, p);
// 	// 		if (std::isfinite(Q) && Q > 0.0)
// 	// 		{
// 	// 			f_MU = Q;
// 	// 			Qsum += Q;
// 	// 		}
// 	// 	}

// 	// 	if (do_PBF)
// 	// 	{
// 	// 		// Placeholder until implemented:
// 	// 		// const double Q = Q_PBF_erg_cm3_s(T_local, rho_g_cm3, ...);
// 	// 		// if (std::isfinite(Q) && Q > 0.0) { f_PBF = Q; Qsum += Q; }
// 	// 		// For now, leave zero.
// 	// 	}

// 	// 	const double w = wV_e2nu[i];
// 	// 	if (!(w > 0.0) || !std::isfinite(w))
// 	// 		return 0.0;

// 	// 	// Convert per-channel to fully-weighted integrands for trapezoid accumulation
// 	// 	// (units: Q * w, integrated over dr in km => "km^3" factor to be converted later).
// 	// 	f_DU *= w;
// 	// 	f_MU *= w;
// 	// 	f_PBF *= w;

// 	// 	return Qsum * w;
// 	// };

// 	// -------------------------------------------------
// 	auto eval_f_MU_PBF = [&](std::size_t i, double &f_MU, double &f_PBF) -> double
// 	{
// 		f_MU = 0.0;
// 		f_PBF = 0.0;

// 		const double T_local = d.Tinf_K * e_minus_nu[i];
// 		if (!(T_local > 0.0) || !std::isfinite(T_local))
// 			return 0.0;

// 		const double rho_g_cm3 = (*rho_g_cm3_col)[i];
// 		if (!(rho_g_cm3 > 0.0) || !std::isfinite(rho_g_cm3))
// 			return 0.0;

// 		const double w = wV_e2nu[i];
// 		if (!(w > 0.0) || !std::isfinite(w))
// 			return 0.0;

// 		const auto p = MakeTPowers9(T_local);

// 		double Qsum = 0.0;

// 		if (do_MU)
// 		{
// 			const double Q = Q_ModifiedUrca_from_T9powers(rho_g_cm3, p);
// 			if (std::isfinite(Q) && Q > 0.0)
// 			{
// 				f_MU = Q * w;
// 				Qsum += Q;
// 			}
// 		}

// 		if (do_PBF)
// 		{
// 			// placeholder; same pattern:
// 			// const double Q = Q_PBF_from_...(rho_g_cm3, p, ...);
// 			// if (finite && >0) { f_PBF = Q*w; Qsum += Q; }
// 		}

// 		return Qsum * w; // weighted total (MU+PBF)
// 	};
// 	// -------------------------------------------------

// 	// // Trapezoid accumulators (km^3 in your current convention).
// 	// double I_total_km3 = 0.0;
// 	// double I_DU_km3 = 0.0;
// 	// double I_MU_km3 = 0.0;
// 	// double I_PBF_km3 = 0.0;

// 	// // Evaluate at i=0 to seed trapezoid.
// 	// double f0_DU, f0_MU, f0_PBF;
// 	// double f0 = eval_f(0, f0_DU, f0_MU, f0_PBF);

// 	// for (std::size_t i = 0; i + 1 < N; ++i)
// 	// {
// 	// 	const double r0 = R_km[i];
// 	// 	const double r1 = R_km[i + 1];
// 	// 	const double dr = r1 - r0;

// 	// 	// // Skip non-monotonic or degenerate steps.
// 	// 	// // (If your grid is guaranteed monotonic, you can remove this branch.)
// 	// 	// if (!(dr > 0.0) || !std::isfinite(dr))
// 	// 	// {
// 	// 	// 	// Re-seed next point in case of weird grid.
// 	// 	// 	double tmpDU, tmpMU, tmpPBF;
// 	// 	// 	f0 = eval_f(i + 1, tmpDU, tmpMU, tmpPBF);
// 	// 	// 	f0_DU = tmpDU;
// 	// 	// 	f0_MU = tmpMU;
// 	// 	// 	f0_PBF = tmpPBF;
// 	// 	// 	continue;
// 	// 	// }

// 	// 	double f1_DU, f1_MU, f1_PBF;
// 	// 	const double f1 = eval_f(i + 1, f1_DU, f1_MU, f1_PBF);

// 	// 	// Standard trapezoid: 0.5*(f0+f1)*dr
// 	// 	const double wtrap = 0.5 * dr;

// 	// 	I_total_km3 += (f0 + f1) * wtrap;

// 	// 	// If you want per-channel bookkeeping, accumulate them too.
// 	// 	// This costs only a few adds/mults and avoids extra integrations.
// 	// 	if (do_DU)
// 	// 		I_DU_km3 += (f0_DU + f1_DU) * wtrap;
// 	// 	if (do_MU)
// 	// 		I_MU_km3 += (f0_MU + f1_MU) * wtrap;
// 	// 	if (do_PBF)
// 	// 		I_PBF_km3 += (f0_PBF + f1_PBF) * wtrap;

// 	// 	// Advance
// 	// 	f0 = f1;
// 	// 	f0_DU = f1_DU;
// 	// 	f0_MU = f1_MU;
// 	// 	f0_PBF = f1_PBF;
// 	// }

// 	// -------------------------------------------------
// 	double I_MU_PBF_km3 = 0.0;
// 	double I_MU_km3 = 0.0;
// 	double I_PBF_km3 = 0.0;

// 	double f0_MU, f0_PBF;
// 	double f0 = eval_f_MU_PBF(0, f0_MU, f0_PBF);

// 	for (std::size_t i = 0; i + 1 < N; ++i)
// 	{
// 		const double dr = R_km[i + 1] - R_km[i];
// 		const double wtrap = 0.5 * dr;

// 		double f1_MU, f1_PBF;
// 		const double f1 = eval_f_MU_PBF(i + 1, f1_MU, f1_PBF);

// 		I_MU_PBF_km3 += (f0 + f1) * wtrap;
// 		if (do_MU)
// 			I_MU_km3 += (f0_MU + f1_MU) * wtrap;
// 		if (do_PBF)
// 			I_PBF_km3 += (f0_PBF + f1_PBF) * wtrap;

// 		f0 = f1;
// 		f0_MU = f1_MU;
// 		f0_PBF = f1_PBF;
// 	}
// 	// -------------------------------------------------

// 	// -------------------------------------------------
// 	auto eval_f_DU = [&](std::size_t i) -> double
// 	{
// 		const double T_local = d.Tinf_K * e_minus_nu[i];
// 		if (!(T_local > 0.0) || !std::isfinite(T_local))
// 			return 0.0;

// 		const double rho_g_cm3 = (*rho_g_cm3_col)[i];
// 		if (!(rho_g_cm3 > 0.0) || !std::isfinite(rho_g_cm3))
// 			return 0.0;

// 		const double w = wV_e2nu[i];
// 		if (!(w > 0.0) || !std::isfinite(w))
// 			return 0.0;

// 		const auto p = MakeTPowers9(T_local);

// 		const double Q = Q_DirectUrca_from_T9powers(rho_g_cm3, p);
// 		if (!std::isfinite(Q) || Q <= 0.0)
// 			return 0.0;

// 		return Q * w;
// 	};

// 	double I_DU_km3 = 0.0;

// 	if (do_DU && durca_last >= 0)
// 	{
// 		const std::size_t M = static_cast<std::size_t>(durca_last);

// 		double f0_du = eval_f_DU(0);

// 		// integrate over [0..M], so loop i=0..M-1
// 		for (std::size_t i = 0; i < M; ++i)
// 		{
// 			const double dr = R_km[i + 1] - R_km[i];
// 			const double wtrap = 0.5 * dr;

// 			const double f1_du = eval_f_DU(i + 1);

// 			I_DU_km3 += (f0_du + f1_du) * wtrap;

// 			f0_du = f1_du;
// 		}
// 	}
// 	// -------------------------------------------------

// 	// -------------------------------------------------
// 	// Convert km^3 -> cm^3, apply global scale
// 	d.L_nu_DU_inf_erg_s = gscale * (I_DU_km3 * KM3_TO_CM3);
// 	d.L_nu_MU_inf_erg_s = gscale * (I_MU_km3 * KM3_TO_CM3);
// 	d.L_nu_PBF_inf_erg_s = gscale * (I_PBF_km3 * KM3_TO_CM3);

// 	d.L_nu_inf_erg_s = gscale * ((I_MU_PBF_km3 + I_DU_km3) * KM3_TO_CM3);

// 	// -------------------------------------------------
// 	// 6) Cooling rate
// 	d.dTinf_dt_K_s = -d.L_nu_inf_erg_s / d.C_eff_erg_K;
// 	d.dLnTinf_dt_1_s = d.dTinf_dt_K_s / d.Tinf_K;
// 	// -------------------------------------------------

// 	return d;
// }

// ------------------------------------------------------
// Helper for safe pow with small integer exponents.
static inline double PowInt(double x, int p) noexcept
{
	// p expected in {6,8} here
	if (p == 6)
	{
		const double x2 = x * x;
		const double x4 = x2 * x2;
		return x4 * x2;
	}
	if (p == 8)
	{
		const double x2 = x * x;
		const double x4 = x2 * x2;
		return x4 * x4;
	}
	return std::pow(x, static_cast<double>(p));
}

// -----------------------------------------------------------------------------
// ComputeDerived (cache-driven)
// -----------------------------------------------------------------------------
NeutrinoCooling_Details NeutrinoCooling_Details::ComputeDerived(const NeutrinoCooling &drv,
																const Evolution::StateVector &Y,
																const Evolution::DriverContext &ctx)
{
	NeutrinoCooling_Details d;

	// ------------------------------------------------------------
	// 1) Read evolved temperature DOF
	// ------------------------------------------------------------
	const auto &thermal = Y.GetThermal();
	if (thermal.NumComponents() == 0)
	{
		d.ok = false;
		d.message = "ThermalState has zero components.";
		return d;
	}

	d.Tinf_K = thermal.Tinf();
	if (!(d.Tinf_K > 0.0) || !std::isfinite(d.Tinf_K))
	{
		d.ok = false;
		d.message = "Tinf <= 0 or non-finite.";
		return d;
	}

	// ------------------------------------------------------------
	// 2) Resolve options / fast disable
	// ------------------------------------------------------------
	const auto &opt = drv.GetOptions();
	const double gscale = opt.global_scale;

	// Still return ok=true with zeros when disabled.
	if (!(gscale > 0.0) ||
		(!opt.include_direct_urca && !opt.include_modified_urca && !opt.include_pair_breaking))
	{
		d.ok = true;
		d.message = "cooling disabled by options/global_scale.";
		d.L_nu_inf_erg_s = 0.0;
		d.L_nu_DU_inf_erg_s = 0.0;
		d.L_nu_MU_inf_erg_s = 0.0;
		d.L_nu_PBF_inf_erg_s = 0.0;
		d.dTinf_dt_K_s = 0.0;
		d.dLnTinf_dt_1_s = 0.0;
		return d;
	}

	// ------------------------------------------------------------
	// 3) Heat capacity policy 
	// ------------------------------------------------------------
	// Keep whatever policy you currently use; the cache only accelerates L_nu.
	// d.C_eff_erg_K = 1.0e40;

	// ------------------------------------------------------------
	// 3) Heat capacity (GR-integrated, CompOSE-based)
	// ------------------------------------------------------------
	if (!ctx.thermo)
	{
		Z_LOG_ERROR("NeutrinoCooling requires ctx.thermo (CompOSE_Thermo) but it is nullptr.");
		d.ok = false;
		d.message = "ctx.thermo == nullptr (CompOSE_Thermo required for heat capacity).";
		return d;
	}

	// Convert evolved Tinf(K) to Tinf(MeV)
	const double Tinf_MeV = d.Tinf_K * MEV_PER_K;
	if (!(Tinf_MeV > 0.0) || !std::isfinite(Tinf_MeV))
	{
		d.ok = false;
		d.message = "Tinf_MeV <= 0 or non-finite after conversion.";
		return d;
	}

	// Require the star context BEFORE dereferencing it. The guard used to sit ~12 lines
	// below this call, where it could never fire (Phase 2A engineering correction).
	if (!ctx.star)
	{
		d.ok = false;
		d.message = "ctx.star == nullptr (StarContext required).";
		return d;
	}

	// ADR-0003 §13: never compute a luminosity from a mismatched profile/geometry pair.
	// Fail closed through the driver's existing diagnostic mechanism — no termination.
	if (ctx.geo && !ctx.geo->Matches(*ctx.star))
	{
		d.ok = false;
		d.message = "GeometryCache provenance does not match ctx.star (different StarProfile "
					"or revision). Rebuild the GeometryCache from the current context "
					"(ADR-0003).";
		return d;
	}

	// Star-integrated heat capacity cache: C(Tinf)
	d.C_eff_erg_K = ctx.star->HeatCapacityStar_Tinf(Tinf_MeV, *ctx.thermo, ctx.geo);

	if (!(d.C_eff_erg_K > 0.0) || !std::isfinite(d.C_eff_erg_K))
	{
		d.ok = false;
		d.message = "C_eff <= 0 or non-finite.";
		return d;
	}

	// ------------------------------------------------------------
	// 4) (Star context was already required above, before first use.)
	// ------------------------------------------------------------

	// ------------------------------------------------------------
	// 5) Fetch cache payload (rebuilt only if profile version changed)
	// ------------------------------------------------------------
	const auto &payload = drv.Cache_(ctx);

	if (!payload.valid)
	{
		d.ok = false;
		d.message = "NeutrinoCooling cache payload invalid.";
		return d;
	}

	d.ok = true;
	d.has_structure = true;
	d.n_zones = payload.n_zones;

	// ------------------------------------------------------------
	// 6) Evaluate luminosities using cached coefficients
	// ------------------------------------------------------------
	double L_DU = 0.0;
	double L_MU = 0.0;
	double L_PBF = 0.0;

	// DUrca ~ T^6
	if (opt.include_direct_urca && payload.K_DU_erg_s_K6 > 0.0)
		L_DU = payload.K_DU_erg_s_K6 * PowInt(d.Tinf_K, 6);

	// MUrca ~ T^8
	if (opt.include_modified_urca && payload.K_MU_erg_s_K8 > 0.0)
		L_MU = payload.K_MU_erg_s_K8 * PowInt(d.Tinf_K, 8);

	// PBF (hook; keep zero unless implemented)
	if (opt.include_pair_breaking && payload.K_PBF > 0.0)
	{
		// If/when we implement, define exponent and compute here.
		// Example: L_PBF = payload.K_PBF * PowInt(d.Tinf_K, 7);
	}

	// Apply global scaling at the end (keeps payload purely structural).
	L_DU *= gscale;
	L_MU *= gscale;
	L_PBF *= gscale;

	// Guard against non-finite
	if (!std::isfinite(L_DU))
		L_DU = 0.0;
	if (!std::isfinite(L_MU))
		L_MU = 0.0;
	if (!std::isfinite(L_PBF))
		L_PBF = 0.0;

	d.L_nu_DU_inf_erg_s = L_DU;
	d.L_nu_MU_inf_erg_s = L_MU;
	d.L_nu_PBF_inf_erg_s = L_PBF;

	d.L_nu_inf_erg_s = L_DU + L_MU + L_PBF;

	// ------------------------------------------------------------
	// 7) Cooling rate
	// ------------------------------------------------------------
	d.dTinf_dt_K_s = -d.L_nu_inf_erg_s / d.C_eff_erg_K;
	d.dLnTinf_dt_1_s = d.dTinf_dt_K_s / d.Tinf_K;

	return d;
}

// -----------------------------------------------------------------------------
//  Diagnostics interface
// -----------------------------------------------------------------------------
void Diagnose(const NeutrinoCooling &self,
			  double t,
			  const Evolution::StateVector &Y,
			  const Evolution::DriverContext &ctx,
			  Evolution::Diagnostics::DiagnosticPacket &out)
{
	PROFILE_FUNCTION();

	out.SetProducer(self.DiagnosticsName());
	out.SetTime(t);

	const NeutrinoCooling_Details d = NeutrinoCooling_Details::ComputeDerived(self, Y, ctx);

	if (!d.ok)
		out.AddWarning("NeutrinoCooling details not OK: " + d.message);
	else if (!d.message.empty())
		out.AddNote(d.message);

	out.AddScalar("Tinf_K", d.Tinf_K, "K",
				  "Redshifted internal temperature (evolved DOF)", "state");

	out.AddScalar("L_nu_inf_erg_s", d.L_nu_inf_erg_s, "erg/s",
				  "Total neutrino luminosity at infinity", "computed");

	out.AddScalar("L_nu_DU_inf_erg_s", d.L_nu_DU_inf_erg_s, "erg/s",
				  "Direct Urca neutrino luminosity at infinity", "computed",
				  Evolution::Diagnostics::Cadence::OnChange);

	out.AddScalar("L_nu_MU_inf_erg_s", d.L_nu_MU_inf_erg_s, "erg/s",
				  "Modified Urca neutrino luminosity at infinity", "computed",
				  Evolution::Diagnostics::Cadence::OnChange);

	out.AddScalar("L_nu_PBF_inf_erg_s", d.L_nu_PBF_inf_erg_s, "erg/s",
				  "Pair breaking/formation neutrino luminosity at infinity", "computed",
				  Evolution::Diagnostics::Cadence::OnChange);

	out.AddScalar("dTinf_dt_K_s", d.dTinf_dt_K_s, "K/s",
				  "NeutrinoCooling contribution to dTinf/dt", "computed");

	out.AddScalar("dLnTinf_dt_1_s", d.dLnTinf_dt_1_s, "1/s",
				  "NeutrinoCooling contribution to d/dt ln(Tinf/Tref)", "computed");

	out.ValidateBasic();
}
// -----------------------------------------------------------------------------
} // namespace Detail

// -----------------------------------------------------------------------------
//  Cache accessor
const NeutrinoCoolingCachePayload &
NeutrinoCooling::Cache_(const Evolution::DriverContext &ctx) const
{
	// ADR-0003: the payload depends on the geometry as well as on the profile. The generic
	// ProfileVersionedCache tracks profile provenance only, so the geometry half of the
	// validity condition is enforced here: if the geometry snapshot backing the cached
	// coefficients is no longer the one supplied, drop the payload before asking for it.
	const Evolution::ProfileProvenance geo_prov =
		ctx.geo ? ctx.geo->Provenance() : Evolution::ProfileProvenance{};
	if (geo_prov != cached_geo_prov_)
	{
		cache_.Invalidate();
		cached_geo_prov_ = geo_prov;
	}

	// Builder lambda lives in this .cpp.
	return cache_.Get(*ctx.star, [&](const Evolution::StarContext &sc,
									 NeutrinoCoolingCachePayload &out)
					  {
						  Detail::BuildNeutrinoCoolingCache(sc, ctx, out); // free function in NeutrinoCooling_Details.cpp
					  });
}

// -----------------------------------------------------------------------------
} // namespace CompactStar::Physics::Driver::Thermal