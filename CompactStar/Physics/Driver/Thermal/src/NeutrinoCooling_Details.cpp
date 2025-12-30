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
namespace CompactStar::Physics::Driver::Thermal::Detail
{

namespace
{

// -----------------------------------------------------------------------------
// Unit conversions
// -----------------------------------------------------------------------------
constexpr double KM3_TO_CM3 = 1.0e15; // (1e5 cm)^3

// -----------------------------------------------------------------------------
// Convert eps [km^-2] -> rho [g/cm^3] via MeV fm^-3 bridge
// -----------------------------------------------------------------------------
// inline double RhoFromEps_km2(double eps_km2) noexcept
// {
// 	// eps(MeV fm^-3) = eps(km^-2) / (MeV fm^-3 -> km^-2)
// 	const double eps_MeV_fm3 = eps_km2 / Zaki::Physics::MEV_FM3_2_INV_KM2;
// 	// rho(g/cm^3) = eps(MeV fm^-3) * (MeV fm^-3 -> g/cm^3)
// 	return eps_MeV_fm3 * Zaki::Physics::MEV_FM3_2_G_CM3;
// }
// -----------------------------------------------------------------------------
// /**
//  * @brief Minimal placeholder emissivity model for wiring tests.
//  *
//  * This exists only so the base structure (catalog/packet/timeseries) can run
//  * before you plug in real emissivities.
//  *
//  * Replace with: DUrca/MUrca/PBF emissivities integrated over proper volume.
//  *
//  * @param Tinf_K Redshifted interior temperature [K].
//  * @return A tiny luminosity scale [erg/s] for non-zero RHS during smoke tests.
//  */
// double PlaceholderLuminosity(double Tinf_K)
// {
// 	// Deliberately small-ish for stability; choose any smooth T dependence.
// 	// Typical neutrino luminosities can be huge; this is not physical.
// 	// L ~ 1e30 * (T/1e8)^6 erg/s as a harmless placeholder.
// 	const double T8 = Tinf_K / 1.0e8;
// 	return 1.0e30 * std::pow(T8, 6.0);
// }

// -----------------------------------------------------------------------------
// Minimal baseline emissivities (replace coefficients with your chosen fits)
// Returns Q in erg cm^-3 s^-1.
// -----------------------------------------------------------------------------
double Q_ModifiedUrca_erg_cm3_s(double T_local_K, double rho_g_cm3) noexcept
{
	// MUrca ~ T^8; include mild density dependence so profiles matter.
	const double T9 = T_local_K / 1.0e9;
	const double rho15 = rho_g_cm3 / 1.0e15;

	// Replace this coefficient with your preferred Yakovlev/Pethick fit.
	// Typical order for MUrca is ~ 1e21 * (rho/1e15) * T9^8 erg cm^-3 s^-1.
	const double Q0 = 1.0e21;

	// Fast T^8 via squaring
	const double T9_2 = T9 * T9;
	const double T9_4 = T9_2 * T9_2;
	const double T9_8 = T9_4 * T9_4;

	// return Q0 * rho15 * std::pow(T9, 8.0);
	return Q0 * rho15 * T9_8;
}
// -----------------------------------------------------------------------------
double Q_DirectUrca_erg_cm3_s(double T_local_K, double rho_g_cm3) noexcept
{
	// DUrca ~ T^6; much stronger normalization.
	const double T9 = T_local_K / 1.0e9;
	const double rho15 = rho_g_cm3 / 1.0e15;

	// Replace with your DUrca fit + threshold logic (proton fraction etc).
	const double Q0 = 1.0e27;

	// Fast T^6 = T^2 * T^4
	const double T9_2 = T9 * T9;
	const double T9_4 = T9_2 * T9_2;
	const double T9_6 = T9_4 * T9_2;

	// return Q0 * rho15 * std::pow(T9, 6.0);
	return Q0 * rho15 * T9_6;
}
// -----------------------------------------------------------------------------
struct TPowers9
{
	double T9; // optional
	double T9_2;
	double T9_4;
	double T9_6;
	double T9_8;
};
// -----------------------------------------------------------------------------
inline TPowers9 MakeTPowers9(double T_local_K) noexcept
{
	const double T9 = T_local_K * 1.0e-9;
	const double T9_2 = T9 * T9;
	const double T9_4 = T9_2 * T9_2;

	TPowers9 p;
	p.T9 = T9;
	p.T9_2 = T9_2;
	p.T9_4 = T9_4;
	p.T9_6 = T9_4 * T9_2;
	p.T9_8 = T9_4 * T9_4;
	return p;
}
// -----------------------------------------------------------------------------
inline double Q_ModifiedUrca_from_T9powers(double rho_g_cm3, const TPowers9 &p) noexcept
{
	constexpr double Q0 = 1.0e21;
	const double rho15 = rho_g_cm3 * 1.0e-15;
	return Q0 * rho15 * p.T9_8;
}
// -----------------------------------------------------------------------------
inline double Q_DirectUrca_from_T9powers(double rho_g_cm3, const TPowers9 &p) noexcept
{
	constexpr double Q0 = 1.0e27;
	const double rho15 = rho_g_cm3 * 1.0e-15;
	return Q0 * rho15 * p.T9_6;
}

// -----------------------------------------------------------------------------
} // namespace

// -----------------------------------------------------------------------------
//  ComputeDerived
// -----------------------------------------------------------------------------
NeutrinoCooling_Details ComputeDerived(const NeutrinoCooling &drv,
									   const Evolution::StateVector &Y,
									   const Evolution::DriverContext &ctx)
{
	PROFILE_FUNCTION();

	NeutrinoCooling_Details d;

	// ---------------------------------------------------------------------
	// 1) Extract evolved DOF: T_inf
	// ---------------------------------------------------------------------
	const auto &thermal = Y.GetThermal();
	if (thermal.NumComponents() == 0)
	{
		d.ok = false;
		d.message = "ThermalState has zero components.";
		return d;
	}

	d.Tinf_K = thermal.Tinf();
	if (!(d.Tinf_K > 0.0))
	{
		d.ok = false;
		d.message = "Tinf <= 0; neutrino cooling ill-defined.";
		return d;
	}

	// ---------------------------------------------------------------------
	// 2) Resolve/validate options used in RHS
	// ---------------------------------------------------------------------
	// Your NeutrinoCooling::Options currently includes global_scale and channel toggles.
	// You also need an effective heat capacity; choose one policy:
	//   A) store C_eff inside NeutrinoCooling::Options (recommended for symmetry w/ PhotonCooling),
	//   B) fetch from ctx (if you already have a thermal capacity model),
	//   C) reuse the same C_eff as photon cooling in a shared thermal config object.
	//
	// For now, we assume NeutrinoCooling has a C_eff accessor or that you will add it.
	//
	// If you don't have it yet, set a temporary constant here and move it later.
	//
	// IMPORTANT: make this match your driver header/cpp choices.
	//
	// 2) Heat capacity policy (keep your requested constant for now)
	d.C_eff_erg_K = 1.0e40;
	if (!(d.C_eff_erg_K > 0.0))
	{
		d.ok = false;
		d.message = "C_eff <= 0.";
		return d;
	}

	if (!(drv.GetOptions().global_scale > 0.0))
	{
		// Disabled but valid: no cooling contribution.
		d.ok = true;
		d.message = "cooling disabled: global_scale <= 0.";
		d.L_nu_inf_erg_s = 0.0;
		d.dTinf_dt_K_s = 0.0;
		d.dLnTinf_dt_1_s = 0.0;
		return d;
	}

	// If all channels disabled, treat as disabled-but-valid.
	if (!drv.GetOptions().include_direct_urca &&
		!drv.GetOptions().include_modified_urca &&
		!drv.GetOptions().include_pair_breaking)
	{
		d.ok = true;
		d.message = "cooling disabled: all neutrino channels disabled by options.";
		d.L_nu_inf_erg_s = 0.0;
		d.dTinf_dt_K_s = 0.0;
		d.dLnTinf_dt_1_s = 0.0;
		return d;
	}

	// ---------------------------------------------------------------------
	// 3) Compute neutrino luminosity at infinity (placeholder vs real)
	// ---------------------------------------------------------------------
	//
	// This is where the real work will go:
	//   - loop over radial zones
	//   - compute local T(r) (isothermal redshifted => T_local(r) = Tinf * e^{-nu(r)})
	//   - compute emissivities Q_nu (erg cm^-3 s^-1) per process via microphysics
	//   - multiply by proper volume element dV_proper (cm^3)
	//   - apply redshift factors to get L_inf
	//
	// For now, implement a minimal placeholder so the infrastructure runs.
	//
	const auto &R_km = ctx.geo->R();
	const auto &e_minus_nu = ctx.geo->ExpMinusNu();
	const auto &wV_e2nu = ctx.geo->WVExp2Nu(); // 4π r^2 e^Λ e^{2ν}

	const std::size_t N = R_km.Size();
	if (N < 2)
	{
		d.ok = true;
		d.n_zones = N;
		d.L_nu_inf_erg_s = 0.0;
		d.dTinf_dt_K_s = 0.0;
		d.dLnTinf_dt_1_s = 0.0;
		return d;
	}
	if (e_minus_nu.Size() != N || wV_e2nu.Size() != N)
	{
		d.ok = false;
		d.message = "GeometryCache size mismatch (R vs ExpMinusNu vs WVExp2Nu).";
		return d;
	}

	// 4) Require structure for eps(r)
	if (!ctx.star)
	{
		d.ok = false;
		d.message = "ctx.star == nullptr (StarContext required for eps(r)).";
		return d;
	}

	// const auto *eps_km2_col = ctx.star->EnergyDensity(); // eps [km^-2]
	// if (!eps_km2_col || eps_km2_col->Size() != N)
	// {
	// 	d.ok = false;
	// 	d.message = "StarContext::EnergyDensity missing or size mismatch with geometry grid.";
	// 	return d;
	// }
	const auto *rho_g_cm3_col = ctx.star->MassDensity_gcm3(); // rho [g/cm^3]
	if (!rho_g_cm3_col || rho_g_cm3_col->Size() != N)
	{
		d.ok = false;
		d.message = "StarContext::MassDensity_gcm3 missing or size mismatch with geometry grid.";
		return d;
	}

	d.ok = true;
	d.has_structure = true;
	d.n_zones = N;

	// 5) Build integrands in (erg cm^-3 s^-1) * (km^2)?? no:
	// integrand[i] = Q(erg/cm^3/s) * WVExp2Nu(km^2?) * ... then dr(km) => km^3; convert km^3->cm^3 later
	// ---------------------------------------------------------------------
	// 5) Single-pass trapezoid integration over enabled channels
	// ---------------------------------------------------------------------

	// Hoist options locally (cuts repeated virtual/inline calls and branches).
	const auto &opt = drv.GetOptions();
	const bool do_DU = opt.include_direct_urca;
	const bool do_MU = opt.include_modified_urca;
	const bool do_PBF = opt.include_pair_breaking;
	const double gscale = opt.global_scale;

	// Fast exit if globally disabled (you already handled <=0 earlier, but keep it safe).
	if (!(gscale > 0.0))
	{
		d.L_nu_DU_inf_erg_s = 0.0;
		d.L_nu_MU_inf_erg_s = 0.0;
		d.L_nu_PBF_inf_erg_s = 0.0;
		d.L_nu_inf_erg_s = 0.0;
		d.dTinf_dt_K_s = 0.0;
		d.dLnTinf_dt_1_s = 0.0;
		return d;
	}

	// Helper lambda: compute emissivity-weighted integrand at zone i.
	// Returns f_total(i) = (Q_DU + Q_MU + Q_PBF) * wV_e2nu[i]
	// and also provides optional per-channel f's for bookkeeping.
	auto eval_f = [&](std::size_t i, double &f_DU, double &f_MU, double &f_PBF) -> double
	{
		f_DU = 0.0;
		f_MU = 0.0;
		f_PBF = 0.0;

		// Local temperature (isothermal interior w/ redshift)
		const double T_local = d.Tinf_K * e_minus_nu[i];
		if (!(T_local > 0.0) || !std::isfinite(T_local))
			return 0.0;

		// const double rho_g_cm3 = RhoFromEps_km2((*eps_km2_col)[i]);
		const double rho_g_cm3 = (*rho_g_cm3_col)[i];
		if (!(rho_g_cm3 > 0.0) || !std::isfinite(rho_g_cm3))
			return 0.0;

		const auto p = MakeTPowers9(T_local);

		// Compute only enabled channels.
		// Multiply by geometric weight once at the end.
		double Qsum = 0.0;

		if (do_DU)
		{
			const double Q = Q_DirectUrca_from_T9powers(rho_g_cm3, p);
			if (std::isfinite(Q) && Q > 0.0)
			{
				f_DU = Q;
				Qsum += Q;
			}
		}

		if (do_MU)
		{
			const double Q = Q_ModifiedUrca_from_T9powers(rho_g_cm3, p);
			if (std::isfinite(Q) && Q > 0.0)
			{
				f_MU = Q;
				Qsum += Q;
			}
		}

		if (do_PBF)
		{
			// Placeholder until implemented:
			// const double Q = Q_PBF_erg_cm3_s(T_local, rho_g_cm3, ...);
			// if (std::isfinite(Q) && Q > 0.0) { f_PBF = Q; Qsum += Q; }
			// For now, leave zero.
		}

		const double w = wV_e2nu[i];
		if (!(w > 0.0) || !std::isfinite(w))
			return 0.0;

		// Convert per-channel to fully-weighted integrands for trapezoid accumulation
		// (units: Q * w, integrated over dr in km => "km^3" factor to be converted later).
		f_DU *= w;
		f_MU *= w;
		f_PBF *= w;

		return Qsum * w;
	};

	// Trapezoid accumulators (km^3 in your current convention).
	double I_total_km3 = 0.0;
	double I_DU_km3 = 0.0;
	double I_MU_km3 = 0.0;
	double I_PBF_km3 = 0.0;

	// Evaluate at i=0 to seed trapezoid.
	double f0_DU, f0_MU, f0_PBF;
	double f0 = eval_f(0, f0_DU, f0_MU, f0_PBF);

	for (std::size_t i = 0; i + 1 < N; ++i)
	{
		const double r0 = R_km[i];
		const double r1 = R_km[i + 1];
		const double dr = r1 - r0;

		// // Skip non-monotonic or degenerate steps.
		// // (If your grid is guaranteed monotonic, you can remove this branch.)
		// if (!(dr > 0.0) || !std::isfinite(dr))
		// {
		// 	// Re-seed next point in case of weird grid.
		// 	double tmpDU, tmpMU, tmpPBF;
		// 	f0 = eval_f(i + 1, tmpDU, tmpMU, tmpPBF);
		// 	f0_DU = tmpDU;
		// 	f0_MU = tmpMU;
		// 	f0_PBF = tmpPBF;
		// 	continue;
		// }

		double f1_DU, f1_MU, f1_PBF;
		const double f1 = eval_f(i + 1, f1_DU, f1_MU, f1_PBF);

		// Standard trapezoid: 0.5*(f0+f1)*dr
		const double wtrap = 0.5 * dr;

		I_total_km3 += (f0 + f1) * wtrap;

		// If you want per-channel bookkeeping, accumulate them too.
		// This costs only a few adds/mults and avoids extra integrations.
		if (do_DU)
			I_DU_km3 += (f0_DU + f1_DU) * wtrap;
		if (do_MU)
			I_MU_km3 += (f0_MU + f1_MU) * wtrap;
		if (do_PBF)
			I_PBF_km3 += (f0_PBF + f1_PBF) * wtrap;

		// Advance
		f0 = f1;
		f0_DU = f1_DU;
		f0_MU = f1_MU;
		f0_PBF = f1_PBF;
	}

	// Convert km^3 -> cm^3, apply global scale
	d.L_nu_DU_inf_erg_s = gscale * (I_DU_km3 * KM3_TO_CM3);
	d.L_nu_MU_inf_erg_s = gscale * (I_MU_km3 * KM3_TO_CM3);
	d.L_nu_PBF_inf_erg_s = gscale * (I_PBF_km3 * KM3_TO_CM3);

	// Total luminosity (either use I_total or sum channels; both should agree if PBF implemented)
	// Using I_total is cheapest and avoids any numerical drift from per-channel guarding.
	d.L_nu_inf_erg_s = gscale * (I_total_km3 * KM3_TO_CM3);

	// If you want a consistency check in debug builds, you can do:
	// const double sum_channels = d.L_nu_DU_inf_erg_s + d.L_nu_MU_inf_erg_s + d.L_nu_PBF_inf_erg_s;
	// assert(std::abs(d.L_nu_inf_erg_s - sum_channels) / std::max(1.0, d.L_nu_inf_erg_s) < 1e-12);
	// 6) Cooling rate
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

	const NeutrinoCooling_Details d = ComputeDerived(self, Y, ctx);

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
} // namespace CompactStar::Physics::Driver::Thermal::Detail