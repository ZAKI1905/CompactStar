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
 * @file PhotonCooling_Details.cpp
 * @brief Out-of-line implementation of PhotonCooling shared computations.
 *
 * Rationale:
 *  - Avoid heavy header-only code that forces includes everywhere.
 *  - Allow dereferencing ctx.geo (GeometryCache) without incomplete-type issues.
 *  - Keep physics and diagnostics consistent by centralizing shared math.
 */

#include "CompactStar/Physics/Driver/Thermal/PhotonCooling_Details.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp"

#include "CompactStar/Physics/Driver/Thermal/Boundary/EnvelopePotekhin1997.hpp" // optional
#include "CompactStar/Physics/Driver/Thermal/Boundary/EnvelopePotekhin2003.hpp" // or 1997 if we want
#include "CompactStar/Physics/Driver/Thermal/Boundary/SurfaceGravity.hpp"
#include "CompactStar/Physics/Driver/Thermal/Boundary/TbDefinition.hpp"

#include <cmath>
#include <string>

#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

// Diagnostics packet API (whatever you already created)
#include "CompactStar/Physics/Evolution/Diagnostics/DiagnosticPacket.hpp"

#include <Zaki/Util/Instrumentor.hpp>
#include <Zaki/Util/Logger.hpp>

namespace CompactStar::Physics::Driver::Thermal::Detail
{

namespace
{
/// Stefan–Boltzmann constant in cgs units: erg cm^-2 s^-1 K^-4.
constexpr double SigmaSB_cgs = 5.670374419e-5;

/// km -> cm conversion for radius.
constexpr double KM_TO_CM = 1.0e5;
} // namespace

// -----------------------------------------------------------------------------
//  ComputeDerived
// -----------------------------------------------------------------------------
PhotonCooling_Details ComputeDerived(const PhotonCooling &drv,
									 const Evolution::StateVector &Y,
									 const Evolution::DriverContext &ctx)
{
	PhotonCooling_Details d;

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
		d.message = "Tinf <= 0; photon cooling ill-defined.";
		return d;
	}

	// ---------------------------------------------------------------------
	// 2) Surface temperature mapping
	// ---------------------------------------------------------------------
	switch (drv.GetOptions().surface_model)
	{
	case PhotonCooling::Options::SurfaceModel::DirectTSurf:
		d.Tsurf_K = (thermal.T_surf > 0.0) ? thermal.T_surf : d.Tinf_K;
		break;

	// case PhotonCooling::Options::SurfaceModel::EnvelopeTbTs:
	// {
	// 	// We need star context to locate rho_b and compute Tb, and we need star/geo to compute g14.
	// 	if (!ctx.star)
	// 	{
	// 		d.ok = false;
	// 		d.message = "SurfaceModel::EnvelopeTbTs selected but ctx.star == nullptr.";
	// 		return d;
	// 	}

	// 	// 1) Define the Tb boundary policy from driver options.
	// 	Boundary::TbDefinition def;
	// 	def.rho_b = drv.GetOptions().rho_b;		 // g/cm^3
	// 	def.assume_isothermal_redshifted = true; // Tb derived from T_inf via redshift
	// 	def.prefer_geometry_cache = true;		 // use geo->ExpNu if present
	// 	// def.prefer_geometry_cache remains true even if ctx.geo==nullptr; ComputeTb will fall back to star.Nu()

	// 	// 2) Compute local Tb at rho_b.
	// 	const double Tb_K = Boundary::ComputeTb(*ctx.star, ctx.geo, d.Tinf_K, def);
	// 	if (!(Tb_K > 0.0))
	// 	{
	// 		d.ok = false;
	// 		d.message = "EnvelopeTbTs: computed Tb <= 0.";
	// 		return d;
	// 	}

	// 	// 3) Compute g14 at the surface (dimensionless).
	// 	// If you want EnvelopeTbTs to REQUIRE geo, make SurfaceGravity_g14 strict internally.
	// 	const double g14 = Boundary::SurfaceGravity_g14(*ctx.star, ctx.geo);
	// 	if (!(g14 > 0.0))
	// 	{
	// 		d.ok = false;
	// 		d.message = "EnvelopeTbTs: computed g14 <= 0.";
	// 		return d;
	// 	}

	// 	// 4) Apply the selected envelope fit.
	// 	const double xi = drv.GetOptions().envelope_xi;

	// 	switch (drv.GetOptions().envelope)
	// 	{
	// 	case PhotonCooling::EnvelopeModel::Iron:
	// 		d.Tsurf_K = Boundary::EnvelopePotekhin2003_Iron{}.Ts_from_Tb(Tb_K, g14, xi);
	// 		break;

	// 	case PhotonCooling::EnvelopeModel::Accreted:
	// 		d.Tsurf_K = Boundary::EnvelopePotekhin2003_Accreted{}.Ts_from_Tb(Tb_K, g14, xi);
	// 		break;

	// 	case PhotonCooling::EnvelopeModel::Custom:
	// 	default:
	// 		d.ok = false;
	// 		d.message = "EnvelopeTbTs: EnvelopeModel::Custom selected but no custom mapping is wired.";
	// 		return d;
	// 	}

	// 	break;
	// }
	case PhotonCooling::Options::SurfaceModel::EnvelopeTbTs:
	{
		if (!ctx.star)
		{
			d.ok = false;
			d.message = "SurfaceModel::EnvelopeTbTs selected but ctx.star == nullptr.";
			return d;
		}

		// 1) Tb boundary policy (rho_b in g/cm^3)
		Boundary::TbDefinition def;
		def.rho_b = drv.GetOptions().rho_b;
		def.assume_isothermal_redshifted = true;
		def.prefer_geometry_cache = true;

		// 2) Compute local Tb at rho_b
		const double Tb_K = Boundary::ComputeTb(*ctx.star, ctx.geo, d.Tinf_K, def);
		if (!(Tb_K > 0.0))
		{
			d.ok = false;
			d.message = "EnvelopeTbTs: computed Tb <= 0.";
			return d;
		}
		d.Tb_K = Tb_K;

		// 3) Compute g14 at the surface
		const double g14 = Boundary::SurfaceGravity_g14(*ctx.star, ctx.geo);
		if (!(g14 > 0.0))
		{
			d.ok = false;
			d.message = "EnvelopeTbTs: computed g14 <= 0.";
			return d;
		}
		d.g14 = g14;

		// 4) Apply envelope fit
		const double xi = drv.GetOptions().envelope_xi;

		switch (drv.GetOptions().envelope)
		{
		case PhotonCooling::EnvelopeModel::Iron:
		{
			static const Boundary::EnvelopePotekhin2003_Iron model{};
			d.Tsurf_K = model.Ts_from_Tb(Tb_K, g14, xi);
			break;
		}
		case PhotonCooling::EnvelopeModel::Accreted:
		{
			static const Boundary::EnvelopePotekhin2003_Accreted model{};
			d.Tsurf_K = model.Ts_from_Tb(Tb_K, g14, xi);
			break;
		}
		case PhotonCooling::EnvelopeModel::Custom:
		default:
			d.ok = false;
			d.message = "EnvelopeTbTs: EnvelopeModel::Custom selected but no custom mapping is wired.";
			return d;
		}

		break;
	}

	case PhotonCooling::Options::SurfaceModel::ApproxFromTinf:
	default:
		d.Tsurf_K = d.Tinf_K;
		break;
	}

	if (!(d.Tsurf_K > 0.0))
	{
		d.ok = false;
		d.message = "Tsurf <= 0 after mapping.";
		return d;
	}

	// ---------------------------------------------------------------------
	// 3) Validate options used in the RHS
	//
	// The heat capacity is NOT a driver option (ADR-0002): it is the canonical
	// C_*(T_inf) obtained from StarContext below, after the disable checks, so a
	// deliberately disabled driver still needs no thermodynamic data.
	// ---------------------------------------------------------------------

	// if (!(drv.GetOptions().radiating_fraction > 0.0))
	// {
	// 	d.ok = false;
	// 	d.message = "radiating_fraction <= 0 (cooling disabled).";
	// 	return d;
	// }
	if (!(drv.GetOptions().radiating_fraction > 0.0))
	{
		// Cooling explicitly disabled: this is a valid configuration,
		// not a failure condition.
		d.ok = true;
		d.message = "radiating_fraction <= 0 (photon cooling disabled).";

		// Ensure zero contribution
		d.L_gamma_inf_erg_s = 0.0;
		d.dLnTinf_dt_1_s = 0.0;
		d.dTinf_dt_K_s = 0.0;
		d.A_inf_cm2 = 0.0;
		d.A_eff_inf_cm2 = 0.0;
		return d;
	}

	// Disabled but valid: no cooling contribution.
	if (!(drv.GetOptions().global_scale > 0.0))
	{
		d.ok = true;
		d.message = "cooling disabled: global_scale <= 0.";
		d.L_gamma_inf_erg_s = 0.0;
		d.dLnTinf_dt_1_s = 0.0;
		d.dTinf_dt_K_s = 0.0;
		d.A_inf_cm2 = 0.0;
		d.A_eff_inf_cm2 = 0.0;
		return d;
	}

	// ---------------------------------------------------------------------
	// 4) Geometry / redshifted area at infinity (STRICT)
	// ---------------------------------------------------------------------
	if (!ctx.geo)
	{
		d.ok = false;
		d.message = "ctx.geo == nullptr (GeometryCache required for A_inf).";
		return d;
	}

	const auto &R = ctx.geo->R();
	const auto &e2nu = ctx.geo->Exp2Nu();

	if (R.Size() == 0 || e2nu.Size() == 0)
	{
		d.ok = false;
		d.message = "GeometryCache arrays empty (R or Exp2Nu size is 0).";
		return d;
	}

	// Surface values (last grid point)
	d.R_surf_km = R[-1];
	d.exp2nu_surf = e2nu[-1];

	if (!(d.R_surf_km > 0.0))
	{
		d.ok = false;
		d.message = "Invalid surface radius: R_surf_km <= 0.";
		return d;
	}

	if (!(d.exp2nu_surf > 0.0))
	{
		d.ok = false;
		d.message = "Invalid surface redshift factor: exp2nu_surf <= 0.";
		return d;
	}

	d.R_surf_cm = d.R_surf_km * KM_TO_CM;
	d.A_inf_cm2 = 4.0 * M_PI * d.R_surf_cm * d.R_surf_cm * d.exp2nu_surf;
	d.A_eff_inf_cm2 = drv.GetOptions().radiating_fraction * d.A_inf_cm2;
	// ---------------------------------------------------------------------
	// 5) Luminosity and RHS term
	// ---------------------------------------------------------------------
	const double T4 = std::pow(d.Tsurf_K, 4.0);

	d.L_gamma_inf_erg_s =
		drv.GetOptions().global_scale * d.A_eff_inf_cm2 * SigmaSB_cgs * T4;

	// ---------------------------------------------------------------------
	//  Canonical heat capacity C_*(T_inf) — ADR-0002.
	//
	//  The photon channel divides by exactly the same GR-integrated, EOS-based
	//  stellar heat capacity that NeutrinoCooling uses, so both contributions sum
	//  into one physically coherent energy equation. This is ADR-0002 Pattern A:
	//  additive drivers sharing one denominator. Verified in
	//  docs/validation/HEAT_CAPACITY_V1.md (V1 VERIFIED).
	//
	//  Fails closed: an actively cooling driver must have a valid star and
	//  thermodynamics. Reached only after the disable checks above.
	// ---------------------------------------------------------------------
	if (!ctx.star)
	{
		d.ok = false;
		d.message = "ctx.star == nullptr (StarContext required for C_*(T_inf)).";
		return d;
	}
	if (!ctx.thermo)
	{
		d.ok = false;
		d.message = "ctx.thermo == nullptr (CompOSE_Thermo required for C_*(T_inf)).";
		return d;
	}

	// Same K -> MeV convention as NeutrinoCooling_Details.cpp (INV-02 records the
	// duplicated-constant debt; consolidating it is separate, later work).
	constexpr double MEV_PER_K = 8.617333262145e-11;
	const double Tinf_MeV = d.Tinf_K * MEV_PER_K;
	if (!(Tinf_MeV > 0.0) || !std::isfinite(Tinf_MeV))
	{
		d.ok = false;
		d.message = "Tinf_MeV <= 0 or non-finite after conversion.";
		return d;
	}

	d.C_star_erg_K = ctx.star->HeatCapacityStar_Tinf(Tinf_MeV, *ctx.thermo, ctx.geo);
	if (!(d.C_star_erg_K > 0.0) || !std::isfinite(d.C_star_erg_K))
	{
		d.ok = false;
		d.message = "C_*(T_inf) <= 0 or non-finite.";
		return d;
	}

	d.dTinf_dt_K_s = -d.L_gamma_inf_erg_s / d.C_star_erg_K;

	// Convert physical dT/dt to log-variable RHS (ODE variable):
	// d/dt ln(Tinf/Tref) = (1/Tinf) dTinf/dt.  (Tref constant cancels)
	d.dLnTinf_dt_1_s = d.dTinf_dt_K_s / d.Tinf_K;

	return d;
}
// -----------------------------------------------------------------------------
//  Diagnostics interface
// -----------------------------------------------------------------------------
void Diagnose(const PhotonCooling &self,
			  double t,
			  const Evolution::StateVector &Y,
			  const Evolution::DriverContext &ctx,
			  Evolution::Diagnostics::DiagnosticPacket &out)
{
	PROFILE_FUNCTION();

	// Reset packet identity/time (packet may be reused by observer).
	out.SetProducer(self.DiagnosticsName()); // or "PhotonCooling" if you want it fixed
	out.SetTime(t);

	// Compute the shared derived quantities (physics-consistent).
	const PhotonCooling_Details d = ComputeDerived(self, Y, ctx);

	// Status / messages
	if (!d.ok)
	{
		out.AddWarning("PhotonCooling details not OK: " + d.message);
		// We still dump what we have, but caller can decide what to do.
	}
	else if (!d.message.empty())
	{
		out.AddNote(d.message);
	}

	// If geometry cache missing, make it explicit (helps debugging wiring).
	if (!ctx.geo)
	{
		out.AddWarning("ctx.geo == nullptr (GeometryCache required; photon cooling not computed).");
	}

	// ---------------------------------------------------------------------
	// Scalars (key + value + unit + description + source)
	// ---------------------------------------------------------------------
	out.AddScalar("Tinf_K", d.Tinf_K, "K",
				  "Redshifted internal temperature (evolved DOF)", "state");

	out.AddScalar("Tsurf_K", d.Tsurf_K, "K",
				  "Surface temperature used in photon luminosity", "computed");

	out.AddScalar("Tb_K", d.Tb_K, "K",
				  "Local temperature at the base of the envelope (Tb) used for Tb→Ts mapping",
				  "computed");

	out.AddScalar("g14", d.g14, "",
				  "Surface gravity in units of 1e14 cm s^-2 used in envelope fit",
				  "computed",
				  Evolution::Diagnostics::Cadence::OncePerRun);

	out.AddScalar("R_surf_km", d.R_surf_km, "km",
				  "Surface radius from GeometryCache (last grid point)", "cache",
				  Evolution::Diagnostics::Cadence::OncePerRun);

	out.AddScalar("R_surf_cm", d.R_surf_cm, "cm",
				  "Surface radius converted to cgs", "computed",
				  Evolution::Diagnostics::Cadence::OncePerRun);

	out.AddScalar("exp2nu_surf", d.exp2nu_surf, "",
				  "exp(2 nu) at surface from GeometryCache", "cache",
				  Evolution::Diagnostics::Cadence::OncePerRun);

	out.AddScalar("A_inf_cm2", d.A_inf_cm2, "cm^2",
				  "Redshifted emitting area at infinity", "computed",
				  Evolution::Diagnostics::Cadence::OncePerRun);

	out.AddScalar("A_eff_inf_cm2", d.A_eff_inf_cm2, "cm^2",
				  "Effective emitting area at infinity (radiating_fraction * A_inf)", "computed",
				  Evolution::Diagnostics::Cadence::OncePerRun);

	out.AddScalar("C_star_erg_K", d.C_star_erg_K, "erg/K",
				  "Canonical GR-integrated stellar heat capacity used as the PhotonCooling denominator");

	out.AddScalar("L_gamma_inf_erg_s", d.L_gamma_inf_erg_s, "erg/s",
				  "Photon luminosity at infinity", "computed");

	out.AddScalar("dTinf_dt_K_s", d.dTinf_dt_K_s, "K/s",
				  "PhotonCooling contribution to dTinf/dt", "computed");

	out.AddScalar("dLnTinf_dt_1_s", d.dLnTinf_dt_1_s, "1/s",
				  "PhotonCooling contribution to d/dt ln(Tinf/Tref)", "computed");

	// ---------------------------------------------------------------------
	// Contracts: either the observer attaches UnitContract(), or do it here.
	// I recommend: observer attaches contracts for all drivers uniformly.
	// If you want it here:
	// for (const auto &line : self.UnitContract().Lines()) out.AddContractLine(line);
	// ---------------------------------------------------------------------

	out.ValidateBasic();
}
// -----------------------------------------------------------------------------
} // namespace CompactStar::Physics::Driver::Thermal::Detail