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
 * @file NeutrinoCooling.cpp
 * @brief Implementation of NeutrinoCooling driver (core neutrino cooling).
 *
 * @ingroup PhysicsDriver
 *
 * This file implements:
 * - NeutrinoCooling::AccumulateRHS
 * - IDriverDiagnostics interface methods
 *
 * ## Design policy
 * This driver mirrors PhotonCooling's architecture:
 * - The *only* physics done here is to call `Detail::ComputeDerived()` and add
 *   the derived RHS contribution.
 * - Diagnostics are produced by calling `Detail::Diagnose()` to prevent drift.
 *
 * ## Why derived helpers matter
 * Neutrino cooling will evolve quickly (DUrca thresholds, composition masks,
 * superfluid suppression, PBF, etc.). If you compute `L_nu` one way in physics
 * and re-implement it differently in diagnostics, they will diverge.
 *
 * Centralizing the model in the details layer ensures:
 * - a single definition of `L_nu_inf`, `dT/dt`, and `d/dt ln(Tinf/Tref)`,
 * - consistent output between RHS and packets,
 * - easier maintenance as microphysics matures.
 */

#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp"

// Diagnostics schema/packet types
#include "CompactStar/Physics/Evolution/Diagnostics/DiagnosticCatalog.hpp"
#include "CompactStar/Physics/Evolution/Diagnostics/DiagnosticPacket.hpp"

// Details helpers (ComputeDerived, Diagnose)
#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling_Details.hpp"

#include <cmath>
#include <limits>

#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/RHSAccumulator.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

#include <Zaki/Util/Instrumentor.hpp> // PROFILE_FUNCTION
#include <Zaki/Util/Logger.hpp>		  // Z_LOG_INFO/WARNING/ERROR

namespace CompactStar::Physics::Driver::Thermal
{
// -----------------------------------------------------------------------------
NeutrinoCooling::NeutrinoCooling(const Options &opts)
	: opts_(opts)
{
}
// -----------------------------------------------------------------------------
NeutrinoCooling::NeutrinoCooling(Options &&opts)
	: opts_(std::move(opts))
{
}

// -----------------------------------------------------------------------------
//  NeutrinoCooling::AccumulateRHS
// -----------------------------------------------------------------------------
void NeutrinoCooling::AccumulateRHS(double t,
									const Evolution::StateVector &Y,
									Evolution::RHSAccumulator &dYdt,
									const Evolution::DriverContext &ctx) const
{
	PROFILE_FUNCTION();

	// The current implementation has no explicit time dependence.
	(void)t;

	// -------------------------------------------------------------------------
	//  Centralized computation (shared with diagnostics)
	// -------------------------------------------------------------------------
	// Compute all quantities through the details helper to prevent drift.
	const auto d = Detail::NeutrinoCooling_Details::ComputeDerived(*this, Y, ctx);

	// If the helper cannot compute a well-defined contribution, do nothing.
	// This should not be treated as an error in infrastructure-first runs.
	if (!d.ok)
	{
		// Optional debug:
		// Z_LOG_WARNING("NeutrinoCooling skipped: " + d.message);
		return;
	}

	// Informational notes can be logged if desired (kept silent by default).
	// if (!d.message.empty()) { Z_LOG_INFO("NeutrinoCooling note: " + d.message); }

	// -------------------------------------------------------------------------
	//  Accumulate RHS contribution
	// -------------------------------------------------------------------------
	// CompactStar evolves x = ln(Tinf/Tref). The Details layer should provide:
	//   dLnTinf_dt_1_s = dx/dt  [1/s]
	if (Y.GetThermal().Size() == 0)
		return;

	// Defensive: avoid poisoning RHS with NaN/Inf.
	if (!std::isfinite(d.dLnTinf_dt_1_s))
		return;

	dYdt.AddTo(State::StateTag::Thermal, 0, d.dLnTinf_dt_1_s);
}

// -----------------------------------------------------------------------------
//  IDriverDiagnostics interface
// -----------------------------------------------------------------------------
std::string NeutrinoCooling::DiagnosticsName() const
{
	return "NeutrinoCooling";
}

// -----------------------------------------------------------------------------
//  UnitContract interface
// -----------------------------------------------------------------------------
Evolution::Diagnostics::UnitContract NeutrinoCooling::UnitContract() const
{
	Evolution::Diagnostics::UnitContract c;
	// Fill using your UnitContract API (left as project-specific).
	// Example intent:
	// - Tinf in K
	// - L_nu_inf in erg/s
	// - dLnTinf_dt in 1/s
	return c;
}

// -----------------------------------------------------------------------------
//  DiagnosticsCatalog interface
// -----------------------------------------------------------------------------
Evolution::Diagnostics::ProducerCatalog NeutrinoCooling::DiagnosticsCatalog() const
{
	using namespace CompactStar::Physics::Evolution::Diagnostics;

	ProducerCatalog pc;
	pc.producer = DiagnosticsName();

	// -------------------------------------------------------------------------
	// Scalars emitted by Detail::Diagnose() for this driver
	// -------------------------------------------------------------------------
	{
		ScalarDescriptor sd;
		sd.key = "Tinf_K";
		sd.unit = "K";
		sd.description = "Redshifted internal temperature (evolved DOF)";
		sd.source_hint = "state";
		sd.default_cadence = Cadence::Always;
		sd.required = true;
		pc.scalars.push_back(sd);
	}

	{
		ScalarDescriptor sd;
		sd.key = "L_nu_inf_erg_s";
		sd.unit = "erg/s";
		sd.description = "Total neutrino luminosity at infinity";
		sd.source_hint = "computed";
		sd.default_cadence = Cadence::Always; // or OnChange if you prefer
		sd.required = false;
		pc.scalars.push_back(sd);
	}

	// Optional partitions (keep if your Details layer emits them)
	{
		ScalarDescriptor sd;
		sd.key = "L_nu_DU_inf_erg_s";
		sd.unit = "erg/s";
		sd.description = "Direct Urca neutrino luminosity at infinity";
		sd.source_hint = "computed";
		sd.default_cadence = Cadence::OnChange;
		sd.required = false;
		pc.scalars.push_back(sd);
	}

	{
		ScalarDescriptor sd;
		sd.key = "L_nu_MU_inf_erg_s";
		sd.unit = "erg/s";
		sd.description = "Modified Urca neutrino luminosity at infinity";
		sd.source_hint = "computed";
		sd.default_cadence = Cadence::OnChange;
		sd.required = false;
		pc.scalars.push_back(sd);
	}

	{
		ScalarDescriptor sd;
		sd.key = "L_nu_PBF_inf_erg_s";
		sd.unit = "erg/s";
		sd.description = "Pair breaking/formation neutrino luminosity at infinity";
		sd.source_hint = "computed";
		sd.default_cadence = Cadence::OnChange;
		sd.required = false;
		pc.scalars.push_back(sd);
	}

	{
		ScalarDescriptor sd;
		sd.key = "dTinf_dt_K_s";
		sd.unit = "K/s";
		sd.description = "NeutrinoCooling contribution to dTinf/dt";
		sd.source_hint = "computed";
		sd.default_cadence = Cadence::Always;
		sd.required = false;
		pc.scalars.push_back(sd);
	}

	{
		ScalarDescriptor sd;
		sd.key = "dLnTinf_dt_1_s";
		sd.unit = "1/s";
		sd.description = "NeutrinoCooling contribution to d/dt ln(Tinf/Tref)";
		sd.source_hint = "computed";
		sd.default_cadence = Cadence::Always;
		sd.required = false;
		pc.scalars.push_back(sd);
	}

	// -------------------------------------------------------------------------
	// Profiles (schema-driven time series)
	// -------------------------------------------------------------------------
	ProducerCatalog::Profile p;
	p.name = "timeseries_default";
	p.keys = {
		// "Tinf_K",
		"L_nu_inf_erg_s",
		// keep partitions if you want them in the default CSV:
		// "L_nu_DU_inf_erg_s",
		// "L_nu_MU_inf_erg_s",
		// "L_nu_PBF_inf_erg_s",
		"dLnTinf_dt_1_s"};
	pc.profiles.push_back(std::move(p));

	return pc;
}

// -----------------------------------------------------------------------------
//  DiagnoseSnapshot interface
// -----------------------------------------------------------------------------
void NeutrinoCooling::DiagnoseSnapshot(double t,
									   const Evolution::StateVector &Y,
									   const Evolution::DriverContext &ctx,
									   Evolution::Diagnostics::DiagnosticPacket &out) const
{
	// Overwrite/fill `out` deterministically using the shared details path.
	Detail::Diagnose(*this, t, Y, ctx, out);
}
// -----------------------------------------------------------------------------
} // namespace CompactStar::Physics::Driver::Thermal