#pragma once
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
 * @file NeutrinoCooling.hpp
 * @brief Core neutrino-cooling driver (DUrca, MUrca, PBF hooks) for thermal evolution.
 *
 * @ingroup PhysicsDriver
 *
 * ## Physical intent
 * This driver contributes a **neutrino luminosity loss** term to the redshifted
 * thermal degree of freedom (typically \(T_\infty\)). CompactStar evolves a
 * **logarithmic temperature variable**:
 *
 * \f[
 *   x \equiv \ln\!\left(\frac{T_\infty}{T_{\rm ref}}\right),
 * \f]
 *
 * so the driver must accumulate \(dx/dt\) into the thermal RHS:
 *
 * \f[
 *   \frac{dx}{dt} = \frac{1}{T_\infty}\,\frac{dT_\infty}{dt}.
 * \f]
 *
 * In a full neutrino-cooling model one would compute:
 *
 * \f[
 *   L_{\nu,\infty} = \int dV_{\rm proper}\; Q_\nu(r)\times(\text{redshift factors}),
 * \f]
 * \f[
 *   \frac{dT_\infty}{dt} = -\frac{L_{\nu,\infty}}{C_{\rm eff}},
 * \f]
 * \f[
 *   \frac{dx}{dt} = -\frac{L_{\nu,\infty}}{C_{\rm eff}\,T_\infty}.
 * \f]
 *
 * ## Implementation pattern
 * This driver follows the same architecture as PhotonCooling:
 * - Physics and diagnostics share the same computation via `NeutrinoCooling_Details`.
 * - `AccumulateRHS()` calls `Detail::ComputeDerived()` and adds only the derived
 *   \(dx/dt\) contribution.
 * - `DiagnoseSnapshot()` calls `Detail::Diagnose()` so diagnostics cannot drift
 *   from the physics implementation.
 *
 * ## Units policy
 * This header documents expected units; the implementation must remain consistent:
 * - \(T_\infty\): Kelvin [K]
 * - \(L_{\nu,\infty}\): [erg/s]
 * - \(C_{\rm eff}\): [erg/K]
 * - \(dT_\infty/dt\): [K/s]
 * - \(dx/dt\): [1/s]
 */

#include <string>
#include <vector>

#include "CompactStar/Physics/Driver/IDriver.hpp"
#include "CompactStar/Physics/State/Tags.hpp"

// IDriverDiagnostics interface (same base as PhotonCooling)
#include "CompactStar/Physics/Driver/Diagnostics/DriverDiagnostics.hpp"

#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling_Cache.hpp"
#include "CompactStar/Physics/Evolution/ProfileProvenance.hpp"

namespace CompactStar::Physics::Driver::Thermal::Detail
{
struct NeutrinoCooling_Details;
} // namespace CompactStar::Physics::Driver::Thermal::Detail

namespace CompactStar::Physics::Driver::Thermal
{

/**
 * @class NeutrinoCooling
 * @brief Adds a neutrino-cooling contribution to the thermal RHS (redshifted frame).
 *
 * **Depends on:** Thermal
 * **Updates:**    Thermal
 *
 * ### What this driver does
 * 1. Reads the evolved thermal DOF \(T_\infty\) from the ThermalState in `Y`.
 * 2. Computes derived neutrino cooling quantities through `Detail::ComputeDerived()`.
 * 3. Accumulates the derived \(dx/dt = d/dt\ln(T_\infty/T_{\rm ref})\) into the RHS.
 *
 * ### What this driver does not do (yet)
 * - It does not implement the full microphysics emissivity set by default.
 * - It does not compute \(C_{\rm eff}\) from EOS/microphysics unless you wire it.
 * - It does not modify structure (TOV/geometry profiles).
 *
 * ### Determinism contract
 * - The driver is pure with respect to `Y` and `ctx` (read-only).
 * - For fixed `(t, Y, ctx)` it produces deterministic output.
 */
class NeutrinoCooling final : public IDriver,
							  public Driver::Diagnostics::IDriverDiagnostics
{
  public:
	/**
	 * @struct Options
	 * @brief Configuration knobs controlling which neutrino channels are enabled.
	 *
	 * These options are intended to map to future microphysics channel selection.
	 * `global_scale` provides a convenient dimensionless knob for debugging/tests.
	 *
	 * @note For a physically meaningful model, you will eventually also want a policy
	 *       for \(C_{\rm eff}\) (either in Options or supplied by `ctx`).
	 */
	struct Options
	{
		/// Enable direct Urca contribution (if available in your details implementation).
		bool include_direct_urca = true;

		/// Enable modified Urca contribution (if available in your details implementation).
		bool include_modified_urca = true;

		/// Enable pair-breaking/formation contribution (future superfluid hook).
		bool include_pair_breaking = false;

		/// Dimensionless multiplicative scale applied to the net cooling rate.
		double global_scale = 1.0;
	};

	/**
	 * @brief Default constructor.
	 *
	 * Uses default Options values.
	 *
	 * @note This avoids the toolchain pitfall where `Options opts = {}` combined with
	 *       default member initializers in nested structs can fail under some configurations.
	 */
	NeutrinoCooling() = default;

	/**
	 * @brief Construct with explicit options (copy).
	 * @param opts Configuration options.
	 */
	explicit NeutrinoCooling(const Options &opts);

	/**
	 * @brief Construct with explicit options (move).
	 * @param opts Configuration options.
	 */
	explicit NeutrinoCooling(Options &&opts);

	// -------------------------------------------------------------------------
	// IDriver interface
	// -------------------------------------------------------------------------

	/// Stable driver identifier.
	[[nodiscard]] std::string Name() const override { return "NeutrinoCooling"; }

	/// State blocks this driver reads.
	[[nodiscard]] const std::vector<State::StateTag> &DependsOn() const override
	{
		static const std::vector<State::StateTag> deps{State::StateTag::Thermal};
		return deps;
	}

	/// State blocks this driver updates (accumulates into RHS).
	[[nodiscard]] const std::vector<State::StateTag> &Updates() const override
	{
		static const std::vector<State::StateTag> ups{State::StateTag::Thermal};
		return ups;
	}

	/**
	 * @brief Accumulate neutrino-cooling contribution into the thermal RHS.
	 *
	 * This method must **accumulate** into `dYdt` (never overwrite).
	 *
	 * ### Algorithm outline
	 * 1. Compute derived neutrino quantities via `Detail::ComputeDerived(*this, Y, ctx)`.
	 * 2. If the derived bundle is not valid (`ok=false`), skip contribution.
	 * 3. Otherwise add:
	 *    \f[
	 *      \frac{dx}{dt} \mathrel{+}= \left(\frac{d}{dt}\ln(T_\infty/T_{\rm ref})\right)_{\nu}
	 *    \f]
	 *    to the Thermal state component 0.
	 *
	 * @param t    Current simulation time (seconds by convention).
	 * @param Y    Composite state vector (read-only).
	 * @param dYdt RHS accumulator (write-only additive).
	 * @param ctx  Driver context (geometry, EOS, microphysics handles, etc.).
	 */
	void AccumulateRHS(double t,
					   const Evolution::StateVector &Y,
					   Evolution::RHSAccumulator &dYdt,
					   const Evolution::DriverContext &ctx) const override;

	// -------------------------------------------------------------------------
	// Options accessors
	// -------------------------------------------------------------------------

	/// Get current configuration options.
	[[nodiscard]] const Options &GetOptions() const { return opts_; }

	/// Replace configuration options.
	void SetOptions(const Options &o) { opts_ = o; }

	// -------------------------------------------------------------------------
	// IDriverDiagnostics interface
	// -------------------------------------------------------------------------

	/**
	 * @brief Producer name used in diagnostics packets and the catalog.
	 *
	 * This must remain stable and match the producer string used by observers.
	 */
	[[nodiscard]] std::string DiagnosticsName() const override;

	/**
	 * @brief Unit contract describing this driver's assumptions.
	 *
	 * Use this to document unit conventions and interpretability constraints.
	 */
	[[nodiscard]] Evolution::Diagnostics::UnitContract UnitContract() const override;

	/**
	 * @brief Diagnostics schema catalog for this driver.
	 *
	 * Declares which scalar keys this driver may emit, and defines profile(s)
	 * such as `timeseries_default` used by the TimeSeriesObserver.
	 */
	[[nodiscard]] Evolution::Diagnostics::ProducerCatalog DiagnosticsCatalog() const override;

	/**
	 * @brief Emit diagnostics scalars at a snapshot in time.
	 *
	 * This function should be deterministic. It delegates to `Detail::Diagnose()`
	 * so the emitted diagnostics stay consistent with the physics computation.
	 *
	 * @param t   Current time.
	 * @param Y   Current composite state (read-only).
	 * @param ctx Driver context (read-only).
	 * @param out Packet object to populate.
	 */
	void DiagnoseSnapshot(double t,
						  const Evolution::StateVector &Y,
						  const Evolution::DriverContext &ctx,
						  Evolution::Diagnostics::DiagnosticPacket &out) const override;

  private:
	/// Stored configuration.
	Options opts_{};

	/// Allow the Details bundle to access private members if needed.
	friend struct CompactStar::Physics::Driver::Thermal::Detail::NeutrinoCooling_Details;

	// friend NeutrinoCooling_Details Detail::ComputeDerived(const NeutrinoCooling &, ...);
	/**
	 * @brief Profile-versioned cache for neutrino luminosity coefficients.
	 *
	 * Marked mutable because AccumulateRHS/DiagnoseSnapshot are const by design
	 * but may rebuild cached coefficients when the StarProfile version changes.
	 */
	mutable NeutrinoCoolingProfileCache cache_;

	/**
	 * @brief Geometry provenance the cached payload was built against (ADR-0003).
	 *
	 * The payload genuinely depends on the GeometryCache (it integrates the geometry
	 * weights), but @c ProfileVersionedCache deliberately knows profile provenance ONLY —
	 * consumer-specific dependencies stay with the consumer. This driver therefore carries
	 * the geometry half of its own validity condition here.
	 */
	mutable Evolution::ProfileProvenance cached_geo_prov_{};

	/**
	 * @brief Access cached coefficients, rebuilding if the profile version changed.
	 *
	 * Implementation lives in NeutrinoCooling_Details.cpp to keep this header light.
	 */
	const NeutrinoCoolingCachePayload &Cache_(const Evolution::DriverContext &ctx) const;
};

} // namespace CompactStar::Physics::Driver::Thermal