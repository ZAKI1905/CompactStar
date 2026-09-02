// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2026
 * Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file AngularVelocity.hpp
 * @brief The public scientific spin input, and the single owner of the physical <-> geometric
 *        angular-velocity conversion.
 *
 * Governed by **ADR-0006** (ACCEPTED 2026-09-02), Q1 = A + D and P1/P2:
 *
 *   - the public scientific spin input is a **physical angular velocity in rad s^-1**;
 *   - it is carried by an **explicit typed quantity**, never a bare `double`, so that the unit is
 *     part of the type rather than of a comment;
 *   - the internal representation of every rotation *solution* stays geometric (km^-1, INV-02),
 *     and **exactly one** conversion turns one into the other.
 *
 * WHY A TYPE AT ALL. Before ADR-0006 the rotation entry points took a geometric `omega_bar`
 * *seed* in km^-1 that was not an angular velocity at all, while `HartleResult::Omega` was
 * annotated `[s^-1]` and stored km^-1, and the legacy exports advertised `(1/s)` for
 * seed-normalized values (`docs/validation/PHASE4_ROTATION_ENTRY.md` §9). Every one of those
 * defects is a `double` crossing a boundary with its unit recorded only in prose. A one-member
 * value type removes the whole class of defect at compile time and costs nothing at run time.
 *
 * THE CONVERSION. With `c` the speed of light in km s^-1,
 *
 *     Omega_geom [km^-1] = Omega_phys [rad s^-1] / c
 *
 * `Zaki::Physics::LIGHT_C_KM_S` is the authoritative constant, per the owner's recorded intent
 * that generic physical constants live in ZakiLib while domain conversions live in CompactStar
 * (`PHASE3_CLOSEOUT.md` §12, and the k_B precedent of Phase 3C). This header therefore
 * deliberately includes `<Zaki/Physics/Constants.hpp>` rather than re-deriving `c`: a second,
 * local `c` would recreate exactly the dual-authority defect Phase 3C spent an increment
 * eliminating. This is the one respect in which it differs from its top-level neighbours
 * `Units.hpp` (includes nothing) and `Geometry.hpp` (standard library only).
 *
 * LAYERING. Beyond that constant it includes only `<cmath>`, `<stdexcept>` and `<string>`. It
 * depends on **no** CompactStar component — not `Core`, not `EOS`, not `Physics`, not
 * `Extensions` — so it adds no edge to the internal dependency graph and may be included from
 * any layer. That is why it sits at the top level beside `Units.hpp` and `Geometry.hpp`: the
 * first-order consumer is `Core` today and the rotochemical consumer will be
 * `Physics/Evolution` tomorrow, and neither should have to reach through the other.
 *
 * NOT OWNED HERE, each for a stated reason:
 *   - The moment of inertia, angular momentum and the first-order profiles. Those are rotation
 *     *solutions*; they live with the solver that produces them.
 *   - Any cgs conversion of `J` or `I`. Those need `G`, whose authority
 *     (`Zaki::Physics::SUN_M_KM` vs `GSL_CONST_CGSM_SOLAR_MASS`, ~6.2e-5 apart) is still
 *     unadjudicated. ADR-0006 §10 defers them deliberately; `c` alone is enough for `Omega`.
 *   - The arbitrary internal `omega_bar` seed. It is not an angular velocity and has no
 *     physical meaning (ADR-0006 Q2); it never appears in this header.
 *
 * SIGN. Negative values are **valid**. The spin-down driver writes
 * `Omega_dot = -K sign(Omega) |Omega|^n` (`Physics/Driver/Spin/MagneticDipole.hpp:17,50`), which
 * is meaningful only if `Omega` may be negative. No positive-only convention is imposed here,
 * because no governing document establishes one.
 */

#ifndef CompactStar_AngularVelocity_H
#define CompactStar_AngularVelocity_H

#include <cmath>
#include <stdexcept>
#include <string>

#include <Zaki/Physics/Constants.hpp>

namespace CompactStar
{

namespace angular_velocity_detail
{
/// Formats a domain violation uniformly, so every fail-closed path is recognizable in a log.
/// Mirrors the convention of `CompactStar/Geometry.hpp`.
inline std::string DomainError(const char *what, double value)
{
	return std::string("CompactStar::AngularVelocity: non-representable ") + what + " (" +
		   std::to_string(value) +
		   "). Governed by ADR-0006 (physical spin input); this value has no physical reading "
		   "and is not silently regularized.";
}
} // namespace angular_velocity_detail

/**
 * @class AngularVelocity
 * @brief A physical angular velocity, stored in rad s^-1, convertible to geometric km^-1.
 *
 * Value semantics, one `double`, no state beyond it, no hierarchy, no unit framework. The
 * default-constructed value is **exactly zero**, which is a valid, well-defined non-rotating
 * spin (ADR-0006 P5).
 *
 * Construction is only through the named factories, so a bare number can never be mistaken for
 * a particular unit at a call site: `AngularVelocity::FromHz(716.0)` cannot be confused with
 * `AngularVelocity::FromRadPerSecond(716.0)`, and neither can be confused with a geometric
 * value — there is no factory that accepts km^-1 (ADR-0006 P1).
 */
class AngularVelocity
{
  public:
	/// Exactly zero spin. `Omega = 0` is valid and physically meaningful (ADR-0006 P5).
	AngularVelocity() = default;

	/// Physical angular velocity in rad s^-1 — the canonical scientific input (ADR-0006 Q1).
	/// @throws std::runtime_error if @p omega_rad_per_s is not finite.
	static AngularVelocity FromRadPerSecond(double omega_rad_per_s)
	{
		if (!std::isfinite(omega_rad_per_s))
			throw std::runtime_error(
				angular_velocity_detail::DomainError("angular velocity [rad/s]", omega_rad_per_s));
		return AngularVelocity(omega_rad_per_s);
	}

	/// Spin frequency in Hz, with `Omega = 2 pi f`. The factor `2 pi` is applied here and
	/// nowhere else; there is deliberately no implicit Hz <-> rad/s conversion anywhere.
	/// @throws std::runtime_error if @p frequency_hz is not finite.
	static AngularVelocity FromHz(double frequency_hz)
	{
		if (!std::isfinite(frequency_hz))
			throw std::runtime_error(
				angular_velocity_detail::DomainError("spin frequency [Hz]", frequency_hz));
		return AngularVelocity(2.0 * M_PI * frequency_hz);
	}

	/// Spin period in seconds, with `Omega = 2 pi / P`.
	/// @throws std::runtime_error if @p period_s is not finite or is zero — unlike a zero
	///         *frequency*, a zero *period* is not a slow limit but an undefined quantity.
	static AngularVelocity FromPeriodSeconds(double period_s)
	{
		if (!std::isfinite(period_s))
			throw std::runtime_error(
				angular_velocity_detail::DomainError("spin period [s]", period_s));
		if (period_s == 0.0)
			throw std::runtime_error(angular_velocity_detail::DomainError(
				"spin period [s] (a zero period is not a limit of a physical spin; pass a zero "
				"angular velocity or frequency instead)",
				period_s));
		return AngularVelocity(2.0 * M_PI / period_s);
	}

	/// The physical angular velocity in rad s^-1.
	[[nodiscard]] double RadPerSecond() const noexcept { return omega_rad_per_s_; }

	/// The geometric angular velocity in km^-1: `Omega_geom = Omega_phys / c`.
	///
	/// **This is the single authoritative physical -> geometric conversion for angular
	/// quantities on the governed first-order rotation path** (ADR-0006 P2). Nothing else may
	/// divide or multiply an angular velocity by the speed of light.
	[[nodiscard]] double GeomKmInverse() const noexcept
	{
		return omega_rad_per_s_ / Zaki::Physics::LIGHT_C_KM_S;
	}

	/// True for exactly zero spin. Provided so callers need not compare a `double` to `0.0`.
	[[nodiscard]] bool IsZero() const noexcept { return omega_rad_per_s_ == 0.0; }

  private:
	explicit AngularVelocity(double omega_rad_per_s) noexcept
		: omega_rad_per_s_(omega_rad_per_s)
	{
	}

	/// The one stored quantity. Physical, because that is what callers supply and what
	/// `SpinState::Omega()` evolves; the geometric form is derived on demand so that no
	/// duplicate representation can drift (ADR-0006 Q3).
	double omega_rad_per_s_ = 0.0;
};

/**
 * @brief The single authoritative **geometric -> physical** angular conversion, `Omega_phys =
 *        Omega_geom * c`. The exact inverse of `AngularVelocity::GeomKmInverse()`.
 *
 * Together with that member function this is the *only* place on the governed first-order
 * rotation path where an angular velocity meets the speed of light (ADR-0006 P2).
 *
 * It deliberately returns a plain `double` in rad s^-1 rather than an `AngularVelocity`. This
 * is a **rendering** of an already-computed geometric result for display or for a physical
 * accessor — it is *not* a spin input, and ADR-0006 P1 forbids requesting a spin in km^-1. By
 * not producing an `AngularVelocity`, the function cannot be used to smuggle a geometric number
 * back into the scientific input path.
 */
[[nodiscard]] inline double
AngularVelocityGeomToRadPerSecond(double omega_geom_km_inverse) noexcept
{
	return omega_geom_km_inverse * Zaki::Physics::LIGHT_C_KM_S;
}

} // namespace CompactStar

#endif /* CompactStar_AngularVelocity_H */
