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
 * @file Geometry.hpp
 * @brief Single mathematical owner for the static-metric radial factor and the proper-volume
 *        measure of a spherically symmetric star.
 *
 * Governed by **ADR-0004** (ACCEPTED 2026-09-01), which resolves INV-04's "single owner"
 * requirement into three distinct roles:
 *
 *   A. **Mathematical owner — this header.** Owns `f`, `Λ`, `e^Λ`, `w_V` *and the domain and
 *      failure semantics*. Knows nothing about profiles, caches, species or units beyond km.
 *   B. **Cached-representation owner — `Physics::Evolution::GeometryCache`.** Owns the
 *      precomputed `ExpLambda`/`WV`/`WVExpNu`/`WVExp2Nu` arrays and their ADR-0003 provenance.
 *      It obtains the *formula* from here; it keeps its own `DataColumn` composition.
 *   C. **Consumer integrands — `NStar`, the thermal drivers, ….** Own their own physics factor
 *      (`n_B`, `c_V`, `Q_ν`) and their own unit conversions. They must never re-derive `w_V`.
 *
 * DEFINITIONS (ADR-0004 §5; INV-03 for the metric convention). With `r` and `m` both in
 * geometric length units (km, `m[km] = GM/c²`, per `CURRENT_ARCHITECTURE.md` §1):
 *
 *     f(r,m)   = 1 − 2m/r                  [metric denominator, dimensionless]
 *     Λ(r,m)   = −½ ln f                   [g_rr = e^{2Λ}]
 *     e^Λ      = exp(Λ)                    [derived from the canonical Λ, NOT from f^{-1/2}]
 *     w_V(r,m) = 4π r² e^Λ                 [proper-volume weight, km² per km of coordinate r]
 *
 * `e^Λ` is deliberately `exp(−½ ln f)` rather than `1/sqrt(f)`. The two agree to ≤ 2 ULP
 * (ADR-0004 §6.2) but are not bitwise equal at high compactness, and the logarithmic form is
 * what every authenticated profile and cache already stores.
 *
 * DOMAIN AND FAILURE SEMANTICS (ADR-0004 §0-Q3 — the *hybrid physical-domain contract*).
 * Before this contract the repository held **six mutually inconsistent behaviors** for the same
 * degenerate input, disagreeing by a factor of 3.16e7 at `f = 0` and including a silent
 * divisor substitution inherited from `Zaki::Vector::DataColumn::operator/=` (ADR-0004 §11).
 * This header replaces them, on the canonical path, with one rule:
 *
 *     r == 0 and m == 0   →  regular-centre limit: f = 1, Λ = 0, e^Λ = 1, w_V = 0
 *     r <  0              →  throw
 *     r == 0 and m != 0   →  throw
 *     r >  0 and f <= 0   →  throw   (horizon / outside the static stellar metric's domain)
 *     r >  0 and f >  0   →  evaluate, with NO clamp and NO epsilon
 *     r or m non-finite   →  throw
 *
 * There is no tolerance band around `r = 0`: the regular-centre case is the exact `r == 0,
 * m == 0` case. The former `denom <= 0 → 1e-15` clamp is **not** reproduced here; it was
 * unreachable on every validated path (measured `max 2m/r = 0.481` across the four
 * authenticated stars, ADR-0004 §11) and it silently manufactured a finite, physically
 * meaningless weight from invalid input.
 *
 * Failures are `std::runtime_error`, the project's ordinary convention. They are never
 * `assert`, never a log-and-continue, and never a process abort.
 *
 * LAYERING. This header includes only `<cmath>`, `<stdexcept>` and `<string>` — standard
 * library only. No `Core`, no `EOS`, no `Physics`, no `Zaki`, no Confind. It therefore adds no
 * edge to the dependency graph and may be included from any layer, which is precisely why it
 * sits at the top level beside `CompactStar/Units.hpp` rather than inside a subsystem. It holds
 * no state, no registry, no class hierarchy and no templates: four free scalar functions.
 *
 * NOT OWNED HERE, each for a stated reason:
 *   - The redshift factors `e^{ν}`, `e^{2ν}` and the composed measures `w_V e^{ν}`,
 *     `w_V e^{2ν}`. This primitive does not know about `ν`. `GeometryCache` composes them from
 *     the single `WV` array (ADR-0004 §12) so that no second measure can come into existence.
 *   - The baryon-number unit conversion `FM3_TO_KM3 = 1e54`. That belongs to INV-14 and to the
 *     consumer, not to `dV`. Putting it here would silently convert the units of every other
 *     consumer of `w_V` (ADR-0004 §14).
 *   - Any per-species factor `Y_i`. ADR-0001 and INV-14 govern that.
 *   - The Hartle ODE coefficients, the TOV mass ODE, the `EvaluateNu` boundary condition and
 *     `SurfaceGravity`. These share the token `4πr²` but are not this measure (ADR-0004 §4.4).
 */

#ifndef CompactStar_Geometry_H
#define CompactStar_Geometry_H

#include <cmath>
#include <stdexcept>
#include <string>

namespace CompactStar::Geometry
{

namespace detail
{
/// Formats a domain violation uniformly, so every fail-closed path is recognizable in a log.
inline std::string DomainError(const char *what, double r_km, double m_km)
{
	return std::string("CompactStar::Geometry: ") + what + " (r = " + std::to_string(r_km) +
		   " km, m = " + std::to_string(m_km) +
		   " km). Governed by ADR-0004 §0-Q3 (hybrid physical-domain contract); this input is "
		   "outside the domain of the static stellar metric and is not regularized.";
}

/// Validates `(r, m)` and reports whether this is the exact regular-centre case.
/// Throws on every input the accepted contract declares invalid.
inline bool ValidateAndIsCentre(double r_km, double m_km)
{
	if (!std::isfinite(r_km) || !std::isfinite(m_km))
		throw std::runtime_error(detail::DomainError("non-finite radius or mass", r_km, m_km));

	if (r_km < 0.0)
		throw std::runtime_error(detail::DomainError("negative radius", r_km, m_km));

	if (r_km == 0.0)
	{
		// The regular centre is the EXACT (r, m) = (0, 0) point. A non-zero enclosed mass at
		// zero radius is a point singularity, not a stellar centre.
		if (m_km != 0.0)
			throw std::runtime_error(
				detail::DomainError("zero radius with non-zero enclosed mass", r_km, m_km));
		return true;
	}

	return false;
}
} // namespace detail

/// Metric denominator `f = 1 − 2m/r`. Returns exactly 1 at the regular centre.
/// @throws std::runtime_error on any input outside the accepted domain (ADR-0004 §0-Q3).
inline double MetricDenominator(double r_km, double m_km)
{
	if (detail::ValidateAndIsCentre(r_km, m_km))
		return 1.0;

	// Expression order is load-bearing: it reproduces `NStar.cpp` and
	// `GeometryCache::DeriveLambdaFromMR_` bit-for-bit on the ordinary domain.
	const double f = 1.0 - 2.0 * m_km / r_km;

	if (f <= 0.0)
		throw std::runtime_error(detail::DomainError(
			"metric denominator 1 - 2m/r is not positive (at or inside the horizon)", r_km,
			m_km));

	return f;
}

/// Metric exponent `Λ = −½ ln f`. Exactly 0 at the regular centre.
/// @throws std::runtime_error on any input outside the accepted domain.
inline double Lambda(double r_km, double m_km)
{
	const double f = MetricDenominator(r_km, m_km);

	// `f == 1.0` at the centre gives std::log(1.0) == 0.0 exactly, so the centre needs no
	// special case here; it is kept implicit rather than branched to avoid a second definition.
	return -0.5 * std::log(f);
}

/// `e^Λ = exp(−½ ln f)`, the radial metric factor `sqrt(g_rr)`. Exactly 1 at the regular centre.
/// Derived from the canonical `Lambda`, deliberately not from `f^{-1/2}` (see file header).
/// @throws std::runtime_error on any input outside the accepted domain.
inline double ExpLambda(double r_km, double m_km)
{
	return std::exp(Lambda(r_km, m_km));
}

/// Proper-volume weight `w_V = 4π r² e^Λ`, in km² per km of coordinate radius, so that
/// `dV = w_V dr`. Exactly 0 at the regular centre (`r² = 0`).
///
/// This is the measure INV-04 governs. It carries **no** baryon density, **no** species
/// fraction and **no** unit conversion; the consumer supplies those.
/// @throws std::runtime_error on any input outside the accepted domain.
inline double ProperVolumeWeight(double r_km, double m_km)
{
	const double el = ExpLambda(r_km, m_km);

	// Associativity matches `GeometryCache::Build_`: (4π · r²) · e^Λ.
	return (4.0 * M_PI) * (r_km * r_km) * el;
}

} // namespace CompactStar::Geometry

#endif /* CompactStar_Geometry_H */
