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
 * @file NeutrinoCooling_Cache.hpp
 * @brief Version-keyed cache payload and wiring for NeutrinoCooling.
 *
 * This header defines the *payload* cached against a StarProfile version
 * and is intentionally small so NeutrinoCooling can own the cache as a
 * normal member without pulling in heavy microphysics headers.
 *
 * @ingroup PhysicsDriver
 */

#include <cstddef>
#include <cstdint>

#include "CompactStar/Physics/Evolution/ProfileCache.hpp" // ProfileVersionedCache<>

/* Forward declarations only; keep this header light. */
namespace CompactStar::Physics::Evolution
{
class StarContext;
} // namespace CompactStar::Physics::Evolution

namespace CompactStar::Physics::Driver::Thermal
{

/**
 * @struct NeutrinoCoolingCachePayload
 * @brief Cached geometric/matter coefficients for fast neutrino luminosity evaluation.
 *
 * ## Intent
 * In the isothermal interior, local temperature is
 * \f[
 *   T(r) = T_\infty\, e^{-\nu(r)}.
 * \f]
 * Many emissivity channels can be written schematically as:
 * \f[
 *   Q_\nu(r) = \mathcal{A}(r)\,T(r)^n,
 * \f]
 * so the redshifted luminosity becomes
 * \f[
 *   L_{\nu,\infty} = \int dV\; e^{2\nu(r)}\, \mathcal{A}(r)\,[T_\infty e^{-\nu(r)}]^n
 *                = T_\infty^n \int dV\; \mathcal{A}(r)\, e^{(2-n)\nu(r)}.
 * \f]
 *
 * Therefore, for a fixed star profile (fixed geometry/composition), one can cache:
 * \f[
 *   K_n \equiv \int dV\; \mathcal{A}(r)\, e^{(2-n)\nu(r)},
 * \f]
 * and evaluate \f$L_{\nu,\infty} = K_n\,T_\infty^n\f$ extremely quickly.
 *
 * ## What is cached here
 * This payload stores channel-specific coefficients (e.g., DUrca uses n=6,
 * MUrca uses n=8) in cgs luminosity units such that multiplying by
 * \f$T_\infty^n\f$ yields \f$L_{\nu,\infty}\f$ in [erg/s].
 *
 * ## Versioning
 * The cache is rebuilt automatically when the underlying StarProfile version changes.
 */
struct NeutrinoCoolingCachePayload
{
	/// True if the payload has been built successfully at least once.
	bool valid = false;

	/// StarProfile::Version() this payload was built against.
	std::uint64_t built_version = 0;

	/// Number of radial zones in the source profile at build time.
	std::size_t n_zones = 0;

	/// Cached DUrca boundary index (last allowed), copied from StarContext for convenience.
	/// -1 means DUrca nowhere allowed or unavailable.
	long durca_last_allowed = -1;

	// ---------------------------------------------------------------------
	// Cached coefficients for fast evaluation (units chosen by builder).
	// ---------------------------------------------------------------------

	/**
	 * @brief DUrca coefficient for \f$L_{\nu,\infty}^{DU} = K_{DU}\,T_\infty^{6}\f$.
	 *
	 * Units must be [erg/s/K^6].
	 */
	double K_DU_erg_s_K6 = 0.0;

	/**
	 * @brief MUrca coefficient for \f$L_{\nu,\infty}^{MU} = K_{MU}\,T_\infty^{8}\f$.
	 *
	 * Units must be [erg/s/K^8].
	 */
	double K_MU_erg_s_K8 = 0.0;

	/**
	 * @brief PBF coefficient for \f$L_{\nu,\infty}^{PBF} = K_{PBF}\,T_\infty^{n}\f$ (future hook).
	 *
	 * Units depend on chosen exponent. For now keep placeholder = 0.
	 */
	double K_PBF = 0.0;
};

/**
 * @brief Cache type used by NeutrinoCooling (payload keyed by StarProfile version).
 */
using NeutrinoCoolingProfileCache =
	CompactStar::Physics::Evolution::ProfileVersionedCache<NeutrinoCoolingCachePayload>;

} // namespace CompactStar::Physics::Driver::Thermal