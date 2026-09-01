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
#ifndef CompactStar_Physics_Evolution_ProfileCache_H
#define CompactStar_Physics_Evolution_ProfileCache_H

#include <cstdint>
#include <utility>

#include "CompactStar/Physics/Evolution/ProfileProvenance.hpp"

namespace CompactStar::Physics::Evolution
{

class StarContext;

/**
 * @brief Return the StarProfile version currently bound to a StarContext.
 *
 * This function is defined out-of-line in ProfileCache.cpp so that ProfileCache.hpp
 * does not need to include StarContext.hpp (keeps headers light and avoids
 * dependency fan-out).
 *
 * @param sc StarContext bound to a StarProfile.
 * @return Monotonic StarProfile version (0 if context is invalid by policy).
 */
std::uint64_t ProfileVersion(const StarContext &sc);

/**
 * @brief Return the full provenance — (profile identity, version) — bound to @p sc.
 *
 * ADR-0003. Defined out-of-line for the same header-weight reason as ProfileVersion().
 */
ProfileProvenance ProfileProvenanceOf(const StarContext &sc);

/**
 * @file ProfileCache.hpp
 * @brief Lightweight, reusable caching utility keyed by @c StarProfile::Version().
 *
 * This header defines a small generic cache object that stores a user-defined
 * payload (e.g., precomputed integrals, masks, boundary indices, lookup tables)
 * and rebuilds it automatically when the underlying star profile changes.
 *
 * ## Design intent
 * Many evolution drivers (thermal, chemical heating, neutrino cooling, etc.)
 * want to cache derived quantities computed from a fixed equilibrium star
 * profile (radial columns and their derived caches). Your core system already
 * provides a mutation/version counter through @c StarProfile::Version() and
 * @c StarContext as the read-only adapter to that profile.
 *
 * This utility:
 * - tracks the last-seen profile version,
 * - rebuilds the payload exactly once per new version,
 * - avoids repeating the "check version / invalidate / rebuild" plumbing in
 *   every consumer class.
 *
 * ## Required contract
 * The cache assumes that @c StarContext can provide the current profile version.
 * The recommended interface is a public method:
 * @code
 * std::uint64_t StarContext::ProfileVersion() const;
 * @endcode
 *
 * (This can simply forward to the bound @c StarProfile::Version().)
 *
 * ## Thread-safety
 * This cache is *not* thread-safe. It uses mutable state and is intended for
 * the single-threaded RHS/diagnostics flow typical in your evolution loop.
 * If you later run per-zone physics in parallel, do not share a single instance
 * of this cache across threads without external synchronization.
 */

/**
 * @class ProfileVersionedCache
 * @brief Generic cache keyed by the bound @c StarProfile::Version().
 *
 * @tparam Payload Arbitrary user-defined payload type that stores cached data.
 *
 * The cache owns a @p Payload instance and guarantees:
 * - The payload is (re)built when either:
 *   - the cache has never been built, or
 *   - the current star-profile version differs from the last-seen version.
 * - Otherwise, the stored payload is returned without calling the builder.
 *
 * ## Builder contract
 * The caller supplies a builder functor/lambda with signature:
 * @code
 * void builder(const StarContext& sc, Payload& out);
 * @endcode
 * The builder must fully populate @p out for the current profile state.
 *
 * ## Typical usage
 * @code
 * struct MyPayload { double A = 0.0; long idx = -1; };
 * mutable ProfileVersionedCache<MyPayload> cache;
 *
 * const MyPayload& p = cache.Get(star_ctx, [&](const StarContext& sc, MyPayload& out){
 *     out.A = ComputeSomething(sc);
 *     out.idx = sc.DirectUrcaLastAllowedIndex();
 * });
 * // p.A and p.idx are now valid and cached against sc.ProfileVersion()
 * @endcode
 *
 * ## Notes on invalidation
 * Invalidation is version-driven, not pointer-driven. If the profile is mutated,
 * it should call @c StarProfile::Touch() (or use @c StarProfile::EditScope),
 * which increments @c StarProfile::Version() and triggers a rebuild on the next
 * access.
 */
template <typename Payload>
class ProfileVersionedCache
{
  public:
	/**
	 * @brief Construct an empty cache in the "not built" state.
	 *
	 * The first call to Get/GetMutable will invoke the builder and store the
	 * resulting payload.
	 */
	ProfileVersionedCache() = default;

	/**
	 * @brief Force a rebuild on the next access.
	 *
	 * This explicitly invalidates the cache regardless of the star-profile
	 * version. The next call to Get/GetMutable will invoke the builder.
	 *
	 * This is useful if a consumer wants to clear internal state even when
	 * the star profile hasn't changed (e.g., toggling options that affect what
	 * is cached, or resetting between runs).
	 */
	void Invalidate() noexcept
	{
		m_built = false;
		m_prov = ProfileProvenance{};
	}

	/**
	 * @brief Access the cached payload (rebuilding if needed).
	 *
	 * @tparam Builder Callable type used to build the payload.
	 * @param sc StarContext bound to a StarProfile that provides a monotonic version.
	 * @param builder Callable invoked as @c builder(sc, payload) when a rebuild occurs.
	 * @return Const reference to the cached payload.
	 *
	 * The builder is invoked if:
	 * - the cache has never been built, or
	 * - @c sc.ProfileVersion() differs from the last version used to build.
	 */
	template <typename Builder>
	const Payload &Get(const StarContext &sc, Builder &&builder) const
	{
		// ADR-0003: the key is (profile identity, version), NOT the version alone. Two
		// independently built profiles routinely share a numeric version, so a
		// version-only key silently served one star's payload for another.
		const ProfileProvenance prov = ProfileProvenanceOf(sc);
		if (!m_built || prov != m_prov)
		{
			builder(sc, m_payload); // customization point: define what is cached
			m_prov = prov;
			m_built = true;
		}
		return m_payload;
	}

	/**
	 * @brief Access the cached payload with mutable reference (rebuilding if needed).
	 *
	 * @tparam Builder Callable type used to build the payload.
	 * @param sc StarContext bound to a StarProfile that provides a monotonic version.
	 * @param builder Callable invoked as @c builder(sc, payload) when a rebuild occurs.
	 * @return Mutable reference to the cached payload.
	 *
	 * This is occasionally useful when the consumer wants to perform small
	 * post-build adjustments or attach extra metadata after the versioned build.
	 *
	 * Note: the cache still keys rebuild decisions on profile version; any
	 * post-build mutations performed by the caller do not change that version.
	 */
	template <typename Builder>
	Payload &GetMutable(const StarContext &sc, Builder &&builder) const
	{
		(void)Get(sc, std::forward<Builder>(builder));
		return m_payload;
	}

	/**
	 * @brief Whether the cache has been built at least once.
	 *
	 * @return True if a payload has been constructed and stored; false otherwise.
	 */
	bool IsBuilt() const noexcept { return m_built; }

	/**
	 * @brief The star-profile version used when the cache was last built.
	 *
	 * @return Last-seen profile version (0 if never built or explicitly invalidated).
	 */
	std::uint64_t BuiltAgainstVersion() const noexcept { return m_prov.version; }

	/**
	 * @brief Full provenance — identity and version — the payload was built against.
	 *
	 * ADR-0003. Identity is what distinguishes two profiles that happen to share a
	 * numeric version.
	 */
	const ProfileProvenance &BuiltAgainst() const noexcept { return m_prov; }

  private:
	// /**
	//  * @brief Obtain the current profile version from the StarContext.
	//  *
	//  * This indirection keeps the cache independent of how StarContext stores
	//  * or exposes the bound profile. The recommended implementation is a public
	//  * @c StarContext::ProfileVersion() accessor.
	//  *
	//  * @param sc StarContext to query.
	//  * @return Current profile version.
	//  */
	// static std::uint64_t ProfileVersion(const StarContext &sc)
	// {
	// 	return sc.ProfileVersion();
	// }

  private:
	/// True if @c m_payload contains valid data built by the provided builder.
	mutable bool m_built = false;

	/// Provenance (profile identity + version) @c m_payload was built against.
	mutable ProfileProvenance m_prov{};

	/// Cached payload owned by the cache object.
	mutable Payload m_payload{};
};

} // namespace CompactStar::Physics::Evolution

#endif /* CompactStar_Physics_Evolution_ProfileCache_H */