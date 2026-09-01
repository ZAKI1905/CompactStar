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
 * @file ProfileProvenance.hpp
 * @brief Runtime provenance of a quantity derived from a @c StarProfile.
 *
 * Governed by **ADR-0003** (ACCEPTED). The owner adjudicated:
 *   - **Q2 = Option A** — provenance is `(source StarProfile object identity, Version())`.
 *
 * @par Lifetime contract (binding)
 * A profile-derived context, cache or snapshot **MUST NOT outlive the @c StarProfile from
 * which it was built.** Pointer identity is therefore meaningful for the whole lifetime of
 * every object that carries it. This is a contract, not a mechanism: if it is violated the
 * derived object is already invalid for reasons no token could detect.
 *
 * @par Deliberately NOT provided
 * No generated UUID, no global counter, no serialized or persistent profile ID, and no
 * copy/move identity semantics. Persistent or cross-process profile identity is a separate
 * architectural requirement and would need its own decision (ADR-0003 §6, Option B).
 *
 * @par Why identity and not version alone
 * @c StarProfile::Version() is the *revision counter of one particular profile*, not a global
 * star identity: two independently constructed profiles routinely both report `1`. Keying a
 * reusable cache on the version alone therefore collides across stars — measured at 85.7 %
 * relative error before this contract existed (`docs/validation/CACHE_CORRECTNESS.md`).
 *
 * This header includes only @c <cstdint> and forward-declares @c StarProfile, so it adds no
 * weight to the headers that carry provenance.
 */

#ifndef CompactStar_Physics_Evolution_ProfileProvenance_H
#define CompactStar_Physics_Evolution_ProfileProvenance_H

#include <cstdint>

namespace CompactStar::Core
{
struct StarProfile;
}

namespace CompactStar::Physics::Evolution
{

/**
 * @struct ProfileProvenance
 * @brief Identity + revision of the @c StarProfile a derived quantity was built from.
 *
 * Comparison is exact: same profile object **and** same revision. Two provenance tokens that
 * compare equal denote the same profile in the same state, so any two quantities derived
 * deterministically from that state are interchangeable — which is why a consumer keys on this
 * token rather than on the address of the derived object itself (ADR-0003 §11).
 */
struct ProfileProvenance
{
	/// Source profile. Non-owning; valid for the lifetime of the source (see lifetime contract).
	const CompactStar::Core::StarProfile *source = nullptr;

	/// @c StarProfile::Version() observed when the derived quantity was built.
	std::uint64_t version = 0;

	/// True once bound to a profile.
	[[nodiscard]] bool IsSet() const noexcept { return source != nullptr; }

	friend bool operator==(const ProfileProvenance &a, const ProfileProvenance &b) noexcept
	{
		return a.source == b.source && a.version == b.version;
	}
	friend bool operator!=(const ProfileProvenance &a, const ProfileProvenance &b) noexcept
	{
		return !(a == b);
	}
};

} // namespace CompactStar::Physics::Evolution

#endif /* CompactStar_Physics_Evolution_ProfileProvenance_H */
