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
 * @file ProfileCache.cpp
 * @brief Out-of-line helpers for ProfileCache to avoid heavy includes in headers.
 */

#include "CompactStar/Physics/Evolution/ProfileCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"

namespace CompactStar::Physics::Evolution
{

std::uint64_t ProfileVersion(const StarContext &sc)
{
	// Policy choice:
	// - If StarContext can be invalid, decide whether to return 0 or throw.
	// - Returning 0 is convenient and keeps caching logic simple; it will still
	//   rebuild if later the context binds to a real profile version != 0.
	//
	// If we prefer fail-fast, replace with:
	//   if (!sc.IsValid()) throw std::runtime_error(...);
	//
	return sc.ProfileVersion();
}

} // namespace CompactStar::Physics::Evolution