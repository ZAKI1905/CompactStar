// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025 Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file Tags.hpp
 * @brief Enumeration of physical state components used by evolution drivers.
 *
 * Provides canonical identifiers (`StateTag`) for each sub-state (spin, thermal,
 * chemical, BNV, etc.) so that drivers and the evolution graph can declare
 * dependencies and update targets in a type-safe manner.
 *
 * @ingroup PhysicsState
 *
 * Example:
 * @code
 * static constexpr std::array<StateTag,2> deps{
 *     StateTag::Spin,
 *     StateTag::Thermal
 * };
 * @endcode
 *
 * @see CompactStar::Physics::Driver::IDriver
 * @see CompactStar::Physics::Evolution::System
 *
 * @author
 *   Mohammadreza Zakeri (M.Zakeri@eku.edu)
 */

#ifndef CompactStar_Physics_State_Tags_H
#define CompactStar_Physics_State_Tags_H

#include <string>
#include <string_view>

namespace CompactStar::Physics::State
{

/**
 * @enum StateTag
 * @brief Enumerates high-level dynamic subsystems ("state blocks") of the model.
 *
 * Each distinct `StateTag` corresponds to one concrete state block
 * (SpinState, ThermalState, ChemState, BNVState, …) in Physics/State,
 * which the evolution System packs into the global StateVector.
 */
enum class StateTag : unsigned int
{
	Spin = 0,  ///< Rotational degree of freedom (Ω, Ω̇, I, etc.)
	Thermal,   ///< Temperature / energy content (T̃, C_v, L_ν, L_γ)
	Chem,	   ///< Chemical imbalance (η_e, η_μ, etc.)
	BNV,	   ///< Baryon-number-violating or exotic particle population
	Structure, ///< Optional: structural/geometric parameters (R, M, etc.)
	Magnetic,  ///< Optional: magnetic-field state (if evolved explicitly)
	Custom	   ///< Placeholder for user-extended or experimental subsystems
};

/**
 * @brief Convert a `StateTag` to a short human-readable name.
 *
 * @param tag  Enumeration value.
 * @return Constant string_view naming the state.
 */
inline constexpr std::string_view ToString(StateTag tag)
{
	switch (tag)
	{
	case StateTag::Spin:
		return "Spin";
	case StateTag::Thermal:
		return "Thermal";
	case StateTag::Chem:
		return "Chem";
	case StateTag::BNV:
		return "BNV";
	case StateTag::Structure:
		return "Structure";
	case StateTag::Magnetic:
		return "Magnetic";
	case StateTag::Custom:
		return "Custom";
	}
	return "Unknown";
}

/**
 * @brief Return a printable `std::string` from a `StateTag`.
 *
 * @param tag Enumeration value.
 * @return Printable string.
 */
inline std::string ToStringCopy(StateTag tag)
{
	return std::string(ToString(tag));
}

} // namespace CompactStar::Physics::State

#endif /* CompactStar_Physics_State_Tags_H */