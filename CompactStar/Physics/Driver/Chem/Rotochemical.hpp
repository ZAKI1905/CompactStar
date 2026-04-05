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
 * @file Rotochemical.hpp
 * @brief Driver for rotochemical heating following Fernandez & Reisenegger (2005).
 *
 * As the star spins down, the centrifugal flattening decreases and the
 * equilibrium composition at each pressure shell shifts. This drives
 * chemical imbalances (delta_mu) that source heating and weak reactions.
 *
 * The rate of change of the chemical imbalance for species i is:
 *
 *   d(eta_i)/dt += Z_i * 2 * Omega * dOmega/dt
 *
 * where Z_i = dN_i/dOmega^2 are the reduced rotochemical coefficients
 * precomputed in RotochemicalCache.
 *
 * Directory-based namespace: CompactStar::Physics::Driver::Chem
 *
 * @ingroup PhysicsDriver
 */

#ifndef CompactStar_Physics_Driver_Chem_Rotochemical_H
#define CompactStar_Physics_Driver_Chem_Rotochemical_H

#include <string>
#include <vector>

#include "CompactStar/Physics/Driver/IDriver.hpp"
#include "CompactStar/Physics/Evolution/RotochemicalCache.hpp"
#include "CompactStar/Physics/State/Tags.hpp"

namespace CompactStar::Physics::Driver::Chem
{

/**
 * @class Rotochemical
 * @brief Evolution driver for rotochemical heating.
 *
 * **Depends on:** Spin (reads Omega, dOmega/dt)
 * **Updates:**    Chem (chemical imbalances)
 *
 * The driver uses precomputed Z_i coefficients from RotochemicalCache.
 * These coefficients encode the full GR (Hartle O(Omega^2)) structural
 * response plus the equilibrium-sequence correction with baryon conservation.
 */
class Rotochemical final : public IDriver
{
  public:
	Rotochemical() = default;

	/// Construct with a precomputed rotochemical cache.
	explicit Rotochemical(const Evolution::RotochemicalCache &cache)
		: cache_(cache)
	{
	}

	// ------------------------------------------------------------------
	//  IDriver interface
	// ------------------------------------------------------------------

	std::string Name() const override { return "Rotochemical"; }

	const std::vector<State::StateTag> &DependsOn() const override
	{
		static const std::vector<State::StateTag> deps{
			State::StateTag::Spin};
		return deps;
	}

	const std::vector<State::StateTag> &Updates() const override
	{
		static const std::vector<State::StateTag> ups{
			State::StateTag::Chem};
		return ups;
	}

	/**
	 * @brief Add rotochemical contribution to dY/dt.
	 *
	 * For each species i with a nonzero Z_i:
	 *   d(eta_i)/dt += Z_i * 2 * Omega * Omega_dot
	 *
	 * where:
	 *   - Omega is read from SpinState
	 *   - Omega_dot is the current spin-down rate (from the Spin block of dY/dt,
	 *     or from a stored value in the DriverContext)
	 *   - Z_i = dN_i/dOmega^2 from the RotochemicalCache
	 *
	 * The mapping from species label to ChemState index depends on
	 * the ChemState layout configured in the EvolutionSystem.
	 */
	void AccumulateRHS(double t,
					   const Evolution::StateVector &Y,
					   Evolution::RHSAccumulator &dYdt,
					   const Evolution::DriverContext &ctx) const override;

	// ------------------------------------------------------------------
	//  Cache access
	// ------------------------------------------------------------------

	const Evolution::RotochemicalCache &GetCache() const { return cache_; }
	void SetCache(const Evolution::RotochemicalCache &c) { cache_ = c; }

  private:
	Evolution::RotochemicalCache cache_;
};

} // namespace CompactStar::Physics::Driver::Chem

#endif /* CompactStar_Physics_Driver_Chem_Rotochemical_H */
