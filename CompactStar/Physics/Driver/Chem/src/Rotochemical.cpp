// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025 Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

// -----------------------------------------------------------------------------
//  Rotochemical.cpp
// -----------------------------------------------------------------------------
//  Implementation of the rotochemical heating driver.
//
//  As the star spins down, the centrifugal deformation decreases and the
//  equilibrium composition at each pressure shell shifts. This drives
//  chemical imbalances (delta_mu) that source heating and weak reactions.
//
//  The rate of change of the chemical imbalance for reaction channel i is:
//
//      d(eta_i)/dt += Z_i * 2 * Omega * dOmega/dt
//
//  where Z_i = dN_i/dOmega^2 are precomputed in RotochemicalCache,
//  combining Hartle O(Omega^2) structural perturbations with
//  equilibrium-sequence corrections under baryon conservation.
//
//  Ref: Fernandez & Reisenegger (2005), Reisenegger (1995)
// -----------------------------------------------------------------------------

#include "CompactStar/Physics/Driver/Chem/Rotochemical.hpp"
#include "CompactStar/Physics/Evolution/RHSAccumulator.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"
#include "CompactStar/Physics/State/ChemState.hpp"
#include "CompactStar/Physics/State/SpinState.hpp"

#include <cmath>

#include <Zaki/Util/Instrumentor.hpp>

namespace CompactStar::Physics::Driver::Chem
{

// -----------------------------------------------------------------------------
//  Rotochemical::AccumulateRHS
// -----------------------------------------------------------------------------
void Rotochemical::AccumulateRHS(double t,
								 const Evolution::StateVector &Y,
								 Evolution::RHSAccumulator &dYdt,
								 const Evolution::DriverContext &ctx) const
{
	PROFILE_FUNCTION();

	(void)t;
	(void)ctx;

	if (!cache_.valid)
		return;

	// -------------------------------------------------------------------------
	//  1. Read Omega from SpinState
	// -------------------------------------------------------------------------
	const auto &spin = Y.GetSpin();
	if (spin.NumComponents() == 0)
		return;

	const double Omega = spin.Omega();
	if (std::abs(Omega) == 0.0)
		return;

	// -------------------------------------------------------------------------
	//  2. Read dOmega/dt from the Spin block of dYdt
	// -------------------------------------------------------------------------
	// The spin-down rate Omega_dot is computed by the MagneticDipole driver
	// (or similar) and accumulated in the Spin block before this driver runs.
	// We read it from the accumulator to get the current spin-down rate.
	//
	// NOTE: This assumes driver ordering: Spin drivers run before Chem drivers.
	// If that's not guaranteed, we may need to store Omega_dot in the
	// DriverContext or compute it here from the SpinState.
	//
	// For now, use a simplified approach: read K and n from context or
	// compute Omega_dot = -K * Omega^n directly.
	// TODO: Wire proper Omega_dot access once the Evolution system
	// guarantees driver ordering.
	//
	// Placeholder: use the accumulated Spin RHS if available.
	// Read the accumulated Spin RHS (Omega_dot) from the Spin block.
	// This requires that spin-down drivers (MagneticDipole) ran first.
	const auto &spin_block = dYdt.Block(State::StateTag::Spin);
	const double Omega_dot = (spin_block.size() > 0) ? spin_block[0] : 0.0;

	if (std::abs(Omega_dot) == 0.0)
		return;

	// -------------------------------------------------------------------------
	//  3. For each species, accumulate d(eta_i)/dt
	// -------------------------------------------------------------------------
	// d(eta_i)/dt += Z_i * 2 * Omega * Omega_dot
	//
	// The mapping from species index to ChemState component depends on
	// the ChemState layout. For a standard npe-mu star:
	//   component 0 = eta_e  (electron channel)
	//   component 1 = eta_mu (muon channel)
	//
	// The Z_i from RotochemicalCache are indexed by species (n, p, e, mu).
	// For the chemical imbalance, we need the *reaction-channel* Z, which
	// combines species contributions:
	//   d(eta_e)/dt  = (Z_n - Z_p - Z_e) * 2 * Omega * Omega_dot
	//   d(eta_mu)/dt = (Z_n - Z_p - Z_mu) * 2 * Omega * Omega_dot
	//
	// This mapping will be finalized when the ChemState layout is defined.
	// For now, we provide the per-species Z_i and leave the reaction-channel
	// combination to be wired at integration time.

	const double factor = 2.0 * Omega * Omega_dot;

	for (std::size_t i = 0; i < cache_.species.size(); i++)
	{
		double deta_dt = cache_.species[i].Z * factor;
		dYdt.AddTo(State::StateTag::Chem, i, deta_dt);
	}
}

} // namespace CompactStar::Physics::Driver::Chem
