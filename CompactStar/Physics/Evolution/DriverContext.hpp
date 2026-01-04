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
 * @file DriverContext.hpp
 * @brief Lightweight read-only context passed to evolution drivers.
 *
 * @ingroup PhysicsEvolution
 *
 * DriverContext is a **non-owning** bundle of pointers to slowly-varying / static
 * objects that drivers may consult when computing RHS contributions.
 *
 * Design goals:
 *  - Keep the driver interface explicit and stable.
 *  - Avoid circular dependencies between drivers and EvolutionSystem.
 *  - Prevent heavy includes in driver-facing headers.
 *
 * Ownership:
 *  - All pointers are **non-owning** and must remain valid for the duration of
 *    the evolution run (typically owned by EvolutionSystem or the caller).
 *
 * Optional components:
 *  - Some pointers (e.g., envelope boundary condition models) may be nullptr.
 *    They are **not** required globally; drivers that require them must validate
 *    their presence based on their own Options/configuration.
 */

#ifndef CompactStar_Physics_Evolution_DriverContext_H
#define CompactStar_Physics_Evolution_DriverContext_H

namespace CompactStar
{
namespace EOS
{
class CompOSE_Thermo; ///< CompOSE thermodynamics table interface (optional).
}
namespace Physics
{

// -----------------------------------------------------------------------------
// Forward declarations (keep this header lightweight).
// -----------------------------------------------------------------------------
namespace Evolution
{
class StarContext;
class GeometryCache;
struct Config;

class StateVector;	  ///< Composite view over sub-states (Spin/Thermal/Chem/…)
class RHSAccumulator; ///< Write-only accumulator for dY/dt components
class StateLayout;	  ///< Layout of the StateVector
} // namespace Evolution

namespace Driver
{
namespace Thermal
{
namespace Boundary
{
class IEnvelope; ///< Tb -> Ts boundary condition interface (optional in context).
} // namespace Boundary
} // namespace Thermal
} // namespace Driver

namespace Evolution
{

/**
 * @struct DriverContext
 * @brief Read-only static context for evolution drivers.
 *
 * DriverContext bundles together pointers to static model objects that drivers
 * may consult while computing dY/dt.
 *
 * It is passed **by const reference** into IDriver::AccumulateRHS().
 *
 * Key principles:
 *  - **Non-owning** pointers
 *  - **Read-only** usage contract for drivers
 *  - **Optional components** allowed (nullptr is valid)
 *
 * Validation policy:
 *  - EvolutionSystem should validate only globally required pointers (e.g. star/geo/cfg
 *    depending on your run mode).
 *  - Feature-specific pointers (e.g. envelope) must be validated by the driver that
 *    uses them, based on driver Options.
 */
struct DriverContext
{
	/**
	 * @brief Stellar structure context (static background star).
	 *
	 * Provides access to structure profiles, metric functions, EOS hooks, etc.
	 *
	 * May be null for drivers that are purely local / toy models, but most physics
	 * drivers will require this.
	 */
	const StarContext *star = nullptr;

	/**
	 * @brief Precomputed geometric cache for spatial integrals and metric factors.
	 *
	 * Provides commonly used quantities such as:
	 *  - proper-volume weights,
	 *  - exp(nu), exp(2 nu), exp(lambda), etc.
	 *
	 * May be null for drivers that do not perform spatial integrals.
	 */
	const GeometryCache *geo = nullptr;

	/**
	 * @brief Optional thermal boundary condition model (Tb -> Ts mapping).
	 *
	 * This is the envelope/blanket interface used to map base-of-envelope temperature
	 * to an effective surface temperature.
	 *
	 * Notes:
	 *  - This pointer is **optional** and may be nullptr.
	 *  - Drivers must only require it when their configured surface model needs it
	 *    (e.g., PhotonCooling with SurfaceModel::EnvelopeTbTs).
	 */
	const Driver::Thermal::Boundary::IEnvelope *envelope = nullptr;

	/**
	 * @brief Optional CompOSE thermodynamics table (for Cv, entropy, etc.).
	 *
	 * Non-owning pointer. May be nullptr if the run does not use thermal microphysics
	 * from CompOSE tables.
	 *
	 * Drivers that require it (e.g. GR heat capacity from CompOSE entropy) must validate
	 * presence based on Options.
	 */
	const CompactStar::EOS::CompOSE_Thermo *thermo = nullptr;

	/**
	 * @brief Global evolution configuration (policy + toggles + numerics).
	 *
	 * Treated as read-only by all drivers.
	 */
	const Config *cfg = nullptr;
};

} // namespace Evolution
} // namespace Physics
} // namespace CompactStar

#endif /* CompactStar_Physics_Evolution_DriverContext_H */